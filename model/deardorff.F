      subroutine deardorff(n,trbe,trbep)
c
c   $Id: deardorff.F,v 1.3 1998/07/24 16:59:31 dwd Exp $
c
c     This routine computes the shear production, buoyant production 
c diffusion, and dissipation of the tke parameter, q.  q is actually 
c the square-root of the turbulent kinetic energy divided by 1/2.
c Then, new values for the turbulent diffusion are calculated.
c
c This routine uses a different notation for dimensions than the
c straight 2-D OCM.
c
c i is left-right, j is front-back, and k is bottom-top
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "moddef.inc"
#include "cyclic.inc"
#include "dentable.inc"
c
      real trbe(0:idim1,0:jdim1,0:kdim1),trbep(0:idim1,0:jdim1,0:kdim1)
      real dudz,dudx,dwdz,dwdx
      real dudy,dwdy,dvdy,dvdx,dvdz,difdy
      real shear,disp,diftke,ftke,coef
      real r1,r2,drhodz,brunts,difdz,difdx,brunt,lambda
      real dels,s
#if defined(cray)
      real testb
#endif

      integer i,j,k,n
c
      real t,tmpt,tmps,y1,y2,y3,y4,ti,si
      real const1,const2,const3,const4,const5,ntests,ntest
#if defined(TKESHEAR)
      real uav(0:kdim1),vav(0:kdim1)
      real dvavdz(0:kdim1),duavdz(0:kdim1)
      integer jj
#endif
      integer kz,index,ii,kk
#if defined(TKESHEAR)
      save uav,vav,dvavdz,duavdz
#endif
c
c empirical constants 
c
      real Cdg,Cel,Ces,Cl,Ckm,Chl,Chs
c
      Cdg = 2.0
      Cel = 0.19
      Ces = 0.51
      Cl  = 0.76
      Ckm = 0.1
      Chl = 1.0
      Chs = 2.0
c
c Setup constants
c
c ds=dx=dy=dz  (not valid if dx, dy, and dz are not equal)
c
      dels = ds
c
      const1 = 0.5*Ces/dels
      const2 = 0.5*Cel
      const3 = 0.5*(Cel+Ces)/dels
      const4 = Chl+Chs
      const5 = 0.5*Ckm
c
c scaled test for lambda evaluation
c
c ntest = Cl*trbe/dels
c
      ntest = 1e-8*Cl/dels
      ntests = ntest*ntest
c
#if defined(TKESHEAR)
      do k=0,kdim
         uav(k) = 0.0
         vav(k) = 0.0
      enddo
      do k=0,kdim
         do jj=1,jdim
            do i=1,idim
               uav(k) = um(i,jj,k)+uav(k)
               vav(k) = vm(i,jj,k)+vav(k)
            enddo
         enddo
      enddo
      do k=1,kdim-1
         duavdz(k) = (uav(k+1)-uav(k-1))*tdz/float(idim*jdim)
         dvavdz(k) = (vav(k+1)-vav(k-1))*tdz/float(idim*jdim)
      enddo
#endif
c     
      do j=1,jdim
         do k=1,kdim-1
            do i=1,idim
         
c 
c First, compute the shear production terms
c
c
c Because of staggered grid, some of the points are averaged
c fdz,fdx,fdy are 1/4*dz, ...
c
               dudz = (um(i,j,k+1)+um(i+1,j,k+1)-um(i,j,k-1)
     *              -um(i+1,j,k-1))*fdz
#if defined(TKESHEAR)
     *              -duavdz(k)
#endif
               dvdz = (vm(i,j,k+1)+vm(i,j+1,k+1)-vm(i,j,k-1)
     *              -vm(i,j+1,k-1))*fdz
#if defined(TKESHEAR)
     *              -dvavdz(k)
#endif
               dwdx = (wm(i+1,j,k)+wm(i+1,j,k+1)-wm(i-1,j,k)
     *              -wm(i-1,j,k+1))*fdx
               dudy = (um(i,j+1,k)+um(i+1,j+1,k)-um(i,j-1,k)
     *              -um(i+1,j-1,k))*fdy
               dvdx = (vm(i+1,j,k)+vm(i+1,j+1,k)-vm(i-1,j,k)
     *              -vm(i-1,j+1,k))*fdx
               dwdy = (wm(i,j+1,k)+wm(i,j+1,k+1)-wm(i,j-1,k)
     *              -wm(i,j-1,k+1))*fdy

               dudx = (um(i+1,j,k)-um(i,j,k))*odx
               dwdz = (wm(i,j,k+1)-wm(i,j,k))*odz
               dvdy = (vm(i,j+1,k)-vm(i,j,k))*ody

               shear = 2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz) +
     *              (dudy+dvdx)**2 + (dudz+dwdx)**2 + (dvdz+dwdy)**2
c
c Move on to the buoyancy term
c
c Compute the vertical density gradient
c We must treat top special because th image point is distorted to get
c correct heat flux into top th 
c
               s=sal(i,j,k-1)
               t=th(i,j,k-1)
               kz=2*k
c
               tmpt=(t-tmin)/dtemp
               tmps=(s-smin)/dsal
c
               ii=ifix(tmpt)+1
               kk=ifix(tmps)+1
c
               index=ii+(kk-1)*tsize+(kz-1)*tssize
               y1=rhotbl(index)
               y2=rhotbl(index+1)
               y3=rhotbl(index+1+tsize)
               y4=rhotbl(index+tsize)
               ti=tmpt-(ii-1)
               si=tmps-(kk-1)
               r2=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
c     
               s=sal(i,j,k+1)
               t=th(i,j,k+1)
               kz=2*k
c     
               tmpt=(t-tmin)/dtemp
               tmps=(s-smin)/dsal
c     
               ii=ifix(tmpt)+1
               kk=ifix(tmps)+1
c     
               index=ii+(kk-1)*tsize+(kz-1)*tssize
               y1=rhotbl(index)
               y2=rhotbl(index+1)
               y3=rhotbl(index+1+tsize)
               y4=rhotbl(index+tsize)
               ti=tmpt-(ii-1)
               si=tmps-(kk-1)
               r1=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
c
               drhodz = (r1-r2)*tdz
c     
               brunts = -grav*drhodz/rhob(k)
c
c Next, compute diffusion of trbe
c tdzs = 1/(dz*dz*2)
c
               difdz = ((kmom(i,j,k+1)+kmom(i,j,k))*
     *              (trbe(i,j,k+1)-trbe(i,j,k)) -
     *              (kmom(i,j,k)+kmom(i,j,k-1))*
     *              (trbe(i,j,k)-trbe(i,j,k-1)))*tdzs
               difdx = ((kmom(i+1,j,k)+kmom(i,j,k))*
     *              (trbe(i+1,j,k)-trbe(i,j,k)) -
     *              (kmom(i,j,k)+kmom(i-1,j,k))*
     *              (trbe(i,j,k)-trbe(i-1,j,k)))*tdxs
               difdy = ((kmom(i,j+1,k)+kmom(i,j,k))*
     *              (trbe(i,j+1,k)-trbe(i,j,k)) -
     *              (kmom(i,j,k)+kmom(i,j-1,k))*
     *              (trbe(i,j,k)-trbe(i,j-1,k)))*tdys
c
               diftke = Cdg*(difdz+difdx+difdy)
c
c Finally, compute the dissipation term and update tke term,q,
c then update diffusion coefficients, kmom and kheat
c
c Check to see if stable, decrease lambda if stable, otherwise
c lambda is dels
c
c dels = (dx*dy*dz) to 1/3 power
c
#if defined(cray)
               testb=brunts-ntests
               lambda=(0.5+sign(0.5,testb))*Cl*trbe(i,j,k)/
     *             sqrt(amax1(ntests,brunts))+
     *              (0.5+sign(0.5,-testb))*dels
               lambda=amin1(dels,lambda)
               lambda=amax1(dels*1e-6,lambda)
               disp = (const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
               coef = Chl+Chs*lambda/dels
               ftke=0.5*kmom(i,j,k)*(-shear + 
     *            coef*brunts)/amax1(1e-8,trbe(i,j,k)) - diftke + disp
               trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
               kmomp(i,j,k) = Ckm*lambda*trbep(i,j,k)+kmb
               kheat(i,j,k) = coef*kmomp(i,j,k)+khb
#else
               if(drhodz.lt.0.0) then
                  brunt = sqrt(brunts)
                  if(brunt.gt.ntest) then
                     lambda = Cl*trbe(i,j,k)/brunt
                  else
                     lambda = dels
                  endif
                  lambda=amin1(dels,lambda)
                  lambda=amax1(dels*1e-6,lambda)
                  disp = (const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
                  coef = Chl+Chs*lambda/dels
                  ftke=0.5*kmom(i,j,k)*(-shear + 
     *                 coef*brunts)/amax1(1e-8,trbe(i,j,k)) - 
     *                 diftke + disp
                  trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
                  kmomp(i,j,k) = Ckm*lambda*trbep(i,j,k)
                  kheat(i,j,k) = coef*kmomp(i,j,k)
               else
                  disp = const3*trbe(i,j,k)*trbe(i,j,k)
                  ftke=0.5*kmom(i,j,k)*(-shear + 
     *                 const4*brunts)/amax1(1e-8,trbe(i,j,k)) - 
     *                 diftke + disp
                  trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
                  kmomp(i,j,k) = Ckm*dels*trbep(i,j,k)
                  kheat(i,j,k) = const4*kmomp(i,j,k)
               endif
               kmomp(i,j,k)=kmomp(i,j,k)+kmb
               kheat(i,j,k)=kheat(i,j,k)+khb
#endif
c
            enddo
         enddo
c
c
         k=kdim
         do i=1,idim
c 
c First, compute the shear production terms
c
c
c Because of staggered grid, some of the points are averaged
c fdz,fdx,fdy are 1/4*dz, ...
c
c upper boundary stress condition which uses similarity
c profile to calculate stress at point kdim-0.5*delx
c
            dudz = 5.0*sqrt(abs(ustr)/rhob(kdim))*odz
            dvdz = 5.0*sqrt(abs(vstr)/rhob(kdim))*odz
            dwdx = (wm(i+1,j,k)+wm(i+1,j,k+1)-wm(i-1,j,k)
     *           -wm(i-1,j,k+1))*fdx
            dudy = (um(i,j+1,k)+um(i+1,j+1,k)-um(i,j-1,k)
     *           -um(i+1,j-1,k))*fdy
            dvdx = (vm(i+1,j,k)+vm(i+1,j+1,k)-vm(i-1,j,k)
     *           -vm(i-1,j+1,k))*fdx
            dwdy = (wm(i,j+1,k)+wm(i,j+1,k+1)-wm(i,j-1,k)
     *           -wm(i,j-1,k+1))*fdy
         
            dudx = (um(i+1,j,k)-um(i,j,k))*odx
            dwdz = (wm(i,j,k+1)-wm(i,j,k))*odz
            dvdy = (vm(i,j+1,k)-vm(i,j,k))*ody

            shear = 2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz) +
     *           (dudy+dvdx)**2 + (dudz+dwdx)**2 + (dvdz+dwdy)**2
c
c Move on to the buoyancy term
c
c Compute the vertical density gradient
c We must treat top special because th image point is distorted to get
c correct heat flux into top th 
c
            s=sal(i,j,k-1)
            t=th(i,j,k-1)
            kz=2*k
c
            tmpt=(t-tmin)/dtemp
            tmps=(s-smin)/dsal
c     
            ii=ifix(tmpt)+1
            kk=ifix(tmps)+1
c
            index=ii+(kk-1)*tsize+(kz-1)*tssize
            y1=rhotbl(index)
            y2=rhotbl(index+1)
            y3=rhotbl(index+1+tsize)
            y4=rhotbl(index+tsize)
            ti=tmpt-(ii-1)
            si=tmps-(kk-1)
            r2=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
c     
            s=sal(i,j,k)
            t=th(i,j,k)
            kz=2*k
c     
            tmpt=(t-tmin)/dtemp
            tmps=(s-smin)/dsal
c     
            ii=ifix(tmpt)+1
            kk=ifix(tmps)+1
c     
            index=ii+(kk-1)*tsize+(kz-1)*tssize
            y1=rhotbl(index)
            y2=rhotbl(index+1)
            y3=rhotbl(index+1+tsize)
            y4=rhotbl(index+tsize)
            ti=tmpt-(ii-1)
            si=tmps-(kk-1)
            r1=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
c     
            drhodz = (r1-r2)/dz
c     
            brunts = -grav*drhodz/rhob(k)
c
c Next, compute diffusion of trbe
c tdzs = 1/(dz*dz*2)
c
            difdz = ((kmom(i,j,k+1)+kmom(i,j,k))*
     *           (trbe(i,j,k+1)-trbe(i,j,k)) -
     *           (kmom(i,j,k)+kmom(i,j,k-1))*
     *           (trbe(i,j,k)-trbe(i,j,k-1)))*tdzs
            difdx = ((kmom(i+1,j,k)+kmom(i,j,k))*
     *           (trbe(i+1,j,k)-trbe(i,j,k)) -
     *           (kmom(i,j,k)+kmom(i-1,j,k))*
     *           (trbe(i,j,k)-trbe(i-1,j,k)))*tdxs
            difdy = ((kmom(i,j+1,k)+kmom(i,j,k))*
     *           (trbe(i,j+1,k)-trbe(i,j,k)) -
     *           (kmom(i,j,k)+kmom(i,j-1,k))*
     *           (trbe(i,j,k)-trbe(i,j-1,k)))*tdys
c     
            diftke = Cdg*(difdz+difdx+difdy)
c
c Finally, compute the dissipation term and update tke term,q,
c then update diffusion coefficients, kmom and kheat
c
c Check to see if stable, decrease lambda if stable, otherwise
c lambda is dels
c
c dels = (dx*dy*dz) to 1/3 power
c
c disp = 3.9*epsilon/2q = 3.9*0.5*C*q^2/lambda
c
#if defined(cray)
            testb=brunts-ntests
            lambda=(0.5+sign(0.5,testb))*Cl*trbe(i,j,k)/
     *           sqrt(amax1(ntests,brunts))+
     *           (0.5+sign(0.5,-testb))*dels
            lambda=amin1(dels,lambda)
            lambda=amax1(dels*1e-6,lambda)
#if defined(RADIATE)
            disp=3.9*(const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
#else
            disp=3.9*(const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
#endif
            coef = Chl+Chs*lambda/dels
            ftke=0.5*kmom(i,j,k)*(-shear + 
     *           coef*brunts)/amax1(1e-8,trbe(i,j,k)) - diftke + disp
            trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
            kmomp(i,j,k) = Ckm*lambda*trbep(i,j,k)+kmb
            kheat(i,j,k) = coef*kmomp(i,j,k)+khb
#else
            if(drhodz.lt.0.0) then
               brunt = sqrt(brunts)
               if(brunt.gt.ntest) then
                  lambda = Cl*trbe(i,j,k)/brunt
               else
                  lambda = dels
               endif
               lambda=amin1(dels,lambda)
               lambda=amax1(dels*1e-6,lambda)
#if defined(RADIATE)
               disp=3.9*(const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
#else
               disp=3.9*(const2/lambda+const1)*trbe(i,j,k)*trbe(i,j,k)
#endif
               coef = Chl+Chs*lambda/dels
               ftke=0.5*kmom(i,j,k)*(-shear + 
     *              coef*brunts)/amax1(1e-8,trbe(i,j,k)) - diftke + disp
               trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
               kmomp(i,j,k) = Ckm*lambda*trbep(i,j,k)
               kheat(i,j,k) = coef*kmomp(i,j,k)
            else
#if defined(RADIATE)
               disp = 3.9*const3*trbe(i,j,k)*trbe(i,j,k)
#else
               disp = 3.9*const3*trbe(i,j,k)*trbe(i,j,k)
#endif
               ftke=0.5*kmom(i,j,k)*(-shear + 
     *            const4*brunts)/amax1(1e-8,trbe(i,j,k)) - diftke + disp
               trbep(i,j,k) = dim(trbe(i,j,k),ftke*delt)
               kmomp(i,j,k) = Ckm*dels*trbep(i,j,k)
               kheat(i,j,k) = const4*kmomp(i,j,k)
            endif
            kmomp(i,j,k)=kmomp(i,j,k)+kmb
            kheat(i,j,k)=kheat(i,j,k)+khb
#endif
         enddo
      enddo
c
c Alldone with interiors, next do boundaries
c
      call image(trbep)
      call image(kmomp)
      call image(kheat)
      call topbotb(trbep,MDTKE,n)
      call topbotb(kmomp,MDKM,n)
      call topbotb(kheat,MDKH,n)
c     
      return
      end
