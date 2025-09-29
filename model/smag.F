      subroutine smag(n)
c
c
c
c     $Id: smag.F,v 1.2 1998/07/24 16:59:54 dwd Exp $
c
c     This routine computes the eddy momentum and heat viscosity terms
c     using the method of Smagorinsky (1963) as described by Lilly (1962).
c
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
c
      real dudz,dudx,dwdz,dwdx
      real dudy,dwdy,dvdy,dvdx,dvdz
      real shear,ri,sheart
      real r1,r2,drhodz,brunts
      real iesfast,dels,s,kmm

      integer i,k,j,n
c
c Setup constants
c
c ds=dx=dy=dz  (not valid if dx, dy, and dz are not equal)
c
c      dels = (dx*dx*dzv(30))**.3333

c 
c
c Do loop through interior of grid first
c
c      dels = sqrt(dx*dy)
c      kmb = 0.001*dels*dels/delt
c      khb = kmb
      do 99  j=1,jdim
        do 100 i=1,idim
         do 101 k=1,kdim
         
            dels = dz
            kmm =0.1*dels*dels/delt
c            kmm =0.8*dzv(k)*dzv(k)/delt
c 
c First, compute the shear production terms
c
c
c Because of staggered grid, some of the points are averaged
c fdz,fdx,fdy are 1/4*dz, ...
c
            dudz = 0.5*(um(i,j,k+1)+um(i+1,j,k+1)-um(i,j,k-1)
     *           -um(i+1,j,k-1))/(dz+dz)
            dwdx = (wm(i+1,j,k)+wm(i+1,j,k+1)-wm(i-1,j,k)
     *           -wm(i-1,j,k+1))*fdx
            dudy = (um(i,j+1,k)+um(i+1,j+1,k)-um(i,j-1,k)
     *           -um(i+1,j-1,k))*fdy
            dvdx = (vm(i+1,j,k)+vm(i+1,j+1,k)-vm(i-1,j,k)
     *           -vm(i-1,j+1,k))*fdx
            dwdy = (wm(i,j+1,k)+wm(i,j+1,k+1)-wm(i,j-1,k)
     *           -wm(i,j-1,k+1))*fdy
            dvdz = 0.5*(vm(i,j,k+1)+vm(i,j+1,k+1)-vm(i,j,k-1)
     *           -vm(i,j+1,k-1))/(dz+dzs)

            dudx = (um(i+1,j,k)-um(i,j,k))*odx
            dwdz = (wm(i,j,k+1)-wm(i,j,k))/dz
            dvdy = (vm(i,j+1,k)-vm(i,j,k))*ody
c
            shear = 2.*(dudx**2+dvdy**2+dwdz**2)+(dudy+dvdx)**2 +
     *             (dudz+dwdx)**2 + (dvdz+dwdy)**2
c
c
c Move on to the buoyancy term
c
c Compute the vertical density gradient
c We must treat top special because th image point is distorted to get
c correct heat flux into top th 
c
c            s=sal(i,j,k-1)
c            r2 = iesfast(s,th(i,j,k-1),2*k)
c            if(k.ne.kdim) then
c               s=sal(i,j,k+1)
c               r1 = iesfast(s,th(i,j,k+1),2*k)
c               drhodz = (r1-r2)*tdz
c            else
c               s=sal(i,j,k)
c               r1 = iesfast(s,th(i,j,k),2*k)
c               drhodz = (r1-r2)/dz
c            endif
c
c        r1 = rhob(k+1)+rhop(i,j,k+1)
c        r2 = rhob(k-1)+rhop(i,j,k-1)
            brunts = grav*(th(i,j,k+1)-th(i,j,k-1))/
     *                (th(i,j,k)*(dz+dz))
c
c
            if(abs(shear).gt.0.01) then
               ri=brunts/shear
            else
               ri = 10.0
            endif
            if(ri.lt.0.0)then
               sheart = sqrt(shear)
               if(ri.gt.0.0) then
                 kmom(i,j,k) = 0.028*dels*dels*sheart*sqrt(1.-4.*ri) 
     *               + kmb
               else
                 kmom(i,j,k)= 0.028*dels*dels*sheart*sqrt(1.-ri) + kmb
               endif
               kheat(i,j,k) = kmom(i,j,k)+khb
            else
               kmom(i,j,k) = kmb
               kheat(i,j,k) = khb
            endif
            kmom(i,j,k) = min(kmom(i,j,k),kmm)
            if(kmom(i,j,k).gt.10.*kmm)then
               write(*,*)"i,j,k,sheart,kmb,ri,km ",i,j,k,sheart,kmb,ri,
     *          kmom(i,j,k)
            endif
 101     continue
         kmm =0.5*dz*dz/delt
         kmom(i,j,kdim) = uistar(i,j)*0.2*dz
         kmom(i,j,kdim) = min(kmom(i,j,kdim),kmm)
         kheat(i,j,kdim) = kmom(i,j,kdim)
c         if(kmom(i,j,1).gt.20000.0) then
c           write(*,*)"kmom,uistar ",kmom(i,j,1),uistar(i,j)
c         endif
 100   continue
  99  continue
c
c Alldone with interiors, next do boundaries
c
      call imag(kmom)
      call topbotb(kmom,MDKM,n)
      call imag(kheat)
      call topbotb(kheat,MDKH,n)
c
      return
      end
