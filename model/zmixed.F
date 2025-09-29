      subroutine zmixed(zi)
      real zi
c
c   $Id: zmixed.F,v 1.4 1998/10/21 21:56:44 dwd Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "communicate.inc"
#include "dentable.inc"
c
c
      real s,t
c
      integer kz,index,kk
      real y1,y2,y3,y4,ti,si
      real tmpt,tmps
c
      real tave(kdim),save(kdim),sigmat(kdim)
      real area,a,b,z,diff
#if !defined(USEMPI)
      integer j
#if defined(cray)
      real davet,daves
#else
      real*8 davet,daves
#endif
#endif
      integer i,k
c
      area = float(idim*jdim)
      if(zimethod.eq.0) then
c
c abs(delta theta) < zicriteria
c
c horizontally average the temperature field
c
#if defined(USEMPI)
         call horzavg(th,tave)
#else
         do k=1,kdim
            davet=0.0D0
            do j=1,jdim
               do i=1,idim
                  davet=davet+th(i,j,k)
               enddo
            enddo
            tave(k)=davet/area
         enddo
#endif
c
c search for |t0-t| > zicriteria
c
         if(iope) then
            do k=kdim-1,1,-1
               diff=tave(k)-tave(kdim)
               if(abs(diff).gt.zicriteria) goto 31
            enddo
c
c model depth
c
            zi=dz*kdim
            return
 31         z=dz*(kdim-k+0.5)
            a=(tave(k)-tave(k+1))/dz
            b=tave(k)-a*z
            zi=(tave(kdim)+sign(zicriteria,diff)-b)/a
            return
         endif
      else if(zimethod.eq.1) then
c
c delta sigmat < zicriteria
c
c horizontally average the temperature and salinity field
c
#if defined(USEMPI)
         call horzavg(th,tave)
#if defined(SNGLSALT)
         call horzavg(sal,save)
#else
         call horzavgdp(sal,save)
#endif
#else
         do k=1,kdim
            davet=0.0D0
            daves=0.0D0
            do j=1,jdim
               do i=1,idim
                  davet=davet+th(i,j,k)
                  daves=daves+sal(i,j,k)
               enddo
            enddo
            tave(k)=davet/area
            save(k)=daves/area
         enddo
#endif
c
c compute sigma-t
c
         if(iope) then
            do kk=1,kdim
               s=save(kk)
               t=tave(kk)
               kz=2*kdim
c     
               tmpt=(t-tmin)/dtemp
               tmps=(s-smin)/dsal
c
               i=ifix(tmpt)+1
               k=ifix(tmps)+1
c
               index=i+(k-1)*tsize+(kz-1)*tssize
               y1=rhotbl(index)
               y2=rhotbl(index+1)
               y3=rhotbl(index+1+tsize)
               y4=rhotbl(index+tsize)
               ti=tmpt-(i-1)
               si=tmps-(k-1)
               sigmat(kk)=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
            enddo
c
c search for sigmat - sigmat0 > zicriteria
c
            do k=kdim-1,1,-1
               diff=sigmat(k)-sigmat(kdim)
               if(diff.gt.zicriteria) goto 40
            enddo
c
c model depth
c
            zi=dz*kdim
            return
 40         z=dz*(kdim-k+0.5)
            a=(sigmat(k)-sigmat(k+1))/dz
            b=sigmat(k)-a*z
            zi=(sigmat(kdim)+sign(zicriteria,diff)-b)/a
            return
         endif
      endif
      return
      end



