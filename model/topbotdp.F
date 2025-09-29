#if !defined(SNGLSALT)
      subroutine topbotdp(phi,iswitch,tstep)
c
c   $Id: topbotdp.F,v 1.4 1998/10/16 21:21:26 dwd Exp $
c
c This subroutine computes the image points for scalars with a
c heat flux imposed at the model top
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
c
      real*8 phi(0:idim1,0:jdim1,0:kdim1)
      real*8 rho
      integer i,j,k,iswitch,tstep
c
      if(iswitch.ne.MDSAL.and.
     *     iswitch.ne.MDNULL)stop 'topbotdp not called for MDSAL'
c
c top and bottom boundary
c
      k=kdim1
      if (tstep.lt.ncool.and.iswitch.eq.MDSAL
     *     .and.tstep.gt.1)then
c
c salt flux at the surface (parts per thousand)
c salinity change in the top model level
c
         if(hfton)then
c
c latent heat flux bc with a mask
c
            do j=0,jdim1
               do i=0,idim1
                  rho=1000.0d0+rhop(i,j,kdim)+rhobm1(kdim)
                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
     *                 (rho+lf_tmask(i,j)*dble(delt*latent*randn(i,j)
     *                 /(lv*dz)))
                  phi(i,j,k)=phi(i,j,kdim)
               enddo
            enddo
         else
c
c latent heat flux bc without a mask
c
            do j=0,jdim1
               do i=0,idim1
                  rho=1000.0d0+rhop(i,j,kdim)+rhobm1(kdim)
                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
     *                 (rho+dble(delt*latent/(lv*dz)))
                  phi(i,j,k)=phi(i,j,kdim)
               enddo
            enddo
         endif
c
c next do the rainfall
c
         if(rainon) then
            do j=0,jdim1
               do i=0,idim1
                  phi(i,j,kdim) = phi(i,j,kdim)-
     *                 rainf*delt*odz*phi(i,j,kdim)
                  phi(i,j,k) = phi(i,j,kdim)
               enddo
            enddo
         endif
      else
c
c zero gradient (no flux of the quantity)
c MDNULL uses this branch
c
         do j=0,jdim1
            do i=0,idim1
               phi(i,j,k) = phi(i,j,kdim)
            enddo
         enddo
      end if 
c
c bottom boundary condition is zero gradient or no flux
c placemark for stress b.c.
c
      k=0
c
c zero gradient (no flux)
c MDNULL uses this branch
c
      do j=0,jdim1
         do i=0,idim1
            phi(i,j,k)=phi(i,j,1)
         enddo
      enddo
c
c all done setting image points
c
      return
      end
#else
      subroutine topbotdp(phi,iswitch,tstep)
      integer iswitch,tstep
      real*8 phi
      end
#endif

