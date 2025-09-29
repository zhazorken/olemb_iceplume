      subroutine topbotb(phi,iswitch,tstep)
c
c   $Id: topbotb.F,v 1.4 1998/10/21 21:56:42 dwd Exp $
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
#include "communicate.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real rwtop,rwbot
#if defined(cray)
      real rho
#else
      real*8 rho
#endif
#if defined(RADIATE)
      integer kk
#endif
      integer i,j,k,iswitch,tstep
c
c top and bottom boundary
c I have inserted placemark for future stress b.c. calculation
c
      k=kdim1
      if (iswitch.eq.MDWVEL) then
c
c rigid lid
c     
#if !defined(RADTOP)
         do j=0,jdim1
            do i=0,idim1
               phi(i,j,k) = 0.0
            enddo
         enddo
#endif
c
      else if (iswitch.eq.MDUVEL) then

c
c compute wind shear at model top based on similarity profile
c
c           phi(i,j,kdim1) = phi(i,j,kdim)
c    *           -sign(1.0,ustr)*dz*sqrt(abs(ustr)/rhob(kdim))
c    *           *odz/0.35
c try new boundary condition for u(kdim)
c
         if(strton)then
            do j=0,jdim1
               do i=0,idim1
                  if(str_tbc(i,j).eq.0)then
c
c no slip bc
c
                     phi(i,j,kdim1)=-phi(i,j,kdim)
                  else
c
c stress bc with mask
c
c                     phi(i,j,kdim) = phi(i,j,kdim)+str_tmask(i,j)*
c     *                    delt*odz*ustr/rhob(kdim)
                     phi(i,j,kdim1) = phi(i,j,kdim)
                  endif
               enddo
            enddo
         else
c
c stress bc without mask
c
            do j=0,jdim1
               do i=0,idim1
c                  phi(i,j,kdim) = phi(i,j,kdim)+delt*odz*ustr/rhob(kdim)
#if defined(RADTOP)
               phi(i,j,kdim1) = phi(i,j,kdim)+
     *                        (phi(i,j,kdim)-phi(i,j,kdim-1))
#else
                  phi(i,j,kdim1) = phi(i,j,kdim)
#endif
               enddo
            enddo
         endif
c
      else if (iswitch.eq.MDVVEL) then
c
c           phi(i,j,kdim1) = phi(i,j,kdim)
c    *           -sign(1.0,vstr)*dz*sqrt(abs(vstr)/rhob(kdim))
c    *           *2.*odz/0.35
c
         if(strton)then
            do j=0,jdim1
               do i=0,idim1
                  if(str_tbc(i,j).eq.0)then
c     
c no slip bc
c     
                     phi(i,j,kdim1)=-phi(i,j,kdim)
                  else
c
c stress bc with mask
c
c                     phi(i,j,kdim) = phi(i,j,kdim)+str_tmask(i,j)*
c     *                    delt*odz*vstr/rhob(kdim)
                     phi(i,j,kdim1) = phi(i,j,kdim)
                  endif
               enddo
            enddo
         else
c
c stress bc without mask
c
            do j=0,jdim1
               do i=0,idim1
c                  phi(i,j,kdim) = phi(i,j,kdim)+delt*odz*vstr/rhob(kdim)
#if defined(RADTOP)
               phi(i,j,kdim1) = phi(i,j,kdim)+
     *                        (phi(i,j,kdim)-phi(i,j,kdim-1))
#else
                  phi(i,j,kdim1) = phi(i,j,kdim)
#endif
               enddo
            enddo
         endif
c     
#if defined(RADTOP)
      else if (iswitch.eq.MDTEMP) then
c
c
         do j=0,jdim1
           do i=0,idim1
             phi(i,j,k) = phi(i,j,kdim) + dtdztop
           enddo
         enddo
      else if (iswitch.eq.MDSAL) then
c
c
         do j=0,jdim1
           do i=0,idim1
             phi(i,j,k) = phi(i,j,kdim) + dsdztop
           enddo
         enddo
c
#else
c
      else if (tstep.lt.ncool.and.iswitch.eq.MDTEMP) then
c
c heat flux at the surface
c
c            phi(i,j,k)=phi(i,j,kdim)+dz*heatf*randn(i,j)/(cp*rho0
c     *           *kheat(i,j,kdim))
c
c temperature change in the top model level
c
         if(hfton)then
c
c heat flux bc with a mask
c if ice model, then heat flux is controled by icemod.F
c
            do j=0,jdim1
               do i=0,idim1
#if !defined(ICE) && !defined(HBOT)
                  phi(i,j,kdim)=phi(i,j,kdim)+
     *                 hf_tmask(i,j)*delt*heatf*randn(i,j)/(cp*
     *                 rho0*dz)
#endif
                  phi(i,j,k) = phi(i,j,kdim)
               enddo
            enddo
         else
c
c heat flux bc without a mask
c if ice model then heat flux is handled by icemod.F
c
#if defined(ICE) || defined(HBOT)
            do j=0,jdim1
              do i=0,idim1
                phi(i,j,k) = phi(i,j,kdim)
              enddo
            enddo
#else
c
c not ice model
c
#if !defined(RADIATE)
            do j=0,jdim1
               do i=0,idim1
                  phi(i,j,kdim)=phi(i,j,kdim)+
     *                 delt*heatf*randn(i,j)/(cp*rho0*dz)
                  phi(i,j,k) = phi(i,j,kdim)
               enddo
            enddo
#else
c            if(my_pe.eq.1) then
c              write(*,*)"heatf,tstep,delt,mod,sw ",heatf,tstep,delt,
c     *         mod(tstep,int(60./delt)),swheatf
c            endif
            if(mod(tstep,int(5./delt)).eq.0) then
              if(swheatf.gt.0.0) then
                 do j=0,jdim1
                    do i=0,idim1

                       do kk=1,kdim
                          phi(i,j,kk)=phi(i,j,kk)+
     *                         5.*swheatf*rad(kk)*randn(i,j)/(cp*
     *                         rho0*dz)
                       enddo
                    enddo
                 enddo
              endif
              do j=0,jdim1
                 do i=0,idim1
                    phi(i,j,kdim)=phi(i,j,kdim)+
     *                   5.*heatf*randn(i,j)/(cp*rho0*dz)
                    phi(i,j,k) = phi(i,j,kdim)
                 enddo
              enddo
            endif
#endif
c ice model endif
#endif
         endif
c
      else if (tstep.lt.ncool.and.iswitch.eq.MDSAL
     *     .and.tstep.gt.1)then
c
c salt flux at the surface (parts per thousand)
c salinity change in the top model level
c
c put in extra if statement to account for ice coverage
c assume zero change in salinity from evaporation if ice covered
c
c
         if(lfton)then
c
c latent heat flux with a mask
c
            do j=0,jdim1
               do i=0,idim1
#if !defined(ICE)
                  rho=1000.0d0+rhop(i,j,kdim)+rhobm1(kdim)
#if defined(cray)
                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
     *                 (rho+lf_tmask(i,j)*delt*latent/(lv*dz))
#else 
c                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
c     *                 (rho+lf_tmask(i,j)*dble(delt*latent*randn(i,j)
c     *                 /(lv*dz)))
                  phi(i,j,kdim) = phi(i,j,kdim) -
     *             lf_tmask(i,j)*phi(i,j,kdim)*latent*delt/(1000.*lv*dz)
#endif
#endif
                  phi(i,j,k)=phi(i,j,kdim)
               enddo
            enddo
         else
c
c latent heat flux without a mask
c
            do j=0,jdim1
               do i=0,idim1
#if !defined(ICE)
                  rho=1000.0d0+rhop(i,j,kdim)+rhobm1(kdim)
#if defined(cray)
c                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
c     *                 (rho+delt*latent/(lv*dz))
#else 
c                  phi(i,j,kdim)=phi(i,j,kdim)*rho/
c     *                 (rho+dble(delt*latent/(lv*dz)))
                  phi(i,j,kdim) = phi(i,j,kdim) -
     *             phi(i,j,kdim)*latent*delt/(1000.*lv*dz)
#endif
#endif
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
c
c RADTOP endif
c
#endif
      else if (iswitch.eq.MDWVEL) then
c
c rigid lid
c
#if !defined(RADTOP)
         do j=0,jdim1
            do i=0,idim1
               phi(i,j,k) = 0.0
            enddo
         enddo
#endif

      else if (iswitch.eq.MDPRES) then
#if !defined(RADTOP)
c         do j=0,jdim1
         do j=0,jdim
            do i=1,idim
c
c zero pressure gradient
c
c     
c #if defined(STOKES)
c                rwtop = -0.5*(um(i,j,kdim)+um(i+1,j,kdim))*tomegay -
c     *           us(kdim1)*tomegay - 
c     *              uave*tomegay -
c     *              us(kdim1)*0.5*(um(i+1,j,kdim)+um(i,j,kdim) -
c     *              um(i+1,j,kdim-1)-um(i,j,kdim-1))*odz
c     *              -vs(kdim1)*0.5*(vm(i,j+1,kdim)+vm(i,j,kdim)-
c     *              vm(i,j+1,kdim-1)-vm(i,j,kdim-1))*odz +
c     *              grav*rhop(i,j,kdim)/rho0 -
c     *              2.*kmom(i,j,kdim)*wmm(i,j,kdim)*odz*odz
c#else
c               rwtop = -0.5*(um(i,j,kdim)+um(i+1,j,kdim))*tomegay -
c     *              uave*tomegay +
c     *              grav*rhop(i,j,kdim)/rho0 -
c     *              2.*kmom(i,j,kdim)*wmm(i,j,kdim)*odz*odz
c#endif
cc
c               phi(i,j,kdim1) = -dz*rwtop
c               phi(i,j,kdim1) = phi(i,j,kdim)-dz*rwtop
                phi(i,j,kdim1) = phi(i,j,kdim)
            enddo
         enddo
#endif
      else if(iswitch.eq.MDRW) then
#if !defined(RADTOP)
c
c rw at the top.  w=0 is assumed, du/dz can be non-zero at top.
c
c         do j=0,jdim1
         do j=0,jdim
            do i=1,idim
c     
#if defined(STOKES)
               phi(i,j,kdim1) = 
     *              -0.5*(um(i,j,kdim)+um(i+1,j,kdim))*tomegay -
     *              us(kdim1)*tomegay - 
     *              us(kdim1)*0.5*(um(i+1,j,kdim)+um(i,j,kdim) -
     *              um(i+1,j,kdim-1)-um(i,j,kdim-1))*odz
     *              -vs(kdim1)*0.5*(vm(i,j+1,kdim)+vm(i,j,kdim)-
     *              vm(i,j+1,kdim-1)-vm(i,j,kdim-1))*odz +
     *              grav*rhop(i,j,kdim)/rho0 -
     *              2.*kmom(i,j,kdim)*wmm(i,j,kdim)*odz*odz
#else
               phi(i,j,kdim1) =
     *              -0.5*(um(i,j,kdim)+um(i+1,j,kdim))*tomegay +
     *              grav*rhop(i,j,kdim)/rho0 -
     *              2.*kmom(i,j,kdim)*wmm(i,j,kdim)*odz*odz
#endif
            enddo
         enddo
#endif
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
      if (iswitch.eq.MDWVEL) then
         do j=0,jdim1
            do i=0,idim1
#if defined(RADIATE)
c
c assume constant gradient for vertical velocity at bottom
c
c               phi(i,j,0) = phi(i,j,1)-(phi(i,j,2)-phi(i,j,1))
               phi(i,j,0) = phi(i,j,1)
c do nothing for now
#else
c
c flat bottom
c
               phi(i,j,k)=0.0
               k=1
               phi(i,j,k)=0.0
#endif
            enddo
         enddo
      else if(iswitch.eq.MDUVEL
     *        .or. iswitch.eq.MDVVEL) then
         do j=0,jdim1
            do i=1,idim
#if defined(RADIATE)
c
c constant gradient of u,v at bottom for radiation conditions
c
c               phi(i,j,0) = phi(i,j,1)-(phi(i,j,2)-phi(i,j,1))
                phi(i,j,0) = phi(i,j,1)
c
c
#else
c
#if defined(CONSGRAD)
c 
               phi(i,j,0) = 2.*phi(i,j,1) - phi(i,j,2) 
#else
c
c no-slip at bottom
c
c set to constant gradient for now
c
c               phi(i,j,0) = -phi(i,j,1)
               phi(i,j,0) = 2.*phi(i,j,1) - phi(i,j,2) 
#endif
#endif
            enddo
         enddo
c
c if radiation conditions, don't mess with the lower boundary pressure
c
      else if (iswitch.eq.MDPRES) then
#if !defined(RADIATE)
c
c zero pressure gradient
c
         do j=0,jdim1
            do i=1,idim 
c               rwbot =  -0.5*(um(i,j,1)+um(i+1,j,1))*tomegay + 
c     *              grav*rhop(i,j,1)/rho0 -
c     *              2.*kmom(i,j,1)*wmm(i,j,2)*odz*odz
c
c               phi(i,j,0)=phi(i,j,1)+dz*rwbot
                phi(i,j,0) = phi(i,j,1)
            enddo
         enddo
#endif
#if defined(RADIATE)
      else if(iswitch.eq.MDTEMP)then
c 
c for radiation conditions,
c advect background temperature and salinity 
c 
         do j=0,jdim1
            do i=0,idim1
c     phi(i,j,0) = phi(i,j,1)-(phi(i,j,2)-phi(i,j,1))
               phi(i,j,0) = phi(i,j,1)-dtdzlow
            enddo
         enddo
      else if(iswitch.eq.MDSAL) then
        do j=0,jdim1
          do i=0,idim1
            phi(i,j,0) = phi(i,j,1)-dsdzlow
          enddo
        enddo
#endif
      else if(iswitch.eq.MDRW)then
c
c if radiation conditions, leave rw alone
c
#if !defined(RADIATE)

c
c rw at the bottom.  w=0 is assumed, du/dz is zero at bottom.
c
         do j=1,jdim
            do i=1,idim
               phi(i,j,1) = -0.5*(um(i,j,1)+um(i+1,j,1))*tomegay + 
     *              grav*rhop(i,j,1)/rho0 -
     *              2.*kmom(i,j,1)*wmm(i,j,2)*odz*odz
            enddo
         enddo
#endif
      else
c
c zero gradient (no flux)
c MDNULL uses this branch
c
         do j=0,jdim1
            do i=0,idim1
               phi(i,j,k)=phi(i,j,1)
            enddo
         enddo
      endif
c
c all done setting image points
c
      return
      end
c
      subroutine topbott(phi,iswitch,tstep)
c
c This subroutine computes the image points for temporary
c u, v and w
c Only needed for smola velocities
c
c This routine is never called as far as I can tell.
c
c
c
#include "param.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k,iswitch,tstep
c
c top and bottom boundary
c I have inserted placemark for future stress b.c. calculation
c
      j=1
c
      do i=0,idim1
         k=kdim1
         if (iswitch.eq.MDWVEL) then
c
c rigid lid
c
            phi(i,j,k) = 0.0
         else
c
c zero gradient (no flux of the quantity)
c MDNULL uses this branch
c
            phi(i,j,k) = phi(i,j,kdim)
         end if 
c
c bottom boundary condition is zero gradient or no flux
c placemark for stress b.c.
c
         k=0
         if (iswitch.eq.MDWVEL) then
c
c flat bottom
c
            phi(i,j,k)=0.0
            k=1
            phi(i,j,k)=0.0
         else
c
c zero gradient (no flux)
c MDNULL uses this branch
c
            phi(i,j,k)=phi(i,j,1)
         endif
      enddo
c
c all done setting image points
c
      return
      end

