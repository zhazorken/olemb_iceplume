      subroutine obcmom(rv)
c
c    $Id: imag.F,v 1.6 1998/10/21 21:56:38 dwd Exp $
c
c this routine sets the lateral image points
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "cyclic.inc"
#include "olemtke.inc"
#include "moddef.inc"
#include "communicate.inc"
#if defined(ICE)
#include "ice.inc"
#endif
#include "finitv.inc"
c
      real rv(0:idim1,0:jdim1,0:kdim1)
      real vblm,vbrm,zpgr,ypgr
      real vr,vl,cd
      integer i,j,k,olev,ierr
c
      if((blockstart(my_pe+1)+jdim).ge.jtotal) then
        j = jdim1
        do i=1,idim
          do k=0,kdim1
           vbrm = v(i,j,k) + maxv/1.2
c
            j = jdim1
            if(vbrm.lt.0.0) then
#if !defined(ICE)
               v(i,j,k) = vgeo
#endif
c               v(i,j,k) = vmm(i,j,k)
c               u(i,j,k) = umm(i,j,k)
c                u(i,j,k) = ugeo
               w(i,j,k) = wmm(i,j,k)
            else
               if(vbrm.lt.vmm(i,j-1,k)) then
                 vbrm = vmm(i,j-1,k)
               endif
               v(i,j,k)=vmm(i,j,k)
     *             -delt*vbrm*(vmm(i,j,k)-vmm(i,j-1,k))*ody
               u(i,jdim,k) = umm(i,jdim,k) - delt*
     *                   (vm(i,jdim,k)+vm(i-1,jdim,k))*
     *                    (umm(i,jdim,k)-umm(i,jdim-1,k))*ody
               w(i,jdim,k) = wmm(i,jdim,k) - delt*
     *                   (vm(i,jdim,k)+vm(i,jdim,k-1))*
     *                    (wmm(i,jdim,k)-wmm(i,jdim-1,k))*ody

            endif
c 
c smooth out boundary velocities
c
c            u(i,j-1,k) = u(i,j-1,k) + 0.1*(um(i+1,j-1,k)+um(i-1,j-1,k) -
c     *                4.*um(i,j,k) + um(i,j,k+1)+um(i,j,k-1))*delt
c            w(i,j-1,k) = w(i,j-1,k) + 0.1*(wm(i+1,j-1,k)+wm(i-1,j-1,k) -
c     *                4.*wm(i,j-1,k) + wm(i,j-1,k+1)+wm(i,j-1,k-1))*delt
c            v(i,j,k) = v(i,j,k) + 0.1*(vm(i+1,j,k)+vm(i-1,j,k) -
c     *                4.*vm(i,j,k) + vm(i,j,k+1)+vm(i,j,k-1))*delt
c
            u(i,j,k) = u(i,j-1,k)
            w(i,j,k) = w(i,j-1,k)
            if(v(i,j-1,k).gt.0.0) then
              kmom(i,j,k) = kmom(i,j-1,k)
              kheat(i,j,k) = kheat(i,j-1,k)
              p(i,j,k) = p(i,j-1,k)
c              p(i,j,k) = p(i,j-1,k) - dx*rv(i,j-1,k)
c     *                -dx*(v(i,j,k)-vmm(i,j,k))/delt
            else
              p(i,j,k) = p(i,j-1,k)
c              p(i,j,k) = p(i,j-1,k) - dx*rv(i,j-1,k)
c     *                -dx*(v(i,j,k)-vmm(i,j,k))/delt

            endif
c
c wall or zero gradient boundary condition
c
            v(i,j,k) = 0.0
            v(i,j-1,k) = 0.0
            u(i,j,k) = u(i,j-1,k)
            w(i,j,k) = w(i,j-1,k)
            if(k.eq.0.or.k.eq.1) then
              w(i,j,k) = 0.0
            endif
            p(i,j,k) = p(i,j-1,k)
          enddo
        enddo
c   
c compute the average outward velocity on the rhs boundary
c
        vr = 0.0
        do i=1,idim
          do k=1,kdim
            vr = vr + v(i,jdim1,k)
          enddo
        enddo
        vr = vr/float(idim*kdim)
        write(*,*)"vr ",vr
      endif
c #if defined(USEMPI)
c
c need to send it to other processors
c
c      do np=1,nprocs
c         if(blockstart(np)+jdim.ge.jtotal) then
c             sproc = np-1
c         endif
c      enddo
c      call MPI_BCAST(vr,1,mpi_real,sproc,comm,ierr)
c #endif

      call MPI_BARRIER(comm,ierr)
c
#if !defined(RECIRC)
      if(blockstart(my_pe+1).eq.1) then
c 
        vl = 0.0
        do i=1,idim
          do k=0,kdim1
c 
c compute the Doppler shifted phase speeds at the left and right 
c boundaries
c 
         j = 1
c
c set outflow velocity using maximum allowable velocity
c
         vblm = v(i,1,k)-maxv/1.2
         if(vblm.gt.0.0 .and. vmm(i,1,k).lt.0.0) then
           vblm = v(i,1,k)
         endif
         vr = vblm
         if(i.eq.29.and.k.eq.kdim-1) then
           write(*,*) "vblm ",vblm
         endif
c
c         if(vblm.gt.vgeo) then
c           vr = vgeo
c         else
c           vr = vblm
c         endif
c         vr = vgeo
         if(vr.gt.0.0) then
           v(i,1,k) = vr
c           u(i,0,k) = ugeo
           u(i,1,k) = u(i,1,k) - delt*vr*(umm(i,1,k)-u(i,0,k))*odx
           w(i,1,k) = w(i,1,k) - delt*vr*(wmm(i,1,k)-w(i,0,k))*odx
           v(i,0,k) = v(i,1,k)
c            u(i,1,k) = u(i,1,k) + 0.1*(um(i+1,1,k)+um(i-1,1,k) -
c    *                4.*um(i,1,k) + um(i,1,k+1)+um(i,1,k-1))*delt
c            w(i,1,k) = w(i,1,k) + 0.1*(wm(i+1,1,k)+wm(i-1,1,k) -
c    *                4.*wm(i,1,k) + wm(i,1,k+1)+wm(i,1,k-1))*delt

           w(i,0,k) = 0.0
           w(i,1,k) = 0.0
c           u(i,1,k) = umm(i,1,k)
c           p(i,0,k) = p(i,1,k) + dx*rv(i,2,k)
c     *                +dx*(v(i,1,k)-vmm(i,1,k))/delt
         else
c             v(i,1,k) = vr
               
           if(vblm.lt.0.0) then
c              if(k.eq.kdist(i,1)) then
c               olev = kdist(i,1)
c                v(i,1,k)=vmm(i,1,k)
c     *           -delt*(vblm*(vmm(i,2,k)-vmm(i,1,k))*ody
c     *           -randn(i,1)*cd*(vmm(i,1,olev)-vice+vave)*
c     *             abs(vmm(i,1,olev)-vice+vave))
c              else
               v(i,1,k)=vmm(i,1,k)
     *           -delt*vblm*(vmm(i,2,k)-vmm(i,1,k))*ody
c              endif
               u(i,1,k) = umm(i,1,k) - delt*
c     *                 0.5*(vmm(i,1,k)+vmm(i-1,1,k))*
     *                  vblm*
     *                  (umm(i,2,k)-umm(i,1,k))*ody
               w(i,1,k) = wmm(i,1,k) - delt*
c     *                 0.5*(vmm(i,1,k)+vmm(i,1,k-1))*
     *                  vblm*
     *                  (wmm(i,2,k)-wmm(i,1,k))*ody

c           else
c             v(i,1,k) = vmm(i,1,k)
c             u(i,1,k) = u(i,2,k)
c             w(i,1,k) = w(i,2,k)
           endif
           v(i,0,k) = v(i,1,k)
           u(i,0,k) = u(i,1,k)
           w(i,0,k) = w(i,1,k)
c
c do wall bc
c
           v(i,0,k) = 0.0
           w(i,0,1) = 0.0
           v(i,1,k) = 0.0
           w(i,1,k) = 0.0
c smooth out boundary velocities
c
c           u(i,1,k) = u(i,1,k) + 0.1*(um(i+1,1,k)+um(i-1,1,k) -
c     *                4.*um(i,j,k) + um(i,1,k+1)+um(i,1,k-1))*delt
c           w(i,1,k) = w(i,1,k) + 0.1*(wm(i+1,1,k)+wm(i-1,1,k) -
c     *                4.*wm(i,1,k) + wm(i,1,k+1)+wm(i,1,k-1))*delt
c           v(i,1,k) = v(i,1,k) + 0.1*(vm(i+1,1,k)+vm(i-1,1,k) -
c     *                4.*vm(i,1,k) + vm(i,1,k+1)+vm(i,1,k-1))*delt
         endif
         p(i,0,k) = p(i,1,k)
c         u(i,0,k) = u(i,1,k)
c         w(i,0,k) = w(i,1,k)
c         v(i,0,k) = v(i,1,k)

        enddo
      enddo
c
      endif
#endif
c 
c
      return
      end
