      subroutine obcmom(ru)
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
      real ru(0:idim1,0:jdim1,0:kdim1)
      real vblm,vbrm,zpgr,ypgr
      real vr,vl,cd,kout,outp
      integer i,j,k,olev,ierr
c
      lh = 3.34e5/cp
      kout = 0.0
      do j=1,jdim
        do i=1,idim
           jj = j+jblkstart(my_pe+1)
           ii = i+iblkstart(my_pe+1)
           rr = sqrt(real(((itotal/2-ii))**2+((jtotal/2-jj))**2))
           if(rr.le.10) then
             kout = kout+1.
           endif
        enddo
      enddo
c      if(iope) write(*,*)"compute nozzle size ",kout
      call MPI_ALLREDUCE(kout,outp,1,MPI_REAL,
     *       MPI_SUM,comm,mpierr)
c      if(iope) write(*,*)"plume size ",outp

      if((iblkstart(my_pe+1)+idim).ge.itotal) then
        i = idim1
        kdist(i,j) = kdist(idim,j)
        do j=0,jdim1
          do k=0,kdim1
c
c wall or zero gradient boundary condition
c
c
c            u(i,j,k) =   0.1*3.*18./(jtotal*kdim)
c            u(i,j,k) =   0.3*8.*37./(jtotal*kdim)
            u(i,j,k) =   1.5*outp*0.5/(jtotal*kdim)
            um(i,j,k) =  1.5*outp*0.5/(jtotal*kdim)
            umm(i,j,k) = 1.5*outp*0.5/(jtotal*kdim)
            u(i-1,j,k) = u(i,j,k)
            um(i-1,j,k) = um(i,j,k)
            umm(i-1,j,k) = umm(i,j,k)
            sal(i,j,k) = sal(i-1,j,k)
            th(i,j,k) = th(i-1,j,k)
            v(i,j,k) = v(i-1,j,k)
            vm(i,j,k) = v(i-1,j,k)
            vmm(i,j,k) = v(i-1,j,k)
            w(i,j,k) = w(i-1,j,k)
            wm(i,j,k) = w(i-1,j,k)
            wmm(i,j,k) = w(i-1,j,k)
            if(k.eq.0.or.k.eq.1.or.k.gt.kdist(i,j)) then
              w(i,j,k) = 0.0
              wm(i,j,k) = 0.0
              wmm(i,j,k) = 0.0
            endif
            p(i,j,k) = p(i-1,j,k)
            ru(i,j,k) = ru(i-1,j,k)
          enddo
        enddo
      endif
c   

      call MPI_BARRIER(comm,ierr)
      
c
      cd = 0.16/(log(1.0/0.01))**2
      if(iblkstart(my_pe+1).eq.1) then
c 
        do j=0,jdim1
#if defined(HBOT)
          kdist(0,j) = kdist(1,j)
          do k=0,kdist(1,j)
#else
          do k=0,kdim
#endif
            jj = jblkstart(my_pe+1)+j
c 
c 
         i = 0
c
c do wall bc
c
           v(0,j,k) = v(1,j,k)
           vm(0,j,k) = vm(1,j,k)
           vmm(0,j,k) = vmm(1,j,k)
            u(i,j,k) =   -1.5*outp*0.5/(jtotal*kdim)
            um(i,j,k) =  -1.5*outp*0.5/(jtotal*kdim)
            umm(i,j,k) = -1.5*outp*0.5/(jtotal*kdim)
c            u(0,j,k) = 0.0
c            um(0,j,k) = 0.0
c            umm(0,j,k) = 0.0
            u(1,j,k) = u(i,j,k)
            um(1,j,k) = um(i,j,k)
            umm(1,j,k) = umm(i,j,k)
            u(2,j,k) = u(i,j,k)
            um(2,j,k) = um(i,j,k)
            umm(2,j,k) = umm(i,j,k)
            sal(i,j,k) = sal(i+1,j,k)
            th(i,j,k) = th(i+1,j,k)
            v(i,j,k) = v(i+1,j,k)
            vm(i,j,k) = v(i+1,j,k)
            vmm(i,j,k) = v(i+1,j,k)
            w(i,j,k) = w(i+1,j,k)
            wm(i,j,k) = w(i+1,j,k)
            wmm(i,j,k) = w(i+1,j,k)
            if(k.eq.0.or.k.eq.1.or.k.gt.kdist(i,j)) then
              w(i,j,k) = 0.0
              wm(i,j,k) = 0.0
              wmm(i,j,k) = 0.0
            endif
            p(i,j,k) = p(i+1,j,k)
            ru(i,j,k) = ru(i+1,j,k)


c           if(k.lt.9.and.jj.gt.232.and.jj.lt.280) then
cc            u(0,j,k) = 0.3
c            u(0,j,k) = 0.0
c            um(0,j,k) = 0.0
c            umm(0,j,k) = 0.0
c            u(1,j,k) = 0.0
c            um(1,j,k) = 0.0
c            umm(1,j,k) = 0.0
c            ax(0,j,k) = 1.0
c            bax(0,j,k) = 1.0
c            ax(1,j,k) = 1.0
c            bax(1,j,k) = 1.0
cc            sal(0,j,k) = 2.0
c            sal(0,j,k) = sal(i,1,k)
cc            th(0,j,k) = 0.0
c            th(0,j,k) = th(1,j,k)
c           else
c            u(1,j,k) = 0.0
c            um(1,j,k) = 0.0
c            umm(1,j,k) = 0.0
c            u(0,j,k) = 0.0
c            um(0,j,k) = 0.0
c            umm(0,j,k) = 0.0
c            th(0,j,k) = th(1,j,k)
c            sal(0,j,k) = sal(1,j,k)
c            ax(1,j,k) = 0.0
c            ax(0,j,k) = 0.0
c            ustarice = sqrt(cd*u(i,1,k)*u(i,1,k)+cd*w(i,1,k)*w(i,1,k))
c            wice =(th(i,1,k)+0.0054*sal(i,1,k))*ustarice*0.0056/lh
c            heatflux = -1.*wice*lh/dx
c            saltflux = -1.*wice*sal(i,1,k)/dx
c            
c            sal(i,1,k) = sal(i,1,k)+saltflux*delt
c            th(i,1,k) = th(i,1,k)+heatflux*delt
c            sal(0,j,k) = sal(i,1,k)
c            th(0,j,k) = th(i,1,k)
c           endif
c           v(i,1,k) = v(0,j,k)
           w(0,j,k) = w(1,j,k)
           wm(0,j,k) = wm(1,j,k)
           wmm(0,j,k) = wmm(1,j,k)
           p(0,j,k) = p(1,j,k)
           ru(0,j,k) = ru(1,j,k)
c           if(k.gt.kdist(0,j).or.k.eq.1) then
c             w(0,j,k) = 0.0
c             wmm(0,j,k) = 0.0
c             wmm(0,j,k) = 0.0
c          endif  
c
c insert a notch
c
c           if(i.gt.150.and.j.lt.jdim-10) then
c             u(i,j,k) = 0.0
c             v(i,j,k) = 0.0
c             w(i,j,k) = 0.0
c             um(i,j,k) = 0.0
c             vm(i,j,k) = 0.0
c             wm(i,j,k) = 0.0
c             umm(i,j,k) = 0.0
c             vmm(i,j,k) = 0.0
c             wmm(i,j,k) = 0.0
c           endif
           w(0,j,1) = 0.0
           wm(0,j,1) = 0.0
           wmm(0,j,1) = 0.0
           w(0,j,kdim1) = 0.0
           wm(0,j,kdim1) = 0.0
           wmm(0,j,kdim1) = 0.0
          enddo
        enddo
c
      endif
      call MPI_BARRIER(comm,ierr)
c 
c
      return
      end
      subroutine imag2p(phi)
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
      real phi(0:idim1,0:jdim1,0:kdim1)

      
      if((iblkstart(my_pe+1)+idim).ge.itotal) then
        i = idim1
        do j=0,jdim1
          do k=0,kdim1
c
c wall or zero gradient boundary condition
c
            phi(i,j,k) = phi(i-1,j,k)
          enddo
        enddo
      endif
c
      call MPI_BARRIER(comm,ierr)

c
      if(iblkstart(my_pe+1).eq.1) then
c 
        do j=0,jdim1
          do k=0,kdim1
c
c do wall bc
c
           phi(0,j,k) = phi(1,j,k)

          enddo
        enddo
c
      endif
      call MPI_BARRIER(comm,ierr)
      return
      end
