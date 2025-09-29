      subroutine obcmom(ru,n)
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
      real udis
      real vr,vl,cd,rr,kout,usum,usuml,masd
      integer i,j,k,olev,ierr,ii,jj,iw,n

c
c
c
c ramp up discharge
c
c#if defined(DISCHARGE)
c      if(n*delt.lt.1200.) then
c         disch = 50.*n*delt/1200.
c      else
c         disch = 50.
c      endif
c#else
c
c      disch = 0.0
c#endif
c      disch = 600.

      lh = 3.34e5/cp
#if defined(SPLUME)
      kout = 0.
      do k=1,kdim
        do j=1,jdim
             jj = j+jblkstart(my_pe+1)
             ii = k
             rr = sqrt(real(((kdim/2-ii))**2+((jtotal/2-jj))**2))
             if(rr.le.5) then
                kout = kout+1
             endif
        enddo
      enddo
#endif

      usuml = 0.0
      if((iblkstart(my_pe+1)+idim).ge.itotal) then
        i = idim1
        kdist(i,j) = kdist(idim,j)
        do j=0,jdim1
          do k=0,kdim1
c
c wall or zero gradient boundary condition
c
c
c            u(i,j,k) =   0.1*9.*10./(jtotal*kdim)
c            um(i,j,k) =  0.1*9.*10./(jtotal*kdim)
c            umm(i,j,k) = 0.1*9.*10./(jtotal*kdim)
c
#if defined(SPLUME)
            u(i,j,k) =   kout*1.5/(jtotal*kdim)
            um(i,j,k) =  kout*1.5/(jtotal*kdim)
            umm(i,j,k) = kout*1.5/(jtotal*kdim)
#else
c            if(k.gt.kdim) then
            if(k.gt.kdim-60) then
              u(i,j,k) =  disch/(dy*dy*jtotal*60)
              um(i,j,k) = disch/(dy*dy*jtotal*60)
              umm(i,j,k) =  disch/(dy*dy*jtotal*60)
            else
              u(i,j,k) = 0.0
              um(i,j,k) = 0.0
              umm(i,j,k) = 0.0
            endif

#endif
c            if(u(idim,j,k) .gt.0.0) then
c               u(i,j,k) = u(idim,j,k)
c               um(i,j,k) = u(idim,j,k)
c               umm(i,j,k) = u(idim,j,k)
c               usuml = usuml+u(i,j,k)
c            endif
        

c            u(i,j,k) = 0.0
c            um(i,j,k) = 0.0
c            umm(i,j,k) = 0.0
c            u(i-1,j,k) = u(i,j,k)
c            um(i-1,j,k) = um(i,j,k)
c            umm(i-1,j,k) = umm(i,j,k)
            if(k.gt.kdim-60) then
            sal(i,j,k) = sal(i-1,j,k)
            th(i,j,k) = th(i-1,j,k)
            v(i,j,k) = v(i-1,j,k)
            vm(i,j,k) = v(i-1,j,k)
            vmm(i,j,k) = v(i-1,j,k)
            w(i,j,k) = w(i-1,j,k)
            wm(i,j,k) = w(i-1,j,k)
            wmm(i,j,k) = w(i-1,j,k)
            ax(i,j,k) = 1.0
            bax(i,j,k) = 1.0
            ax(idim,j,k) = 1.0
            bax(idim,j,k) = 1.0
c            ax(0,j,k) = 1.0
            if(k.eq.0.or.k.eq.1.or.k.gt.kdist(i,j)) then
              w(i,j,k) = 0.0
              wm(i,j,k) = 0.0
              wmm(i,j,k) = 0.0
            endif
            p(i,j,k) = p(i-1,j,k)
            ru(i,j,k) = ru(i-1,j,k)
            endif

          enddo
        enddo
c        call MPI_REDUCE(usuml,usum,1,
c     &        MPI_REAL,MPI_SUM,mastertask,
c     &        comm, ierr)
c        masd = (usum*dxx*dxx-300.)
c        do j=1,jdim
c           do k=1,kdim
c              u(i,j,k) = u(i,j,k)-masd/(jtotal*kdim)
c              um(i,j,k) = u(i,j,k)
c              umm(i,j,k) = u(i,j,k)
c           enddo
c        enddo



      endif
c
c sum up outflow boundary velocities over all processors
c easier than trying to limit sum to just the processors along the edge of the domain
c
c      call MPI_REDUCE(usuml,usum,1,
c     &        MPI_REAL,MPI_SUM,mastertask,
c     &        comm, ierr)
cc
c      if((iblkstart(my_pe+1)+idim).ge.itotal) then
c        i = idim1
c        masd = (usum*dxx*dxx-300.)
c        do j=1,jdim
c           do k=1,kdim
c              u(i,j,k) = u(i,j,k)-masd/(jtotal*kdim)
c              um(i,j,k) = u(i,j,k)
c              umm(i,j,k) = u(i,j,k)
c           enddo
c        enddo
c      endif

c   

      call MPI_BARRIER(comm,ierr)
      
c
      cd = 0.16/(log(1.0/0.02))**2
      if(iblkstart(my_pe+1).eq.1) then
c 
        do j=0,jdim1
#if defined(HBOT)
          kdist(0,j) = kdist(1,j)
          do k=0,kdist(19,j)
#else
          do k=0,kdim
#endif
            jj = jblkstart(my_pe+1)+j
c 
c 
         i = 1
c
c do wall bc
c
           v(0,j,k) = v(1,j,k)
           vm(0,j,k) = vm(1,j,k)
           vmm(0,j,k) = vmm(1,j,k)
c           if(k.lt.0) then
c           if(k.lt.11.and.jj.gt.107.and.jj.lt.117) then
c            u(0,j,k) = 0.1
c            um(0,j,k) = 0.1
c            umm(0,j,k) = 0.1
c            u(1,j,k) = 0.1
c            um(1,j,k) = 0.1
c            umm(1,j,k) = 0.1
#if defined(SPLUME)
             jj = j+jblkstart(my_pe+1)
             ii = k
             rr = sqrt(real(((kdim/2-ii))**2+((jtotal/2-jj))**2))
             if(rr.le.5) then
               u(0,j,k) = 0.5
               um(0,j,k) = 0.5
               umm(0,j,k) = 0.5
               u(1,j,k) = 0.5
               um(1,j,k) = 0.5
               umm(1,j,k) = 0.5
               th(0,j,k) = th(1,j,k)
               sal(0,j,k) = sal(1,j,k)

#else
c
c 100 x 40  = 400 units
c
c           if(k.lt.11.and.jj.gt.780.and.jj.lt.811) then
c           if(k.lt.5.and.jj.gt.745.and.jj.lt.846) then
c           if(k.lt.5.and.jj.gt.940.and.jj.lt.1041) then
c           if(k.lt.5.and.jj.gt.670.and.jj.lt.771) then
           if(k.lt.5.and.jj.gt.430.and.jj.lt.531) then
            udis = disch/(dy*dy*100.*4.)
c            udis = 0.0
c            u(0,j,k) = udis
c            um(0,j,k) = udis
c            umm(0,j,k) = udis
c            u(1,j,k) = udis
c            um(1,j,k) = udis
c            umm(1,j,k) = udis
c            if(jj.gt.395.and.jj.lt.596) then
c              u(0,j,k) = 0.125
c              um(0,j,k) = 0.125
c              umm(0,j,k) = 0.125
c              u(1,j,k) = 0.125
c              um(1,j,k) = 0.125
c              umm(1,j,k) = 0.125
c            endif
c            umm(1,j,k) = 0.375

            sal(0,j,k) = 0.5
            th(0,j,k) = 0.1

#endif
            u(0,j,k) = 0.0
            um(0,j,k) = 0.0
            umm(0,j,k) = 0.0
            u(1,j,k) = 0.0
            um(1,j,k) = 0.0
            umm(1,j,k) = 0.0
cc
c            ax(0,j,k) = 1.0
c            bax(0,j,k) = 1.0
c            ax(1,j,k) = 1.0
c            bax(1,j,k) = 1.0
c            th(0,j,k) = th(1,j,k)
c            sal(0,j,k) = sal(1,j,k)
           else
            u(1,j,k) = 0.0
            um(1,j,k) = 0.0
            umm(1,j,k) = 0.0
            u(0,j,k) = 0.0
            um(0,j,k) = 0.0
            umm(0,j,k) = 0.0
            th(0,j,k) = th(1,j,k)
            sal(0,j,k) = sal(1,j,k)
            ax(1,j,k) = 0.0
            ax(0,j,k) = 0.0
#if !defined(UNDERCUT)
c            iw = 19
c            ustarice = sqrt(cd*v(iw,j,k)*v(iw,j,k)
c     *          +cd*w(iw,j,k)*w(iw,j,k)
c     *          +cd*0.05*0.05)
c            wice =(th(iw,j,k)+0.0054*sal(iw,j,k))*ustarice*0.0056/lh
c            wice = wice-wice*ran2(icrseed)/50.
cc
c            tke(iw,j,k) = wice
c            dcor(iw,j,k) = ustarice
c            dcor(iw+1,j,k) = ustarice
cc
cc melt flux of 4 m/day
cc
cc            wice = 4/86400.;
c            heatflux = -1.*wice*lh/dxx
c            saltflux = -1.*wice*sal(iw,j,k)/dxx
c            
#if !defined(SPLUME)
            if(k.gt.18) then
c            sal(iw,j,k) = sal(iw,j,k)+saltflux*delt*0.5
c            th(iw,j,k) = th(iw,j,k)+heatflux*delt*0.5
            endif
c            sal(0,j,k) = sal(i,1,k)
c            th(0,j,k) = th(i,1,k)
#endif
#endif
           endif
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
