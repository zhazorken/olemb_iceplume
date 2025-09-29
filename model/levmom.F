c
c     =================================================
      subroutine levmom(n,q,qout,ivar)
c     =================================================
c
c   $Id: leveq77.F,v 1.7 1998/10/16 22:12:35 dwd Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "ice.inc"
#include "communicate.inc"
#include "finitv.inc"
c
      integer msize,mitot,mjtot,mktot
#if defined(USEMPI)
      integer ier
      parameter ( msize = idim+4 )
#else
      parameter ( msize = jdim+4 )
#endif
      parameter ( mitot = idim+4,mjtot = jdim+4, mktot = kdim+4)
c
      real q(0:idim1,0:jdim1,0:kdim1)
      real qout(0:idim1,0:jdim1,0:kdim1)
      real u1d(msize),v1d(msize),w1d(msize)
      real f(mitot,mjtot,mktot),g(mitot,mjtot,mktot)
      real h(mitot,mjtot,mktot)
      real delta(mitot,mjtot,mktot),deladd(msize),q1d(msize)
      real fadd(msize),gadd(msize,2),hadd(msize,3,2)
      real kright,kleft,kback,ktop,kbot,sgs,b1,b2,b3,kfront
#if defined(ANALY)
      real*8 pel(kdim),zz
      real den(0:idim1,0:jdim1,0:kdim1)
#endif
#if defined(USEMPI)
      real qnorth(0:idim1,0:kdim1),qsouth(0:idim1,0:kdim1)
      real unorth(0:idim1,0:kdim1),usouth(0:idim1,0:kdim1)
      real vnorth(0:idim1,0:kdim1),vsouth(0:idim1,0:kdim1)
      real wnorth(0:idim1,0:kdim1),wsouth(0:idim1,0:kdim1)
#endif
      integer method,mthlim,olev
      integer i,j,k,ivar,n
      common /mth/ method,mthlim
c
c using routine with MC limiter
c
      method = 4
      mthlim = 2
c
c
c
c     # initialize fluxes:
c
      do k=1,mktot
       do j = 1, mjtot
        do i = 1, mitot
          f(i,j,k) = 0.0
          g(i,j,k) = 0.0
          h(i,j,k) = 0.0
          delta(i,j,k) = 0.0
        enddo
       enddo
      enddo
c
c
c     # perform x-sweeps
c     ==================
c
      do 50 k=2,mktot-1
       do 50 j=2,mjtot-1
        do i=2,mitot-1
            q1d(i) = q(i-2,j-2,k-2)
            u1d(i) = u(i-2,j-2,k-2) 
            v1d(i) = v(i-2,j-2,k-2)
            w1d(i) = w(i-2,j-2,k-2)
        enddo
c
c do min boundary points
c
        q1d(1) = q(idim-1,j-2,k-2)
        u1d(1) = u(idim-1,j-2,k-2)
        v1d(1) = v(idim-1,j-2,k-2)
        w1d(1) = w(idim-1,j-2,k-2)
c
c do max boundary points
c
        q1d(mitot) = q(2,j-2,k-2)
        u1d(mitot) = u(2,j-2,k-2)
        v1d(mitot) = v(2,j-2,k-2)
        w1d(mitot) = w(2,j-2,k-2)
        
c
        call fluxes(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
     &              mitot)
c
        do 30 i=3,mitot-1
            delta(i,j,k) = delta(i,j,k) + deladd(i)
            f(i,j,k) = f(i,j,k) + fadd(i)
            g(i,j,k) = g(i,j,k) + gadd(i,1)
            g(i,j+1,k) = g(i,j+1,k) + gadd(i,2)
c
            h(i,j-1,k) = h(i,j-1,k) + hadd(i,1,1)
            h(i,j,k) = h(i,j,k) + hadd(i,2,1)
            h(i,j+1,k) = h(i,j+1,k) + hadd(i,3,1)
c
            h(i,j-1,k+1) = h(i,j-1,k+1) + hadd(i,1,2)
            h(i,j,k+1) = h(i,j,k+1) + hadd(i,2,2)
            h(i,j+1,k+1) = h(i,j+1,k+1) + hadd(i,3,2)
   30   continue
   50 continue
c
c     # perform y sweeps
c     ==================
c
#if defined(USEMPI)
c      write(*,*)"Starting south shifts ",my_pe
      call southshift(q,qsouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift(v,vsouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift(u,usouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift(w,wsouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
c
c      write(*,*)"Starting north shifts ",my_pe
      call northshift(q,qnorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call northshift(u,unorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call northshift(v,vnorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call northshift(w,wnorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
c     write(*,*)"Performed all shifts in ",my_pe
c
#endif
      do 100 k = 2,mktot-1
       do 100 i = 2,mitot-1
        do j=2,mjtot-1
            q1d(j) = q(i-2,j-2,k-2)
            u1d(j) = v(i-2,j-2,k-2)
            v1d(j) = w(i-2,j-2,k-2)
            w1d(j) = u(i-2,j-2,k-2)
        enddo
c
c if using MPI, then we need to get second layer of 
c boundary points from north and south neighbors
c
#if defined(USEMPI)
c
c do min boundary points
c
        q1d(1) = qsouth(i-2,k-2)
        u1d(1) = vsouth(i-2,k-2)
        v1d(1) = wsouth(i-2,k-2)
        w1d(1) = usouth(i-2,k-2)
c
c do max boundary points
c
        q1d(mjtot) = qnorth(i-2,k-2)
        u1d(mjtot) = vnorth(i-2,k-2)
        v1d(mjtot) = wnorth(i-2,k-2)
        w1d(mjtot) = unorth(i-2,k-2)
#else
#if defined(OPENBC)
        q1d(1) = q(i-2,0,k-2)
        u1d(1) = v(i-2,0,k-2)
        v1d(1) = w(i-2,0,k-2)
        w1d(1) = u(i-2,0,k-2)
c
c do max boundary points
c
        q1d(mjtot) = q(i-2,jdim1,k-2)
        u1d(mjtot) = v(i-2,jdim1,k-2)
        v1d(mjtot) = w(i-2,jdim1,k-2)
        w1d(mjtot) = u(i-2,jdim1,k-2)
#else
        q1d(1) = q(i-2,jdim-1,k-2)
        u1d(1) = v(i-2,jdim-1,k-2)
        v1d(1) = w(i-2,jdim-1,k-2)
        w1d(1) = u(i-2,jdim-1,k-2)
c
c do max boundary points
c
        q1d(mjtot) = q(i-2,2,k-2)
        u1d(mjtot) = v(i-2,2,k-2)
        v1d(mjtot) = w(i-2,2,k-2)
        w1d(mjtot) = u(i-2,2,k-2)
#endif
#endif
c
c        # note: u, v, and w are rotated for solving in the y-direction
c        # The roles of fadd, gadd, and hadd are also rotated.
c
         call fluxes(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
     &               mjtot)
c
         do 80 j=3,mjtot-1
           delta(i,j,k) = delta(i,j,k) + deladd(j)
           g(i,j,k) = g(i,j,k) + fadd(j)
           h(i,j,k) = h(i,j,k) + gadd(j,1)
           h(i,j,k+1) = h(i,j,k+1) + gadd(j,2)
c
           f(i,j,k-1) = f(i,j,k-1) + hadd(j,1,1)
           f(i,j,k) = f(i,j,k) + hadd(j,2,1)
           f(i,j,k+1) = f(i,j,k+1) + hadd(j,3,1)
c
           f(i+1,j,k-1) = f(i+1,j,k-1) + hadd(j,1,2)
           f(i+1,j,k) = f(i+1,j,k) + hadd(j,2,2)
           f(i+1,j,k+1) = f(i+1,j,k+1) + hadd(j,3,2)
   80    continue
  100 continue
c
c
c     # perform z sweeps
c     ==================
c
      do 140 j = 2,mjtot-1
       do 140 i = 2,mitot-1
        do k=2,mktot-1
            q1d(k) = q(i-2,j-2,k-2)
            u1d(k) = w(i-2,j-2,k-2)
            v1d(k) = u(i-2,j-2,k-2)
            w1d(k) = v(i-2,j-2,k-2)
        enddo
c
c do min boundary points
c
        if(ivar.eq.MDTEMP) then
           q1d(1) = q(i-2,j-2,0)-dtdzlow
	else if(ivar.eq.MDSAL) then
           q1d(1) = q(i-2,j-2,0)-dsdzlow
        else
           q1d(1) = q(i-2,j-2,1)
        endif
        u1d(1) = w(i-2,j-2,1)
        v1d(1) = u(i-2,j-2,1)
        w1d(1) = v(i-2,j-2,1)
c
c do max boundary points
c
        q1d(mktot) = q(i-2,j-2,kdim1) 
        u1d(mktot) = w(i-2,j-2,kdim1)
        v1d(mktot) = u(i-2,j-2,kdim1)
        w1d(mktot) = v(i-2,j-2,kdim1)
c
c        # note: u, v, and w are rotated for solving in the z-direction
c        # The roles of fadd, gadd, and hadd are also rotated.
c
         call fluxes(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
     &               mktot)
c
         do 120 k=3,mktot-1
           delta(i,j,k) = delta(i,j,k) + deladd(k)
           h(i,j,k) = h(i,j,k) + fadd(k)
           f(i,j,k) = f(i,j,k) + gadd(k,1)
           f(i+1,j,k) = f(i+1,j,k) + gadd(k,2)
c
           g(i-1,j,k) = g(i-1,j,k) + hadd(k,1,1)
           g(i,j,k) = g(i,j,k) + hadd(k,2,1)
           g(i+1,j,k) = g(i+1,j,k) + hadd(k,3,1)
c
           g(i-1,j+1,k) = g(i-1,j+1,k) + hadd(k,1,2)
           g(i,j+1,k) = g(i,j+1,k) + hadd(k,2,2)
           g(i+1,j+1,k) = g(i+1,j+1,k) + hadd(k,3,2)
  120    continue
  140 continue
c
c
      do i=1,idim
        do j=1,jdim
          h(i,j,mktot-1) = 0.0
        enddo
      enddo
c
c     # update q values by differencing fluxes
c     ========================================
c

      do  k = 3,mktot-2
       do  j = 3,mjtot-2
        do  i = 3,mitot-2
#if defined(ICE) || defined(HBOT)
           qout(i-2,j-2,k-2) =
     &                     odx * (ax(i-1,j-2,k-2)*f(i+1,j,k) - 
     &                            ax(i-2,j-2,k-2)*f(i,j,k)
     &                    + ay(i-2,j-1,k-2)*g(i,j+1,k) - 
     &                            ay(i-2,j-2,k-2)*g(i,j,k)
     &                    + h(i,j,k+1) - h(i,j,k) + delta(i,j,k))
#else
          qout(i-2,j-2,k-2) =  
     &                     odx * (f(i+1,j,k) - f(i,j,k)
     &                    + g(i,j+1,k) - g(i,j,k) 
     &                    + h(i,j,k+1) - h(i,j,k) + delta(i,j,k))
#endif
        enddo
       enddo
      enddo
c finish up the processing by mixing through turbulence parameterization.
c Don't do this from the TKE fields.  TKE is diffused in the sgs update
c routine.
c
c      call image(qout)
      call topbotb(qout,MDNULL,n)
c
c
c      call imag(qout)
c
c  **  do compressiblity correction
c
c      do j=1,jdim
c      do i=1,idim
c         do k=1,kdim
c            qout(i,j,k)=qout(i,j,k)+((u(i+1,j,k)-u(i,j,k))*dtdx+
c     *              (w(i,j,k+1)-w(i,j,k))*dtdz
c     *              + (v(i,j+1,k)-v(i,j,k))*dtdy)*q(i,j,k)
c         enddo
c      enddo
c      enddo
c
c
c
c  **  get out of here normally
c
c
      return
      end
c
c
