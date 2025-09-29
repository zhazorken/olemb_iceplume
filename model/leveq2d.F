c
c     =================================================
      subroutine leveq2d(n,q,qout,ivar)
c     =================================================
c
c  $Id: leveq2d.F,v 1.6 1999/08/16 20:22:31 eds Exp $
c
c
c this is a 2-d version of the leveq77 routine
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "ice.inc"
#include "moddef.inc"
#include "communicate.inc"
c
      integer msize,mitot,mjtot
#if defined(USEMPI)
      integer ier
#endif
      parameter ( msize = idim+4 )
      parameter ( mitot = idim+4,mjtot = jdim+4)
c
      real q(0:idim1,0:jdim1)
      real qout(0:idim1,0:jdim1)
      real uu(0:idim1,0:jdim1)
      real vv(0:idim1,0:jdim1)
      real u1d(msize),v1d(msize)
      real f(mitot,mjtot),g(mitot,mjtot)
      real h(mitot,mjtot)
      real delta(mitot,mjtot),deladd(msize),q1d(msize)
      real fadd(msize),gadd(msize,2)
#if defined(USEMPI)
      real qnorth(0:idim1),qsouth(0:idim1)
      real unorth(0:idim1),usouth(0:idim1)
      real vnorth(0:idim1),vsouth(0:idim1)
#endif
      integer method,mthlim
      integer i,j,ivar,n
      common /mth/ method,mthlim
c
c using routine with MC limiter
c
      method = 4
      mthlim = 4
c
c
c
c     # initialize fluxes:
c
      do j = 1, mjtot
        do i = 1, mitot
          f(i,j) = 0.0
          g(i,j) = 0.0
          h(i,j) = 0.0
          delta(i,j) = 0.0
        enddo
      enddo
c
c set velocities to the ice flow velocity for now
c
      do j=0,jdim1
        do i=0,idim1
          uu(i,j) = uice
          vv(i,j) = 0.0
        enddo
      enddo
c
      call imag2d(uu)
      call imag2d(vv)
c
c
c     # perform x-sweeps
c     ==================
c
      do 50 j=2,mjtot-1
        do i=2,mitot-1
            q1d(i) = q(i-2,j-2)
            u1d(i) = uu(i-2,j-2) 
            v1d(i) = vv(i-2,j-2)
        enddo
c
c do min boundary points
c
        q1d(1) = q(idim-1,j-2)
        u1d(1) = uu(idim-1,j-2)
        v1d(1) = vv(idim-1,j-2)
c
c do max boundary points
c
        q1d(mitot) = q(2,j-2)
        u1d(mitot) = uu(2,j-2)
        v1d(mitot) = vv(2,j-2)
        
c
        call flux2d(q1d,u1d,v1d,deladd,fadd,gadd,1,mitot)
c
        do 30 i=3,mitot-1
            delta(i,j) = delta(i,j) + deladd(i)
            f(i,j) = f(i,j) + fadd(i)
            g(i,j) = g(i,j) + gadd(i,1)
            g(i,j+1) = g(i,j+1) + gadd(i,2)
   30   continue
   50 continue
c
c     # perform y sweeps
c     ==================
c
#if defined(USEMPI)
c      write(*,*)"Starting south shifts ",my_pe
      call southshift1(q,qsouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift1(vv,vsouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift1(uu,usouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
c
      call northshift1(q,qnorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call northshift1(uu,unorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call northshift1(vv,vnorth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
c
#endif
      do 100 i = 2,mitot-1
        do j=2,mjtot-1
            q1d(j) = q(i-2,j-2)
            u1d(j) = vv(i-2,j-2)
            v1d(j) = uu(i-2,j-2)
        enddo
c
c if using MPI, then we need to get second layer of 
c boundary points from north and south neighbors
c
#if defined(USEMPI)
c
c do min boundary points
c
        q1d(1) = qsouth(i-2)
        u1d(1) = vsouth(i-2)
        v1d(1) = usouth(i-2)
c
c do max boundary points
c
        q1d(mjtot) = qnorth(i-2)
        u1d(mjtot) = vnorth(i-2)
        v1d(mjtot) = unorth(i-2)
#else
        q1d(1) = q(i-2,jdim-1)
        u1d(1) = vv(i-2,jdim-1)
        v1d(1) = uu(i-2,jdim-1)
c
c do max boundary points
c
        q1d(mjtot) = q(i-2,2)
        u1d(mjtot) = vv(i-2,2)
        v1d(mjtot) = uu(i-2,2)
#endif
c
c        # note: u, v, and w are rotated for solving in the y-direction
c        # The roles of fadd, gadd, and hadd are also rotated.
c
         call flux2d(q1d,u1d,v1d,deladd,fadd,gadd,1,mjtot)
c
         do 80 j=3,mjtot-1
           delta(i,j) = delta(i,j) + deladd(j)
           g(i,j) = g(i,j) + fadd(j)
           f(i,j) = f(i,j) + gadd(j,1)
           f(i+1,j) = f(i+1,j) + gadd(j,2)
   80    continue
  100 continue
c
c
c     # update q values by differencing fluxes
c     ========================================
c

      do  j = 3,mjtot-2
        do  i = 3,mitot-2
          qout(i-2,j-2) = q(i-2,j-2) 
     &                    - dtdx * (f(i+1,j) - f(i,j)
     &                    + g(i,j+1) - g(i,j) + delta(i,j))
        enddo
      enddo
c
c update boundary points
c
      call imag2d(qout)
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
c
c
c     =================================================
      subroutine flux2d(q,u,v,deladd,fadd,gadd,i1,i2)
c     =================================================
#include "param.inc"
#include "olemp.inc"

      integer msize

      parameter ( msize = idim+4 )
      real q(msize),fadd(msize),gadd(msize,2),u(msize),v(msize)
      real deladd(msize),dq(msize),gincr,fincr,fract,fincr3
      real dq2
      integer i1,i2,i,j,jside,method,mthlim
      common /mth/ method,mthlim
c
c     # compute fluxes along a single grid line in the x or y direction
c     # In this routine q is a one-dimensional array that has been set
c     # to the values along one line:   q(i), i running from i1 to i2
c     # u and v have been set to the velocities at the "left" edge of each cell. 
c     
c     # On return, 
c     #    deladd(i) is the increment to add to cell i
c     #    fadd(i) is the new contribution to the flux at the "left" edge,
c     #    gadd(i,1) is the contribution to the flux in the other direction
c     #              in the cell "below" this cell
c     #    gadd(i,2) is the contribution to the flux in the other direction
c     #              in the cell "above" this cell
c
c     # Methods:
c     # --------
c     #   method = 1:  donor-cell upwind
c     #   method = 2:  corner transport upwind
c     #   method = 3:  second order correction terms added 
c     #   method = 4:  with transverse propagation of 2nd order corrections
c     #   method = 5:  third order correction terms added
c     #   method = 6:  with transverse propagation of 3rd order corrections
c     #   method = 7:  Lax-Wendroff
c
c     # initialize fluxes:
c     --------------------
c
      do 20 i=i1,i2
	 deladd(i) = 0.d0
	 fadd(i) = 0.d0
         do 10 jside=1,2
	    gadd(i,jside) = 0.d0
   10       continue
   20    continue
c
c
      do 250 i=i1+1,i2
	 dq(i) = q(i) - q(i-1)
c
 	 if (u(i).gt.0.d0) then
c            # rightward going
             j = i
	    else
c            # leftward going
	     j = i-1
	    endif

c        # Basic first order increment for this cell based on piecewise 
c        # constant dq propagating normal to interface:
c        -------------------------------------------------------------
c
	 deladd(j) = deladd(j) + u(i)*dq(i)
	 if (method.eq.1) go to 250
c
	 if (method.eq.7) then
c           # Lax-Wendroff:
	    gincr = 0.125d0 * dtdx * u(i)*v(i)*dq(i)
	    gadd(i,1) = gadd(i,1) - gincr
	    gadd(i,2) = gadd(i,2) - gincr
	    gadd(i-1,1) = gadd(i-1,1) - gincr
	    gadd(i-1,2) = gadd(i-1,2) - gincr
	    go to 250
	    endif
c
c        # corrections for transverse propagation of dq:
c        -------------------------------------------------
	 if (v(i).gt.0.d0) then
c             # upward
	      gadd(j,2) = gadd(j,2) - 0.5d0 * dtdx * 
     &				u(i)*v(i)*dq(i)
	    else
c             # downward
	      gadd(j,1) = gadd(j,1) - 0.5d0 * dtdx * 
     &				u(i)*v(i)*dq(i)
  	    endif
  250     continue
c
      if (method.le.2) return  
c
c     # modify F fluxes for second order u_{xx} terms:
c     ------------------------------------------------
c
      call limiter(dq,u,i1,i2)
      do 200 i=i1+2,i2-1
	 fincr = 0.5d0*abs(u(i))*(1.d0 - abs(u(i))*dtdx) * dq(i)
	 fadd(i) = fadd(i) + fincr 
	 if (method.eq.3) go to 200
	 if (method.eq.7) go to 200
c
	 if (method.ge.5) then
c           # third order correction:
	    if (u(i) .gt. 0.d0) then
	        dq2 = dq(i) - dq(i-1)
	      else
	        dq2 = dq(i+1) - dq(i)
	      endif
	    fincr3 = -u(i)/6.d0 * (1.d0 - (u(i)*dtdx)**2) * dq2
	    fadd(i) = fadd(i) +  fincr3
	    if (method.eq.6) fincr = fincr + fincr3
	    endif
c
c	 # compute fraction of dq that goes into adjacent row
c        ----------------------------------------------------
c        # (transverse propagation of second order and possibly third 
c        # order correction waves)
c
         fract = abs(v(i))*dtdx  * fincr
	 if (v(i).lt.0.d0) then
c            # downward:
	     gadd(i,1) = gadd(i,1) - fract
	     gadd(i-1,1) = gadd(i-1,1) + fract
	   else
c            # upward:
	     gadd(i,2) = gadd(i,2) + fract
	     gadd(i-1,2) = gadd(i-1,2) - fract
	   endif
  200    continue
      return
      end
c
c
#if defined(NEWLEV)
c
c     =================================================
      subroutine limiter(dq,uu,i1,i2)
c     =================================================
#include "param.inc"
      integer msize
      parameter ( msize = idim+4 )
      real dq(msize),uu(msize)
      integer i1,i2,method,mthlim
      integer i
      real phi,philim,dqm1
      common /mth/ method,mthlim
c
      if (mthlim.eq.0) return
c
      dqm1 = dq(i1+1)
      do 100 i=i1+2,i2-1
        if (dq(i).eq.0.d0) go to 100
        if (uu(i).gt.0.d0) then
          phi =  philim(dq(i), dqm1, mthlim)
        else
          phi =  philim(dq(i), dq(i+1), mthlim)
        endif
        dqm1 = dq(i)
        dq(i) = dq(i)*phi
  100   continue
c
      return
      end
c
c
c
c
c
c
c     =====================================================
      real function philim(a,b,meth)
c     =====================================================
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
      real r,a,b,c
      integer meth
c
      r = b/a

c
      if(meth.eq.1) then
c     --------
c     # minmod
c     --------
         philim = max(0.0, min(1.0, r))
c
      else if (meth.eq.2) then
c
c     ----------
c     # superbee
c     ----------
         philim = max(0.0, min(1.0, 2.0*r), min(2.0, r))
c
      else if (meth.eq.3) then
c
c     ----------
c     # van Leer
c     ----------
         philim = (r + abs(r)) / (1.0 + abs(r))
c
      else
c     ----------------------
c     # monotinized centered 
c     ----------------------
         c = (1.0 + r)/2.0
         philim = max(0.0, min(c, 2.0, 2.0*r))
c
      endif
      return
c
      end
#endif
