#if !defined(SNGLSALT)
c
c     =================================================
      subroutine leveqdp(n,q,qout,ivar)
c     =================================================
c
c   $Id: leveq77dp.F,v 1.6 1998/10/16 22:12:35 dwd Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
c
      integer msize,mitot,mjtot,mktot
      parameter ( msize = idim+4 )
      parameter ( mitot = idim+4,mjtot = jdim+4, mktot = kdim+4)
c
      real*8 q(0:idim1,0:jdim1,0:kdim1)
      real*8 qout(0:idim1,0:jdim1,0:kdim1)
      real u1d(msize),v1d(msize),w1d(msize)
      real*8 f(mitot,mjtot,mktot),g(mitot,mjtot,mktot)
      real*8 h(mitot,mjtot,mktot)
      real*8 delta(mitot,mjtot,mktot),deladd(msize),q1d(msize)
      real*8 fadd(msize),gadd(msize,2),hadd(msize,3,2)
      real*8 kright,kleft,kback,ktop,kbot,sgs,b1,b2,b3,kfront
#if defined(ANALY)
      real*8 pel(kdim),zz
      real den(0:idim1,0:jdim1,0:kdim1)
#endif
      integer method,mthlim
      integer i,j,k,ivar,n
      common /mth/ method,mthlim
c
      if(ivar.ne.MDSAL) then
         write(*,*) "Calling leveq77dp with var other than salt"
         stop
      endif
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
        call fluxesdp(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
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
      do 100 k = 2,mktot-1
       do 100 i = 2,mitot-1
        do j=2,mjtot-1
            q1d(j) = q(i-2,j-2,k-2)
            u1d(j) = v(i-2,j-2,k-2)
            v1d(j) = w(i-2,j-2,k-2)
            w1d(j) = u(i-2,j-2,k-2)
        enddo
c
c do min boundary points
c
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
c
c        # note: u, v, and w are rotated for solving in the y-direction
c        # The roles of fadd, gadd, and hadd are also rotated.
c
         call fluxesdp(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
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
        q1d(1) = q(i-2,j-2,0) - dsdzlow
        u1d(1) = w(i-2,j-2,1)
        v1d(1) = u(i-2,j-2,1)
        w1d(1) = v(i-2,j-2,1)
c
c do max boundary points
c
        q1d(mktot) = q(i-2,j-2,kdim)
        u1d(mktot) = w(i-2,j-2,kdim)
        v1d(mktot) = u(i-2,j-2,kdim)
        w1d(mktot) = v(i-2,j-2,kdim)
c
c        # note: u, v, and w are rotated for solving in the z-direction
c        # The roles of fadd, gadd, and hadd are also rotated.
c
         call fluxesdp(q1d,u1d,v1d,w1d,deladd,fadd,gadd,hadd,1,
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

      do k = 3,mktot-2
         do j = 3,mjtot-2
            do i = 3,mitot-2
               qout(i-2,j-2,k-2) = q(i-2,j-2,k-2) 
     &              - dtdx * (f(i+1,j,k) - f(i,j,k)
     &              + g(i,j+1,k) - g(i,j,k) 
     &              + h(i,j,k+1) - h(i,j,k) + delta(i,j,k))
            enddo
         enddo
      enddo
c finish up the processing by mixing through turbulence parameterization.
c Don't do this from the TKE fields.  TKE is diffused in the sgs update
c routine.
c
c      call imagedp(qout)
      call topbotdp(qout,MDNULL,n)
c
c
c compute change in p.e. from advection
c
#if defined(ANALY)
c
c compute density perturbation, only works for constant salinity case
c
      call densitydp(th,qout,den)
c
         do k=1,kdim
            pel(k) = 0.0
            zz = (k-1)*dz + 0.5*dz
            do j=1,jdim
              do i=1,idim
                pel(k)=pel(k)+(den(i,j,k)+rhobm1(k))*zz*9.81/rhob(k)
              enddo  
            enddo 
         enddo
c 
c compute average rate of change due to advection term
c 
         do k=1,kdim 
            dpea(k) = dpea(k)+(pel(k)-pen(k))/delt 
         enddo 
#endif 
c     call imagdp(qout)
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
      do j=1,jdim
        do i=1,idim
           do k=1,kdim
              kright = kheat(i+1,j,k)+kheat(i,j,k)
              kleft = kheat(i,j,k)+kheat(i-1,j,k)
              b1=(kright*(q(i+1,j,k)-q(i,j,k))-
     *           kleft*(q(i,j,k)-q(i-1,j,k)))*0.5*odx*odx
c
              kback = kheat(i,j+1,k)+kheat(i,j,k)
              kfront = kheat(i,j,k)+kheat(i,j-1,k)
              b2=(kback*(q(i,j+1,k)-q(i,j,k))-
     *             kfront*(q(i,j,k)-q(i,j-1,k)))*0.5*ody*ody
c
              ktop = kheat(i,j,k+1)+kheat(i,j,k)
              kbot = kheat(i,j,k)+kheat(i,j,k-1)
#if defined(RADIATE)
              if(k.eq.1)then
                kbot = 0.0
              endif
#endif
              b3=(ktop*(q(i,j,k+1)-q(i,j,k))-
     *             kbot*(q(i,j,k)-q(i,j,k-1)))*0.5*odz*odz
              sgs=b1+b2+b3
              qout(i,j,k)=qout(i,j,k)+delt*sgs
           enddo
        enddo
c
c update boundary values
c
      enddo
c      call imagedp(qout)
      call topbotdp(qout,ivar,n)
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
      subroutine fluxesdp(q,u1,v1,w1,deladd,fadd,gadd,hadd,i1,i2)
c     =================================================
#include "param.inc"
#include "olemp.inc"
c
      integer msize
      parameter ( msize = idim+4 )
      real*8 q(msize),fadd(msize),gadd(msize,2)
      real u1(msize),v1(msize),w1(msize)
      real*8 deladd(msize),dq(msize),hadd(msize,3,2)
      real*8 corner,fincr,dq2,fincr3,fract,fract2,fractz

      integer method,mthlim
      integer i1,i2,i,j,k,kside,jside
c
      common /mth/ method,mthlim
c
c
c     # initialize fluxes:
c     --------------------
c
      do i=i1,i2
         deladd(i) = 0.d0
         fadd(i) = 0.d0
         do 10 jside=1,2
           gadd(i,jside) = 0.d0
           do 10 kside=1,3
              hadd(i,kside,jside) = 0.d0
   10    continue
      enddo

c
      do 250 i=i1+1,i2
        dq(i) = q(i) - q(i-1)
c
        if (u1(i).gt.0.d0) then
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
        deladd(j) = deladd(j) + u1(i)*dq(i)
        if (method.eq.1) go to 250
c
c
c        # corrections for transverse propagation of dq:
c        -------------------------------------------------
        if (w1(i).gt.0.d0) then
          k = 2
        else
          k = 1
        endif
        corner = 1.d0/3.d0 * u1(i)*abs(v1(i))*w1(i)*(dtdx)**2 * dq(i)
        hadd(j,2,k) = hadd(j,2,k) - 0.5d0*u1(i)*w1(i)*dtdx * dq(i) 
     &       + corner
c
        if (v1(i).gt.0.d0) then
c             # upward
          gadd(j,2) = gadd(j,2) - 0.5d0 * dtdx * 
     &                u1(i)*v1(i)*dq(i)
          hadd(j,3,k) = hadd(j,3,k) - corner 
        else
c             # downward
          gadd(j,1) = gadd(j,1) - 0.5d0 * dtdx * 
     &                u1(i)*v1(i)*dq(i)
          hadd(j,1,k) = hadd(j,1,k) - corner 
        endif
  250 continue
c
      if (method.le.2) return  
c
c     # modify F fluxes for second order u_{xx} terms:
c     ------------------------------------------------
c
      if (mthlim.gt.0) call limiterdp(dq,u1,i1,i2)
      do 200 i=i1+2,i2-1
        if (w1(i).gt.0.d0) then
          k = 2
        else
          k = 1
        endif
        fincr = 0.5d0*abs(u1(i))*(1.d0 - abs(u1(i))*dtdx) * dq(i)
        fadd(i) = fadd(i) + fincr 
        if (method.eq.3) go to 200
c
        if (method.ge.5) then
c           # third order correction:
          if (u1(i) .gt. 0.d0) then
            dq2 = dq(i) - dq(i-1)
          else
            dq2 = dq(i+1) - dq(i)
          endif
          fincr3 = -u1(i)/6.d0 * (1.d0 - (u1(i)*dtdx)**2) * dq2
          fadd(i) = fadd(i) +  fincr3
          if (method.eq.6) fincr = fincr + fincr3
        endif
c
c	 # compute fraction of dq that goes into adjacent row
c        ----------------------------------------------------
c        # (transverse propagation of second order and possibly third 
c        # order correction waves)
c
        fract = abs(v1(i))*dtdx  * fincr
        fract2 = fincr - fract
        fractz = w1(i)*dtdx
        if (v1(i).lt.0.d0) then
c            # downward:
          gadd(i,1) = gadd(i,1) - fract
          gadd(i-1,1) = gadd(i-1,1) + fract
          hadd(i,1,k) = hadd(i,1,k) + fract*fractz
          hadd(i-1,1,k) = hadd(i-1,1,k) - fract*fractz
          hadd(i,2,k) = hadd(i,2,k) + fract2*fractz
          hadd(i-1,2,k) = hadd(i-1,2,k) - fract2*fractz
        else
c            # upward:
          gadd(i,2) = gadd(i,2) + fract
          gadd(i-1,2) = gadd(i-1,2) - fract
          hadd(i,3,k) = hadd(i,3,k) + fract*fractz
          hadd(i-1,3,k) = hadd(i-1,3,k) - fract*fractz
          hadd(i,2,k) = hadd(i,2,k) + fract2*fractz
          hadd(i-1,2,k) = hadd(i-1,2,k) - fract2*fractz
        endif
c
  200 continue
      return
      end
c
c
c
c
c     =================================================
      subroutine limiterdp(dq,uu,i1,i2)
c     =================================================
#include "param.inc"
      integer msize
      parameter ( msize = idim+4 )
      real*8 dq(msize)
      real uu(msize)
      integer i1,i2,method,mthlim
      integer i
      real*8 phi,philimdp,dqm1
      common /mth/ method,mthlim
c
      if (mthlim.eq.0) return
c
      dqm1 = dq(i1+1)
      do 100 i=i1+2,i2-1
        if (dq(i).eq.0.d0) go to 100
        if (uu(i).gt.0.d0) then
          phi =  philimdp(dq(i), dqm1, mthlim)
        else
          phi =  philimdp(dq(i), dq(i+1), mthlim)
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
      real*8 function philimdp(a,b,meth)
c     =====================================================
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
      real*8 r,a,b,c
      integer meth
c
      r = b/a

c
      if(meth.eq.1) then
c     --------
c     # minmod
c     --------
         philimdp = dmax1(0.d0, dmin1(1.d0, r))
c
      else if (meth.eq.2) then
c
c     ----------
c     # superbee
c     ----------
         philimdp = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
c
      else if (meth.eq.3) then
c
c     ----------
c     # van Leer
c     ----------
         philimdp = (r + dabs(r)) / (1.d0 + dabs(r))
c
      else
c     ----------------------
c     # monotinized centered 
c     ----------------------
         c = (1.0 + r)/2.0
         philimdp = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
c
      endif
      return
c
      end
#else
      subroutine leveqdp(n,q,qout,ivar)
      integer n,ivar
      real*8 q,qout
      end
#endif
