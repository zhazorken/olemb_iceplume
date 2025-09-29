      Subroutine conres(n,forc,nerr)
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
c
      integer n
      real forc(0:idim1,0:jdim1,0:kdim1)
#include "olemf.inc"
#include "olemp.inc"
#include "communicate.inc"
c
c This routine is a conjugate-residual solver as outlined in
c
c
      real r(0:idim1,0:jdim1,0:kdim1),ar(0:idim1,0:jdim1,0:kdim1)
      real qr(0:idim1,0:jdim1,0:kdim1)
      real x(0:idim1,0:jdim1,0:kdim1,8)
      real ax(0:idim1,0:jdim1,0:kdim1,8)
      real del(8),axar(8),axarl(8)
      real beta,alpha,pmax,rr,rms_residual,rrold,pbarr
      real rrl,pbarl,eps,dvmx
      real dvmxl,pmaxl
      real*8 top, topl, bot(8), botl(8),epa
      integer mscan,mxscan,ncheck,i,j,k,m,nerr
      integer lord,iter,ll
#include "moddef.inc"
c
c initialize fields and varibles
      mxscan = 200
      lord = 4 
      ncheck = 30
      eps = 1.e-5
      epa = 1.e-30
      pmax = 0.0
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            pmaxl = max(pmax,abs(p(i,j,k)))
            do m=1,lord
              x(i,j,k,m) = 0.0
              ax(i,j,k,m) = 0.0
            enddo
            r(i,j,k) = 0.0
            ar(i,j,k) = 0.0
          enddo
        enddo
      enddo
#if defined(USEMPI)
         call MPI_ALLREDUCE(pmaxl,pmax,1,
     &        MPI_REAL,MPI_MAX,MPI_COMM_WORLD, nerr)
#else
        pmax = pmaxl
#endif
c
      if(iope) then
        write(*,*)"Pmax is ",pmax
      endif
      call operator(r,p)
c
c compute residual, assumes the forc is opposite sign
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            r(i,j,k) = forc(i,j,k) + r(i,j,k)
          enddo
        enddo
      enddo
c
c
      call imag(r)
      call rboun(r)
      call precon(r,qr,1)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            x(i,j,k,1) =  qr(i,j,k)
          enddo
        enddo
      enddo
      call operator(ar,qr)
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            ax(i,j,k,1) = ar(i,j,k)
          enddo
        enddo
      enddo

c
c start iteration
c
      do iter=1,mxscan
      do m=1,lord
c
        topl = 0.0
        botl(m) = 0.0
c
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              topl = topl + r(i,j,k)*ax(i,j,k,m)
              botl(m) = botl(m) + ax(i,j,k,m)*ax(i,j,k,m)
            enddo
          enddo
        enddo
#if defined(USEMPI)
         call MPI_ALLREDUCE(topl,top,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, nerr)
         call MPI_ALLREDUCE(botl(m),bot(m),1, 
     &        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, nerr)
#else
        top = topl
        bot(m)= botl(m)
#endif
c
        bot(m) = dmax1(epa,bot(m))
        beta = -top/bot(m)
        dvmxl = -1.0e15
        do k=1,kdim 
          do j=1,jdim 
            do i=1,idim 
              p(i,j,k) = p(i,j,k) + beta*x(i,j,k,m)
              r(i,j,k) = r(i,j,k) + beta*ax(i,j,k,m)
              dvmxl = amax1(dvmxl,abs(r(i,j,k)))
            enddo 
          enddo 
        enddo
c
#if defined(USEMPI)
         call MPI_ALLREDUCE(dvmxl,dvmx,1,
     &        MPI_REAL,MPI_MAX,MPI_COMM_WORLD, nerr)
#else
        dvmx = dvmxl
#endif

        if(dvmx.le.eps) then
          goto 100
        else
c          write(*,*)"dvmx ",iter,dvmx
        endif
        call imag(r)
        call rboun(r)
c
        call precon(r,qr,1)
        call operator(ar,qr)
c
        do ll=1,m
          axarl(ll)= 0.0
        
          do k=1,kdim 
            do j=1,jdim 
              do i=1,idim 
                axarl(ll) = axarl(ll) + ar(i,j,k)*ax(i,j,k,ll)
              enddo 
            enddo 
          enddo
#if defined(USEMPI)
           call MPI_ALLREDUCE(axarl(ll),axar(ll),1,
     &          MPI_REAL,MPI_SUM,MPI_COMM_WORLD, nerr)
#else
          axar(ll) = axarl(ll)
#endif
          del(ll) = -axar(ll)/bot(ll)
        enddo
        if(m.lt.lord) then
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                x(i,j,k,m+1) = qr(i,j,k)
                ax(i,j,k,m+1) = ar(i,j,k)
              enddo
            enddo
          enddo

          do ll=1,m
          do k=1,kdim 
            do j=1,jdim 
              do i=1,idim 
                x(i,j,k,m+1) = x(i,j,k,m+1) + del(ll)*x(i,j,k,ll)
                ax(i,j,k,m+1) = ax(i,j,k,m+1)+del(ll)*ax(i,j,k,ll)
              enddo 
            enddo 
          enddo
          enddo
        else
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                x(i,j,k,1) = qr(i,j,k) + del(1)*x(i,j,k,1)
                ax(i,j,k,1) = ar(i,j,k)+del(1)*ax(i,j,k,1)
              enddo
            enddo
          enddo
          do ll=2,m
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                x(i,j,k,1) = x(i,j,k,1) + del(ll)*x(i,j,k,ll)
                x(i,j,k,ll) = 0.0
                ax(i,j,k,1) = ax(i,j,k,1)+del(ll)*ax(i,j,k,ll)
                ax(i,j,k,ll) = 0.0
              enddo
            enddo
          enddo
          enddo
        endif
      enddo

c
      enddo
c
100   continue
      if(iope) then
        write(*,*)"Converged at ",iter, " iterations ",dvmx
      endif
      pbarl = 0.0
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            pbarl = pbarl + p(i,j,k)
          enddo
        enddo
      enddo
#if defined(USEMPI)
      call MPI_ALLREDUCE(pbarl,pbarr,1, 
     &        MPI_REAL,MPI_SUM,MPI_COMM_WORLD, nerr)
#else
      pbarr = pbarl
#endif
      pbarr = pbarr/float(kdim*idim*jtotal)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            p(i,j,k) = p(i,j,k) - pbarr
          enddo
        enddo
      enddo
      if(iter.ge.mxscan) then
          write(*,*)' No Convergence ',dvmx
          nerr = 1
#if defined(USEMPI)
          call MPI_ABORT(MPI_COMM_WORLD,1,mpierr)
#endif
c      else
c          if(iope)then
c           write(*,*)"m and dvmx ",m,dvmx
c          endif
      endif
      return
      end
c
c
c Define operator
c
      subroutine operator(r,p)
#include "param.inc"
c
c Define laplacian operator for the incompressible system
c
      real p(0:idim1,0:jdim1,0:kdim1),r(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k
#include "olemp.inc"
#include "finitv.inc"
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
c            r(i,j,k) = (p(i+1,j,k)+p(i-1,j,k)+p(i,j+1,k)+p(i,j-1,k)+
c     *                  p(i,j,k+1)+p(i,j,k-1)-6.*p(i,j,k))*tdxs
            r(i,j,k)=(ax(i+1,j,k)**2*(p(i+1,j,k)-p(i,j,k))*
     *            volu(i+1,j,k)
     *           -ax(i,j,k)**2*(p(i,j,k)-p(i-1,j,k))*volu(i,j,k))*tdxs
     *           +(ay(i,j+1,k)**2*(p(i,j+1,k)-p(i,j,k))*volv(i,j+1,k)
     *           -ay(i,j,k)**2*(p(i,j,k)-p(i,j-1,k))*volv(i,j,k))*tdxs
     *           + ((p(i,j,k+1)-p(i,j,k))*volw(i,j,k+1) - (p(i,j,k)
     *           - p(i,j,k-1))*volw(i,j,k))*tdxs

          enddo
        enddo
      enddo
      return
      end
c
c
      subroutine rboun(r)
#include "param.inc"
#include "communicate.inc"
      real r(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k
c
#if defined(OPENBC)
      if(blockstart(my_pe+1).eq.1) then
        do k=0,kdim1
          do i=0,idim1
            r(i,0,k) = r(i,1,k)
c            r(i,0,k) = 0.0
          enddo
        enddo
      endif
      if((blockstart(my_pe+1)+jdim).ge.jtotal) then
        do k=0,kdim1
          do i=0,idim1
            r(i,jdim1,k) = r(i,jdim,k)
c            r(i,jdim1,k) = 0.0
          enddo
        enddo
      endif
#endif
#if defined(OPENBCX)
      do k=0,kdim1
        do j=0,jdim1
c          r(0,j,k) = r(1,j,k)
          r(0,j,k) = 0.0
c          r(idim1,j,k) = r(idim,j,k)
          r(idim1,j,k) = 0.0
        enddo
      enddo
#endif
      do j=0,jdim1
        do i=0,idim1
           r(i,j,0) = 0.0
           r(i,j,kdim1) = 0.0
c           r(i,j,kdim1) = r(i,j,kdim)
c           r(i,j,0) = r(i,j,1)
        enddo
      enddo
      return
      end
c
c Richardson Preconditioner
c
      subroutine precon(r,e,iflg)
#include "param.inc"
#include "olemp.inc"
#include "finitv.inc"
      real r(0:idim1,0:jdim1,0:kdim1)
      real e(0:idim1,0:jdim1,0:kdim1)
      real ae(0:idim1,0:jdim1,0:kdim1)
      real beta,junk
      real ee(kdim),lfh(kdim)
      integer ntr,i,j,k,it,iflg
c
      if(iflg.eq.0) then
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              e(i,j,k)=r(i,j,k)
            enddo
          enddo
        enddo
      else

c
        beta=0.055*(dx**2+dy**2+dz**2)
c
c number of passes
c
        ntr= 1
  
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              ae(i,j,k)=0.
              e(i,j,k) =  beta*r(i,j,k)
            enddo
          enddo
        enddo
c
        do it=1,ntr
c          call operator(ae,e)
              do k=1,kdim
          do j=1,jdim
            do i=1,idim
c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)
c     *                  +e(i,j+1,k)+e(i,j-1,k)+
c     *                  e(i,j,k+1)+e(i,j,k-1)-6.*e(i,j,k))*tdxs
            ae(i,j,k)=(ax(i+1,j,k)**2*(e(i+1,j,k)-e(i,j,k))*
     *            volu(i+1,j,k)
     *           -ax(i,j,k)**2*(e(i,j,k)-e(i-1,j,k))*volu(i,j,k))*tdxs
     *           +(ay(i,j+1,k)**2*(e(i,j+1,k)-e(i,j,k))*volv(i,j+1,k)
     *           -ay(i,j,k)**2*(e(i,j,k)-e(i,j-1,k))*volv(i,j,k))*tdxs
     *           + ((e(i,j,k+1)-e(i,j,k))*volw(i,j,k+1) - (e(i,j,k)
     *           - e(i,j,k-1))*volw(i,j,k))*tdxs
c
                e(i,j,k)=e(i,j,k)+beta*(ae(i,j,k)+r(i,j,k))

c                lfh(k) = beta*((e(i+1,j,k)+e(i-1,j,k)
c     *                  +e(i,j+1,k)+e(i,j-1,k) +
c     *                  2.*e(i,j,k))*tdxs+r(i,j,k))
c
c              enddo
c              call tridag(lfh,ee,kdim)
c              do k=1,kdim
c                e(i,j,k) = ee(k)
              enddo
            enddo
          enddo
          call imag(e)
          call rboun(e)
        enddo
      endif
  
       
      return
      end
c
c tridiagonal solver
c
      subroutine tridag(r,u,n)
      integer n
#include "param.inc"
      real r(kdim),u(kdim)
      integer k
      real a, bet, gam(kdim)
c
      a = 1./6.
c
      bet = 4./3.
      u(1) = r(1)/bet
      do k=2,kdim
        gam(k) = -a/bet
        bet = 4./3.+a*gam(k)
        if(bet.eq.0.0) write(*,*)"Trouble in Tridag"
        u(k) = (r(k)+a*u(k-1))/bet
      enddo
      do k=n-1,1,-1
        u(k) = u(k) - gam(k+1)*u(k+1)
      enddo
      return
      end
      
