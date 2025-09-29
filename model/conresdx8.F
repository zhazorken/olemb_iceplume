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
c This routine is a conjugate-residual solver as outlined in Smolarkiewicz
c and Margolin, Appl. Math. & Comp. Sci., 4, 527-551.
c
c n is the time step
c forc is the right hand side of the equation
c operator is the function to be solved 
c precon is a user set preconditioner (currently set to a 1 pass Richardson
c solver).
c
c
c
      real*8 r(0:idim1,0:jdim1,0:kdim1),ar(0:idim1,0:jdim1,0:kdim1)
      real*8 qr(0:idim1,0:jdim1,0:kdim1)
      real*8 pp(0:idim1,0:jdim1,0:kdim1)
      real*8 x(0:idim1,0:jdim1,0:kdim1,8)
      real*8 ax(0:idim1,0:jdim1,0:kdim1,8)
      real*8 del(8)
      real*8 axar(8),axarl(8)
      real*8 alpha,pmax,rr,rms_residual,rrold,pbarr
      real*8 rrl,pbarl,eps,dvmx
      real*8 dvmxl,pmaxl
      real*8 top, topl, bot(8), botl(8),epa,beta
      integer mscan,mxscan,ncheck,i,j,k,m,nerr,mpierr
      integer lord,iter,ll
#include "moddef.inc"
c
c initialize fields and varibles
      mxscan = 200
      lord = 4 
      ncheck = 30
c
c eps sets the iteration limit
c
      eps = 0.01
      epa = 1.e-30
      pmax = 0.0
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            pmaxl = max(pmax,abs(p(i,j,k)))
            pp(i,j,k) = p(i,j,k)
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
     &        MPI_REAL8,MPI_MAX,comm, nerr)
#else
        pmax = pmaxl
#endif
c
      if(iope) then
        write(*,*)"Pmax is ",pmax
      endif
      call operator(r,pp)
c
c compute residual, assumes the forc is opposite sign
c
#if defined(RECIRC)
      call imag8(r)
#endif
      dvmxl = -1.e-15
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            r(i,j,k) = forc(i,j,k) + r(i,j,k)
c            tke(i,j,k) = r(i,j,k)
            dvmxl = dmax1(dvmxl,abs(r(i,j,k)))
          enddo
        enddo
      enddo
#if defined(USEMPI)
         call MPI_ALLREDUCE(dvmxl,dvmx,1,
     &        MPI_REAL8,MPI_MAX,comm, nerr)
#else
        dvmx = dvmxl
#endif
        if(iope) write(*,*)"First dvmx is ",dvmx

c
c
      call imag8(r)
#if defined(FRONTO)
      do k=1,kdim
        do j=0,jdim1
          r(0,j,k) = r(1,j,k)
          r(idim1,j,k) = r(idim,j,k)
        enddo
      enddo
#endif
      call rboun(r)
      call precon(r,qr,1)
      call rboun(qr)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            x(i,j,k,1) =  qr(i,j,k)
          enddo
        enddo
      enddo
      call imag8(qr)
      call operator(ar,qr)
#if defined(RECIRC)
      call imag8(ar)
#endif
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
     &        MPI_REAL8,MPI_SUM,comm, nerr)
         call MPI_ALLREDUCE(botl(m),bot(m),1, 
     &        MPI_REAL8,MPI_SUM,comm, nerr)
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
              pp(i,j,k) = pp(i,j,k) + beta*x(i,j,k,m)
              r(i,j,k) = r(i,j,k) + beta*ax(i,j,k,m)
              dvmxl = dmax1(dvmxl,abs(r(i,j,k)))
            enddo 
          enddo 
        enddo
c
#if defined(USEMPI)
         call MPI_ALLREDUCE(dvmxl,dvmx,1,
     &        MPI_REAL8,MPI_MAX,comm, nerr)
#else
        dvmx = dvmxl
#endif
c        if(iope) write(*,*)"dvmx is ",dvmx
        if(dvmx.le.eps) then
          goto 100
c        else
          write(*,*)"dvmx ",iter,dvmx
        endif
        call imag8(r)
#if defined(FRONTO)
      do k=1,kdim
        do j=0,jdim1
          r(0,j,k) = r(1,j,k)
          r(idim1,j,k) = r(idim,j,k)
        enddo
      enddo
#endif
        call rboun(r)
c
        call precon(r,qr,1)
        call imag8(qr)
        call rboun(qr)
        call operator(ar,qr)
c
#if defined(RECIRC)
        call imag8(ar)
#endif
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
     &          MPI_REAL8,MPI_SUM,comm, nerr)
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
        write(*,*)"Converged at ",iter, " iterations ",dvmx, " order ",m
      endif
      pbarl = 0.0
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            p(i,j,k) = pp(i,j,k)
            pbarl = pbarl + pp(i,j,k)
          enddo
        enddo
      enddo
#if defined(USEMPI)
      call MPI_ALLREDUCE(pbarl,pbarr,1, 
     &        MPI_REAL8,MPI_SUM,comm, nerr)
#else
      pbarr = pbarl
#endif
      pbarr = pbarr/float(kdim*itotal*jtotal)
c      do k=1,kdim
c        do j=1,jdim
c          do i=1,idim
c            p(i,j,k) = p(i,j,k) - pbarr
c          enddo
c        enddo
c      enddo
      call imag(p)
c      call rboun(p)
#if defined(FRONTO)
      do k=1,kdim
        do j=0,jdim1
          p(0,j,k) = p(1,j,k)
          p(idim1,j,k) = p(idim,j,k)
        enddo
      enddo
#endif
      do j=0,jdim1
        do i=0,idim1
          p(i,j,0) = p(i,j,1)
          p(i,j,kdim1) = p(i,j,kdim)
        enddo
      enddo
      if(iter.ge.mxscan) then
          write(*,*)' No Convergence ',dvmx
          nerr = 1
#if defined(USEMPI)
          call MPI_ABORT(comm,1,mpierr)
#endif
c      else
c          write(*,*)"m and dvmx ",m,dvmx
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
      real*8 p(0:idim1,0:jdim1,0:kdim1),r(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k
#include "olemp.inc"
#include "finitv.inc"
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim

#if defined(HTOP) 
c             r(i,j,k) = (p(i+1,j,k)+p(i-1,j,k)+p(i,j,k+1)+p(i,j,k-1)-
c     *                  4.*p(i,j,k))*tdxs + (vay(i,j+1,k)*(p(i,j+1,k)-
c     *                  p(i,j,k))-vay(i,j,k)*(p(i,j,k)-p(i,j-1,k)))*tdxs

             r(i,j,k) = (bax(i+1,j,k)*(p(i+1,j,k)-p(i,j,k))-
     *                  bax(i,j,k)*(p(i,j,k)-p(i-1,j,k))
     *                  +way(i,j,k+1)*(p(i,j,k+1)-p(i,j,k))-
     *                   way(i,j,k)*(p(i,j,k)-p(i,j,k-1))
     *                   +ay(i,j+1,k)*(p(i,j+1,k)-
     *                   p(i,j,k))-ay(i,j,k)*(p(i,j,k)-p(i,j-1,k)))*tdxs
#else
            r(i,j,k) = (p(i+1,j,k)+p(i-1,j,k)-2.*p(i,j,k))/
     *                 (dx(i)*dx(i))+(p(i,j+1,k)
     *                  +p(i,j-1,k)+
     *                  p(i,j,k+1)+p(i,j,k-1)-4.*p(i,j,k))*tdxs
#endif

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
#include "finitv.inc"
      real*8 r(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k
c
#if defined(OPENBC)
#if !defined(RECIRC)
      if(iblkstart(my_pe+1).eq.1) then
        do k=0,kdim1
          do j=0,jdim1
            r(0,j,k) = r(1,j,k)
c            r(i,0,k) = 0.0
          enddo
        enddo
      endif
#endif
      if((iblkstart(my_pe+1)+idim).ge.itotal) then
        do k=0,kdim1
          do j=0,jdim1
            r(idim1,j,k) = r(idim,j,k)
c            r(i,jdim1,k) = 0.0
          enddo
        enddo
      endif
#endif
#if defined(CHANNEL)
      if(jblkstart(my_pe+1).eq.1) then
        do k=0,kdim1
          do i=0,idim1
            r(i,0,k) = r(i,1,k)
          enddo
        enddo
      endif
      if((jblkstart(my_pe+1)+jdim).ge.jtotal) then
        do k=0,kdim1
          do i=0,idim1
            r(i,jdim1,k) = r(i,jdim,k)
c            r(i,jdim1,k) = 0.0
          enddo
        enddo
      endif
#endif

      do j=0,jdim1
        do i=0,idim1
#if defined(ICE)
c           do k=kdist(i,j)+1,kdim1
c             r(i,j,k) = 0.0
c           enddo
#endif
#if defined(HTOP)
            do k=1,kdim1
              if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
                  r(i-1,j,k) = r(i,j,k)
              endif
            enddo
#endif
c#if defined(HBOT)
c           do k=0, kdist(i,j)-1
c             r(i,j,k) = 0.0
c           enddo
c#endif
#if defined(RADTOP)
             r(i,j,kdim1) = 0.0
#else
             r(i,j,kdim1) = r(i,j,kdim)
#endif
c            r(i,j,kdim1) = 0.0
c            r(i,j,0) = r(i,j,1)
#if defined(RADIATE)
            r(i,j,0) = 0.0
#else
            r(i,j,0) = r(i,j,1)
#endif
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
      real*8 r(0:idim1,0:jdim1,0:kdim1)
      real*8 e(0:idim1,0:jdim1,0:kdim1)
      real*8 ae(0:idim1,0:jdim1,0:kdim1)
      real*8 beta,junk
      real*8 ee(kdim),lfh(kdim)
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
        beta=0.05*(dxx**2+dy**2+dz**2)
c
c number of passes
c
        ntr= 1
  
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              ae(i,j,k)=0.
              e(i,j,k) =  beta*r(i,j,k)
c              e(i,j,k) =  r(i,j,k)
c              e(i,j,k) =  0.5*r(i,j,k)/(2./dx(i)**2+4./(dy*dy))
            enddo
          enddo
        enddo
c
        do it=1,ntr
c          call operator(ae,e)
c
c do 1 ADI, horizontal first
c
         do k=1,kdim,2
          do j=1,jdim,2
            do i=1,idim,2

c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
c     *                  (dx(i)*dx(i))
c     *                  +(e(i,j+1,k)+e(i,j-1,k)+
c     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
c
c             ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)+e(i,j,k+1)+e(i,j,k-1)-
c     *                  4.*e(i,j,k))*tdxs + (vay(i,j+1,k)*(e(i,j+1,k)-
c     *                  e(i,j,k))-vay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs

#if defined(HTOP)
             ae(i,j,k) = (bax(i+1,j,k)*(e(i+1,j,k)-e(i,j,k))-
     *                  bax(i,j,k)*(e(i,j,k)-e(i-1,j,k))
     *                  +way(i,j,k+1)*(e(i,j,k+1)-e(i,j,k))-
     *                   way(i,j,k)*(e(i,j,k)-e(i,j,k-1))
     *                   +ay(i,j+1,k)*(e(i,j+1,k)-
     *                   e(i,j,k))-ay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs
#else
             ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
     *                  (dx(i)*dx(i))
     *                  +(e(i,j+1,k)+e(i,j-1,k)+
     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
#endif

                beta = 0.5/(2./(dx(i)*dx(i))+4./(dy*dy))

                e(i,j,k)=e(i,j,k)+beta
     *                 *(ae(i,j,k)+r(i,j,k))
         enddo
         enddo
         enddo
         do k=2,kdim,2
          do j=2,jdim,2
            do i=2,idim,2

c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
c     *                  (dx(i)*dx(i))
c     *                  +(e(i,j+1,k)+e(i,j-1,k)+
c     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
c
c             ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)+e(i,j,k+1)+e(i,j,k-1)-
c     *                  4.*e(i,j,k))*tdxs + (vay(i,j+1,k)*(e(i,j+1,k)-
c     *                  e(i,j,k))-vay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs

#if defined(HTOP)
             ae(i,j,k) = (bax(i+1,j,k)*(e(i+1,j,k)-e(i,j,k))-
     *                  bax(i,j,k)*(e(i,j,k)-e(i-1,j,k))
     *                  +way(i,j,k+1)*(e(i,j,k+1)-e(i,j,k))-
     *                   way(i,j,k)*(e(i,j,k)-e(i,j,k-1))
     *                   +ay(i,j+1,k)*(e(i,j+1,k)-
     *                   e(i,j,k))-ay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs
#else
            ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
     *                  (dx(i)*dx(i))
     *                  +(e(i,j+1,k)+e(i,j-1,k)+
     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
#endif


                beta = 0.5/(2./(dx(i)*dx(i))+4./(dy*dy))

                e(i,j,k)=e(i,j,k)+beta
     *                 *(ae(i,j,k)+r(i,j,k))


c
c code segment to do vertical direction implicitly
c
c might help in thin shell cases where dz is much smaller than dx,dy
c
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
          call imag8(e)
#if defined(FRONTO)
      do k=1,kdim
        do j=0,jdim1
          e(0,j,k) = e(1,j,k)
          e(idim1,j,k) = e(idim,j,k)
        enddo
      enddo
#endif

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
      
