      Subroutine sorvdx(n,forc,nerr)
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
      real*8 ae(0:idim1,0:jdim1,0:kdim1)
      real*8 e(0:idim1,0:jdim1,0:kdim1)
      real*8 beta(0:idim1,0:jdim1,0:kdim1)
      real*8 pp(0:idim1,0:jdim1,0:kdim1)
      real*8 del(8)
      real*8 axar(8),axarl(8)
      real*8 alpha,pmax,rr,rms_residual,rrold,pbarr
      real*8 rrl,pbarl,eps,dvmx
      real*8 dvmxl,pmaxl
      real*8 top, topl, bot(8), botl(8),epa, womega
      integer mscan,mxscan,ncheck,i,j,k,m,nerr,mpierr
      integer lord,iter,ll
#include "moddef.inc"
c
c initialize fields and varibles
      mxscan = 300
      lord = 4 
      ncheck = 30
c
c eps sets the iteration limit
c
c      womega = 2./(1.+sin(3.14159/480.))
      womega = 1.4

      eps = 3.0e-4
      eps = 3.0e-3
      epa = 1.e-30
      pmaxl = 0.0
      do k=0,kdim1
        do j=0,jdim1
          do i=1,idim1
c            pmaxl = max(pmaxl,abs(p(i,j,k)))
            pmaxl = pmaxl+abs(p(i,j,k))
            e(i,j,k) = p(i,j,k)
c            beta(i,j,k) = (2./(dx(i)*dx(i))+2./(dy*dy)+2./(dz*dz)
            beta(i,j,k) = 2./(dx(i)*dx(i-1))+2./(dy*dy)+2./(dz*dz)
c            beta(i,j,k) = 6./(dxx*dxx)
          enddo
        enddo
      enddo
#if defined(USEMPI)
         call MPI_ALLREDUCE(pmaxl,pmax,1,
     &        MPI_REAL8,MPI_SUM,comm, nerr)
#else
        pmax = pmaxl
#endif
      pmax = pmax/(itotal*jtotal*kdim)
c
      if(iope) then
        write(*,*)"Pave is ",pmax
      endif
c
      do it=1,mxscan
c
c
         
c         if(mod(it,2).eq.0) then

         dvmxl = 1.0e-15
         do k=1,kdim,1
          do j=1,jdim,1
            do i=1,idim,1

                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k))/(dx(i)*dx(i-1))
     *                  +(e(i,j+1,k)+e(i,j-1,k))*tdxs+
     *                  (e(i,j,k+1)+e(i,j,k-1))*tdxs

                corr = ((ae(i,j,k)+forc(i,j,k))/beta(i,j,k)-e(i,j,k))
                e(i,j,k) = e(i,j,k)+womega*corr
                dvmxl = dmax1(dvmxl,abs(corr))
            enddo
          enddo
         enddo
         call imag8(e)
         do i=0,idim1
          do j=0,jdim1
            e(i,j,0) = e(i,j,1)
            e(i,j,kdim1) = e(i,j,kdim)
          enddo
         enddo

c
c next do even cells with updated values
c
c         do k=2,kdim,2
c          do j=2,jdim,2
c            do i=2,idim,2
c
c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k))/(dx(i)*dx(i))
c     *                  +(e(i,j+1,k)+e(i,j-1,k))*tdxs+
c     *                  (e(i,j,k+1)+e(i,j,k-1))*tdxs
cc                 pp(i,j,k) = ae(i,j,k)-(4.*e(i,j,k)/tdxs+
cc     *                      2.*e(i,j,k)/(dx(i)*dx(i)))
cc
cc                r(i,j,k) = (ae(i,j,k)+dxx*dxx*forc(i,j,k))/6.-e(i,j,k)
cc                e(i,j,k) = (ae(i,j,k)+dxx*dxx*forc(i,j,k))/6.
c                corr = ((ae(i,j,k)+forc(i,j,k))/beta(i,j,k)-e(i,j,k))
c                e(i,j,k) = e(i,j,k)+womega*corr
c                dvmxl = dmax1(dvmxl,abs(corr))
cc
cc                e(i,j,k) = e(i,j,k)+womega*r(i,j,k)
c            enddo
c          enddo
c         enddo
c         call imag8(e)
c         do i=0,idim1
c          do j=0,jdim1
c            e(i,j,0) = e(i,j,1)
c            e(i,j,kdim1) = e(i,j,kdim)
c          enddo
c         enddo
cc
c         else
c
c flip red/black
c
c         do k=2,kdim,2
c          do j=2,jdim,2
c            do i=2,idim,2
c
c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k))/(dx(i)*dx(i))
c     *                  +(e(i,j+1,k)+e(i,j-1,k))*tdxs+
c     *                  (e(i,j,k+1)+e(i,j,k-1))*tdxs
c
c                corr = ((ae(i,j,k)+forc(i,j,k))/beta(i,j,k)-e(i,j,k))
c                e(i,j,k) = e(i,j,k)+womega*corr
c                dvmxl = dmax1(dvmxl,abs(corr))
c            enddo
c          enddo
c         enddo
c         call imag8(e)
c         do i=1,idim1
c          do j=1,jdim1
c            e(i,j,0) = e(i,j,1)
c            e(i,j,kdim1) = e(i,j,kdim)
c          enddo
c         enddo
c
cc
cc next do odd cells with updated values
cc
c         do k=1,kdim,2
c          do j=1,jdim,2
c            do i=1,idim,2
c
c                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k))/(dx(i)*dx(i))
c     *                  +(e(i,j+1,k)+e(i,j-1,k))*tdxs+
c     *                  (e(i,j,k+1)+e(i,j,k-1))*tdxs
cc                 pp(i,j,k) = ae(i,j,k)-(4.*e(i,j,k)/tdxs+
cc     *                      2.*e(i,j,k)/(dx(i)*dx(i)))
cc
cc                r(i,j,k) = (ae(i,j,k)+dxx*dxx*forc(i,j,k))/6.-e(i,j,k)
cc                e(i,j,k) = (ae(i,j,k)+dxx*dxx*forc(i,j,k))/6.
c                corr = ((ae(i,j,k)+forc(i,j,k))/beta(i,j,k)-e(i,j,k))
c                e(i,j,k) = e(i,j,k)+womega*corr
c                dvmxl = dmax1(dvmxl,abs(corr))
cc
cc                e(i,j,k) = e(i,j,k)+womega*r(i,j,k)
c            enddo
c          enddo
c         enddo
c         call imag8(e)
c         do i=0,idim1
c          do j=0,jdim1
c            e(i,j,0) = e(i,j,1)
c            e(i,j,kdim1) = e(i,j,kdim)
c          enddo
c         enddo
c        endif
c
         call MPI_ALLREDUCE(dvmxl,dvmx,1,
     &        MPI_REAL8,MPI_MAX,comm, nerr)
c        if(iope) write(*,*)"First dvmx is ",dvmx

c
c compute residual, assumes the forc is opposite sign
c
c      if(mod(it,5).eq.0) then
#if defined(USEMPI)
c         call MPI_ALLREDUCE(dvmxl,dvmx,1,
c     &        MPI_REAL8,MPI_MAX,comm, nerr)
#else
        dvmx = dvmxl
#endif
c        if(iope) write(*,*)"First dvmx is ",dvmx

c
        if(dvmx.le.eps) then
          iter = it 
          goto 100
        else
          if(iope) write(*,*)"dvmx ",it,dvmx
        endif
c      endif
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
            p(i,j,k) = e(i,j,k)
            pbarl = pbarl + p(i,j,k)
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

#if defined(MPOND) 
c             r(i,j,k) = (p(i+1,j,k)+p(i-1,j,k)+p(i,j,k+1)+p(i,j,k-1)-
c     *                  4.*p(i,j,k))*tdxs + (vay(i,j+1,k)*(p(i,j+1,k)-
c     *                  p(i,j,k))-vay(i,j,k)*(p(i,j,k)-p(i,j-1,k)))*tdxs

             r(i,j,k) = (ax(i+1,j,k)*(p(i+1,j,k)-p(i,j,k))-
     *                  ax(i,j,k)*(p(i,j,k)-p(i-1,j,k))
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
c        beta=0.055*(dxx**2+dy**2+dz**2)
        beta=0.025*(dxx**2+dy**2+dz**2)
c
c number of passes
c
        ntr= 1
  
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              ae(i,j,k)=0.
c              e(i,j,k) =  beta*r(i,j,k)
              e(i,j,k) =  0.5*r(i,j,k)/(2./dx(i)**2+4./(dy*dy))
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

                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
     *                  (dx(i)*dx(i))
     *                  +(e(i,j+1,k)+e(i,j-1,k)+
     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
c
c             ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)+e(i,j,k+1)+e(i,j,k-1)-
c     *                  4.*e(i,j,k))*tdxs + (vay(i,j+1,k)*(e(i,j+1,k)-
c     *                  e(i,j,k))-vay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs

#if defined(MPOND)
             ae(i,j,k) = (ax(i+1,j,k)*(e(i+1,j,k)-e(i,j,k))-
     *                  ax(i,j,k)*(e(i,j,k)-e(i-1,j,k))
     *                  +way(i,j,k+1)*(e(i,j,k+1)-e(i,j,k))-
     *                   way(i,j,k)*(e(i,j,k)-e(i,j,k-1))
     *                   +ay(i,j+1,k)*(e(i,j+1,k)-
     *                   e(i,j,k))-ay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs
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

                ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)-2.*e(i,j,k))/
     *                  (dx(i)*dx(i))
     *                  +(e(i,j+1,k)+e(i,j-1,k)+
     *                  e(i,j,k+1)+e(i,j,k-1)-4.*e(i,j,k))*tdxs
c
c             ae(i,j,k) = (e(i+1,j,k)+e(i-1,j,k)+e(i,j,k+1)+e(i,j,k-1)-
c     *                  4.*e(i,j,k))*tdxs + (vay(i,j+1,k)*(e(i,j+1,k)-
c     *                  e(i,j,k))-vay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs

#if defined(MPOND)
             ae(i,j,k) = (ax(i+1,j,k)*(e(i+1,j,k)-e(i,j,k))-
     *                  ax(i,j,k)*(e(i,j,k)-e(i-1,j,k))
     *                  +way(i,j,k+1)*(e(i,j,k+1)-e(i,j,k))-
     *                   way(i,j,k)*(e(i,j,k)-e(i,j,k-1))
     *                   +ay(i,j+1,k)*(e(i,j+1,k)-
     *                   e(i,j,k))-ay(i,j,k)*(e(i,j,k)-e(i,j-1,k)))*tdxs
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
      
