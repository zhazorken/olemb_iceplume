      subroutine cycelim(a,rr,lk,nmx,nm4,n2,m2)
c
c  $Id: cycelim.F,v 1.4 1998/10/16 21:21:20 dwd Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "communicate.inc"
c
      integer lk,nmx,nm4,n2,m2,k1,k2,i,j,k
      integer nnmx,nn,n
c
      real a(idim,jdim,kdimt2),b(idim,jdim,kdimt2)
     *,c(idim,jdim,kdimt2),a1(idim,jdim,kdimt2),a2(idim,jdim,kdimt2)
      complex
     * r(idim,jdim,kdimt2)
     *,rt(idim,jdim,kdimt2)
     *,rr(idim,jdim,kdim)
     *,frr(idim,jdim,kdim),xx(0:kdim1)
c
c
      nnmx=.5000001*lk
c
      k1=nnmx+1
      k2=nnmx+lk
c *******
c *** special sol for inverse of matrix with det=0 ***
c *** rhs=rr, rr(1,1,2) not even used but sum rr(1,1,:) should=0 ****
c
c for mpi code, we only want to do this if blockstart is 1
c
c this routine assumes idim and jdim are the sam
c rr is actually transposed i.e. rr(jdim,idim,kdim), but this doesn't
c matter unless idim .ne. jdim cause the a matrix is symmetric.
c
c will need to modify if idim.ne.jdim
c
c
      if(blockstart(my_pe+1).eq.1) then
        do k=1,lk+1
          xx(k)=0.0
        enddo
        do k=lk,3,-1
           xx(k-1)=-a(1,1,nnmx+k)*xx(k)-xx(k+1)+rr(1,1,k)
        enddo
        xx(1)=(rr(1,1,1)-xx(2))/a(1,1,nnmx+1)
        do k=1,lk
          rr(1,1,k)=xx(k)
        enddo
      endif
c *******
      do k=1,lk
        do j=1,m2
          do i=1,n2
            frr(i,j,k)=rr(i,j,k)
          enddo
        enddo
      enddo
c  **** frr now replaces rr *****
      do k=1,nm4
        do j=1,m2
          do i=1,n2
            r(i,j,k) = 0.0
            rt(i,j,k) = 0.0
            b(i,j,k) = 0.0
            c(i,j,k) = 0.0
            a1(i,j,k) = 0.0
            a2(i,j,k) = 0.0
          enddo
        enddo
      enddo
      do k=1,lk
        do j=1,m2
          do i=1,n2
            r(i,j,nnmx+k)=frr(i,j,k)
          enddo
        enddo
      enddo
      do k=k1,k2-1
        do j=1,m2
          do i=1,n2
            b(i,j,k)=1.
            c(i,j,k+1)=1.
          enddo
        enddo
      enddo
      nn=1
c
c here is the reduction loop
c
      do n=1,nmx
         do k=k1,k2
           do j=1,m2
             do i=1,n2
               a1(i,j,k)=c(i,j,k)/a(i,j,k-nn)
               a2(i,j,k)=b(i,j,k)/a(i,j,k+nn)
             enddo
           enddo
         enddo
         do k=k1,k2
           do j=1,m2
             do i=1,n2
               a(i,j,k)=a(i,j,k)
     *           -a1(i,j,k)*b(i,j,k-nn)
     *           -a2(i,j,k)*c(i,j,k+nn)
             enddo
           enddo
         enddo
         do k=k2,k1,-1
           do j=1,m2
             do i=1,n2
               c(i,j,k)=-a1(i,j,k)*c(i,j,k-nn)
             enddo
           enddo
         enddo
         do k=k1,k2
           do j=1,m2
             do i=1,n2
               b(i,j,k)=-a2(i,j,k)*b(i,j,k+nn)
             enddo
           enddo
         enddo
c
c  ***** rt needed cause of do loop ****
c
         do k=k1,k2
           do j=1,m2
             do i=1,n2
               rt(i,j,k)=r(i,j,k)
     *           -a1(i,j,k)*r(i,j,k-nn)
     *           -a2(i,j,k)*r(i,j,k+nn)
             enddo
           enddo
         enddo
         do k=k1,k2
           do j=1,m2
             do i=1,n2
               r(i,j,k)=rt(i,j,k)
             enddo
           enddo
         enddo
         nn=nn*2
      enddo
c
c done with the reduction
c
      do k=1,lk
        do j=1,m2
          do i=1,n2
            frr(i,j,k)=r(i,j,nnmx+k)/(a(i,j,nnmx+k)+1.e-20)
          enddo
        enddo
      enddo
c
c need to copy back frr into rr with special  care when dealing
c with the j index for MPI
c
c this is where the transposed fft routines have an effect
c
c
c if on edge of domain then need to be careful about corner point
c
      if(blockstart(my_pe+1).eq.1) then
        do k=1,lk
          do j=1,m2
            do i=2,n2
              rr(i,j,k)=frr(i,j,k)
            enddo
          enddo
        enddo
        do k=1,lk
          do j=2,m2
            rr(1,j,k)=frr(1,j,k)
          enddo
        enddo
      else
        do k=1,lk
          do j=1,m2
            do i=1,n2
              rr(i,j,k)=frr(i,j,k)
            enddo
          enddo
        enddo
      endif
c
c all done
c
      return
      end
