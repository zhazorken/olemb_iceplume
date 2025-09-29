      subroutine intg(var,sum)
c
c   $Id: intg.F,v 1.4 1998/10/21 21:56:38 dwd Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "communicate.inc"
c
      integer i,j,k
#if defined(USEMPI)
      integer ierr
#endif
      real var(0:idim1,0:jdim1,0:kdim1),sum
      real*8 sum8,totsum
c
      sum8=0.0d0
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               sum8=sum8+var(i,j,k)
            enddo
         enddo
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(sum8,totsum,1,     
     &     MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &     comm, ierr)
#else
      totsum = sum8
#endif
      if(iope) then
         sum=totsum/dble(idim*kdim*jtotal)
      endif
      return
      end
#if !defined(SNGLSALT)
      subroutine intgdp(var,sum)
c
c   $Id: intg.F,v 1.4 1998/10/21 21:56:38 dwd Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "communicate.inc"
c
      integer i,j,k
#if defined(USEMPI)
      integer ierr
#endif
      real sum
      real*8 var(0:idim1,0:jdim1,0:kdim1),sum8,totsum
c
      sum8=0.0d0
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               sum8=sum8+var(i,j,k)
            enddo
         enddo
      enddo
#if defined(USEMPI) 
      call MPI_REDUCE(sum8,totsum,1,     
     &     MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &     comm, ierr)
#else 
      totsum = sum8 
#endif 
      if(iope) then
        sum=totsum/dble(idim*kdim*jtotal)
      endif
      return
      end
#endif
