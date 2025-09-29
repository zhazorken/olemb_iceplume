      subroutine read3ds(cdfin,cdfid,fld)
c
c routine for reading in a subdomain field
c does not require MPI
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      real fld(0:idim1,0:jdim1,0:kdim1)
      real phi(idim,jdim,kdim)
      integer cdfin,cdfid,start(4),count(4),rcode


c
      count(1) = idim
      count(2) = jdim
      count(3) = kdim
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1

c
c      start(1) = iblkstart(my_pe+1)
c      start(2) = jblkstart(my_pe+1)
      start(1) = 1
      start(2) = 1
c      write(*,*)"reading input file on ",my_pe,start(1),start(2)
      call ncvgt(cdfin,cdfid,start,count,phi,rcode)
c      write(*,*)"finished reading on ",my_pe
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            fld(i,j,k) = phi(i,j,k)
          enddo
        enddo
      enddo
      call MPI_BARRIER(comm,ier)
      return
      end
