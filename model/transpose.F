      subroutine transpose(xi,xo,comm,my_pe)
c
c  $Id: transpose.F,v 1.3 1998/10/16 21:21:26 dwd Exp $
c
c
c this routine transposes array xi(idim,jtotal/nprocs) to array
c xo(jtotal,idim/nprocs) both internally and across the nodes
c
      include 'mpif.h'
#include "param.inc"

      integer idiml
      parameter(idiml = idim/nprocs)
c
      complex xi(idim,jdim,kdim), xo(jtotal,idim/nprocs,kdim)
      complex buf(jtotal*idiml*kdim)
      complex rbuf(jtotal*idiml*kdim)
      integer ibsiz, offset(nprocs),scount(nprocs)
      integer comm,totc,my_pe
      integer i,j,k,ierr,n

c
c first break up into buffers and do internal transpose
c
      ibsiz = idim/nprocs
c
      totc = 0.0
      do n=1,nprocs
        do k=1,kdim
        do i=1,ibsiz
          do j=1,jdim
            totc = totc+1
            buf(totc) = xi(i+ibsiz*(n-1),j,k)
          enddo
        enddo
        enddo
      
      enddo
      
c      if(my_pe.eq.0) then
c         write(*,*)
c         write(*,100)(buf(i),i=1,idim*jdim)
c         write(*,*)
c100      format(16f4.0)
c      endif
c
c next set up all to all communication routine
c
      do n=1,nprocs
        offset(n) = (n-1)*jdim*ibsiz*kdim
        scount(n) = ibsiz*jdim*kdim
      enddo
c
      call MPI_ALLTOALLV(buf,scount,offset,MPI_COMPLEX,
     $     rbuf,scount,offset,MPI_COMPLEX,comm,ierr)
c
c
c      if(my_pe.eq.0) then
c         write(*,*)
c         write(*,100)(rbuf(i),i=1,idim*jdim)
c         write(*,*)
c      endif
c
c next, unpack the buffer
c
      totc = 0
      do n=1,nprocs
        do k=1,kdim
        do j=1,ibsiz
          do i=1,jdim
            totc = totc + 1
            xo(i+jdim*(n-1),j,k) = rbuf(totc)
          enddo
        enddo
        enddo
      enddo
c
c all done
c
      return
      end
c
c routine for a single level transpose
c

      subroutine transpose1(xi,xo,comm,my_pe)
c
c this routine transposes array xi(idim,jtotal/nprocs) to array
c xo(jtotal,idim/nprocs) both internally and across the nodes
c
      include 'mpif.h'
#include "param.inc"

      integer idiml
      parameter(idiml = idim/nprocs)
c
      complex xi(idim,jdim), xo(jtotal,idim/nprocs)
      complex buf(jtotal*idiml)
      complex rbuf(jtotal*idiml)
      integer ibsiz, offset(nprocs),scount(nprocs)
      integer comm,totc,my_pe
      integer i,j,ierr,n

c
c first break up into buffers and do internal transpose
c
      ibsiz = idim/nprocs
c
      totc = 0.0
      do n=1,nprocs
        do i=1,ibsiz
          do j=1,jdim
            totc = totc+1
            buf(totc) = xi(i+ibsiz*(n-1),j)
          enddo
        enddo
      enddo
      
c      if(my_pe.eq.0) then
c         write(*,*)
c         write(*,100)(buf(i),i=1,idim*jdim)
c         write(*,*)
c100      format(16f4.0)
c      endif
c
c next set up all to all communication routine
c
      do n=1,nprocs
        offset(n) = (n-1)*jdim*ibsiz
        scount(n) = ibsiz*jdim
      enddo
c
      call MPI_ALLTOALLV(buf,scount,offset,MPI_COMPLEX,
     $     rbuf,scount,offset,MPI_COMPLEX,comm,ierr)
c
c
c      if(my_pe.eq.0) then
c         write(*,*)
c         write(*,100)(rbuf(i),i=1,idim*jdim)
c         write(*,*)
c      endif
c
c next, unpack the buffer
c
      totc = 0
      do n=1,nprocs
        do j=1,ibsiz
          do i=1,jdim
            totc = totc + 1
            xo(i+jdim*(n-1),j) = rbuf(totc)
          enddo
        enddo
      enddo
c
c all done
c
      return
      end
