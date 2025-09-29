      subroutine broadcast_string(string, len,root_pe)
c
c  $Id: global.F,v 1.3 1998/10/16 21:21:21 dwd Exp $
c

      implicit none

      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
      integer root_pe, ierr, len, ilen, lnblk, i
      integer istring(80)
      character string*80

      if(iope) then
        ilen = lnblk(string,80)
        do i=1,ilen
          istring(i) = ichar(string(i:i))
        enddo
      endif
      call MPI_BCAST(ilen,1,MPI_INTEGER,mastertask,
     &     comm,ierr)
      call MPI_BCAST(istring, ilen, MPI_INTEGER, mastertask,
     &                       comm, ierr)
      do i=1,ilen
         string(i:i) = char(istring(i))
      enddo
      call MPI_BARRIER(comm, ierr)

      end

c
      subroutine broadcast_scalar_int(scalar, root_pe)

      implicit none

      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"

      integer root_pe, scalar, ierr

      call MPI_BCAST(scalar, 1, MPI_INTEGER, root_pe,
     &                       comm, ierr)
      call MPI_BARRIER(comm, ierr)

      end
c
