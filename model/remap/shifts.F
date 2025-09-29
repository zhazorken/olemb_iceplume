      subroutine northshift(vv,vvnorth)
c
c  $Id: shifts.F,v 1.5 1998/10/21 21:56:42 dwd Exp $
c
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real vv(0:idim1,0:jdim1,0:kdim1),vvnorth(0:idim1,0:kdim1)
      real buf(0:idim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i,k
c
c post receives
c
      call MPI_IRECV(vvnorth(0,0),1,nscont_type_real,nbr_north,
     &               mpitag_nshift,comm,request(2),ierr)

c
c pack data before sending
c
      do k=0,kdim1
        do i=0,idim1
          buf(i,k) = vv(i,2,k)
        enddo
      enddo
c
c post sends
c
      call MPI_ISEND(buf,1,nscont_type_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
      return
      end
      subroutine southshift(vv,vvsouth) 
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c 
      include 'mpif.h' 
c 
#include "param.inc" 
#include "communicate.inc" 
c 
      real vv(0:idim1,0:jdim1,0:kdim1),vvsouth(0:idim1,0:kdim1) 
      real buf(0:idim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i,k
c 
c post receives 
c 
      call MPI_IRECV(vvsouth(0,0),1,nscont_type_real,nbr_south, 
     &               mpitag_sshift,comm,request(2),ierr) 
 
c
c
c pack data before sending 
c 
      do k=0,kdim1 
        do i=0,idim1 
          buf(i,k) = vv(i,jdim-1,k) 
        enddo 
      enddo
c
c post sends
c
      call MPI_ISEND(buf,1,nscont_type_real,nbr_north,
     &               mpitag_sshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
c
      return
      end
c
c next do shifts for 2-d arrays
c
      subroutine northshift1(vv,vvnorth)
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real vv(0:idim1,0:jdim1),vvnorth(0:idim1)
      real buf(0:idim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i
c
c post receives
c
      call MPI_IRECV(vvnorth(0),idim2,MPI_REAL,nbr_north,
     &               mpitag_nshift,comm,request(2),ierr)

c
c pack data before sending
c
        do i=0,idim1
          buf(i) = vv(i,2)
        enddo
c
c post sends
c
      call MPI_ISEND(buf,idim2,MPI_REAL,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
      return
      end
      subroutine southshift1(vv,vvsouth)
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real vv(0:idim1,0:jdim1),vvsouth(0:idim1)
      real buf(0:idim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i
c
c post receives
c
      call MPI_IRECV(vvsouth(0),idim2,MPI_REAL,nbr_south,
     &               mpitag_sshift,comm,request(2),ierr)

c
c
c pack data before sending
c
        do i=0,idim1
          buf(i) = vv(i,jdim-1)
        enddo
c
c post sends
c
      call MPI_ISEND(buf,idim2,MPI_REAL,nbr_north,
     &               mpitag_sshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
c
      return
      end
c
#if !defined(SNGLSALT)
      subroutine northshiftdp(vv,vvnorth)
c
c  $Id: shifts.F,v 1.5 1998/10/21 21:56:42 dwd Exp $
c
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real*8 vv(0:idim1,0:jdim1,0:kdim1),vvnorth(0:idim1,0:kdim1)
      real*8 buf(0:idim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i,k
c
c post receives
c
      call MPI_IRECV(vvnorth(0,0),1,nscont_type_double,nbr_north,
     &               mpitag_nshift,comm,request(2),ierr)

c
c pack data before sending
c
      do k=0,kdim1
        do i=0,idim1
          buf(i,k) = vv(i,2,k)
        enddo
      enddo
c
c post sends
c
      call MPI_ISEND(buf,1,nscont_type_double,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
      return
      end
      subroutine southshiftdp(vv,vvsouth) 
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c 
      include 'mpif.h' 
c 
#include "param.inc" 
#include "communicate.inc" 
c 
      real*8 vv(0:idim1,0:jdim1,0:kdim1),vvsouth(0:idim1,0:kdim1) 
      real*8 buf(0:idim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i,k
c 
c post receives 
c 
      call MPI_IRECV(vvsouth(0,0),1,nscont_type_double,nbr_south, 
     &               mpitag_sshift,comm,request(2),ierr) 
 
c
c
c pack data before sending 
c 
      do k=0,kdim1 
        do i=0,idim1 
          buf(i,k) = vv(i,jdim-1,k) 
        enddo 
      enddo
c
c post sends
c
      call MPI_ISEND(buf,1,nscont_type_double,nbr_north,
     &               mpitag_sshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
c
      return
      end
c
c next do shifts for 2-d arrays
c
      subroutine northshift1dp(vv,vvnorth)
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real*8 vv(0:idim1,0:jdim1),vvnorth(0:idim1)
      real*8 buf(0:idim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i
c
c post receives
c
      call MPI_IRECV(vvnorth(0),idim2,MPI_DOUBLE_PRECISION,nbr_north,
     &               mpitag_nshift,comm,request(2),ierr)

c
c pack data before sending
c
        do i=0,idim1
          buf(i) = vv(i,2)
        enddo
c
c post sends
c
      call MPI_ISEND(buf,idim2,MPI_DOUBLE_PRECISION,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
      return
      end
      subroutine southshift1dp(vv,vvsouth)
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real*8 vv(0:idim1,0:jdim1),vvsouth(0:idim1)
      real*8 buf(0:idim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,i
c
c post receives
c
      call MPI_IRECV(vvsouth(0),idim2,MPI_DOUBLE_PRECISION,nbr_south,
     &               mpitag_sshift,comm,request(2),ierr)

c
c
c pack data before sending
c
        do i=0,idim1
          buf(i) = vv(i,jdim-1)
        enddo
c
c post sends
c
      call MPI_ISEND(buf,idim2,MPI_DOUBLE_PRECISION,nbr_north,
     &               mpitag_sshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
c
      return
      end
#endif
