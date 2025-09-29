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
      call MPI_IRECV(vvnorth(0,0),idim2*kdim2,MPI_REAL,nbr_north,
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
      call MPI_ISEND(buf,idim2*kdim2,MPI_REAL,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
#if defined(CHANNEL)
      if(jblkstart(my_pe+1)+jdim.gt.jtotal)then
        do k=0,kdim1
          do i=0,idim1
c            vvnorth(i,k) = vv(i,jdim,k)
            vvnorth(i,k) = 0.0
          enddo
        enddo
      endif
#endif

      return
      end
c
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
      call MPI_IRECV(vvsouth(0,0),idim2*kdim2,MPI_REAL,nbr_south,
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
      call MPI_ISEND(buf,idim2*kdim2,MPI_REAL,nbr_north,
     &               mpitag_sshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
c
#if defined(CHANNEL)
      if(jblkstart(my_pe+1).eq.1)then
        do k=0,kdim1
          do i=0,idim1
c            vvsouth(i,k) = vv(i,1,k)
            vvsouth(i,k) = 0.0
          enddo
        enddo
      endif
#endif

      return
      end
c
c do east west shifts
c

      subroutine eastshift(vv,vveast)
c
c  $Id: shifts.F,v 1.5 1998/10/21 21:56:42 dwd Exp $
c
c
c This routine returns an array of boundary points from the sub array
c that is east of the on-processor array in physical space.
c
      include 'mpif.h'
c
#include "param.inc"
#include "communicate.inc"
c
      real vv(0:idim1,0:jdim1,0:kdim1),vveast(0:jdim1,0:kdim1)
      real buf(0:jdim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,j,k

c
c post receives
c
      call MPI_IRECV(vveast(0,0),jdim2*kdim2,MPI_REAL,nbr_east,
     &               mpitag_eshift,comm,request(2),ierr)

c
c pack data before sending
c
      do k=0,kdim1
        do j=0,jdim1
          buf(j,k) = vv(2,j,k)
        enddo
      enddo
c
c post sends
c
      call MPI_ISEND(buf,jdim2*kdim2,MPI_REAL,nbr_west,
     &               mpitag_eshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
#if defined(OPENBC)
      if(iblkstart(my_pe+1)+idim.ge.itotal)then
        do k=0,kdim1
          do j=0,jdim1
            vveast(j,k) = vv(idim,j,k)
          enddo
        enddo
      endif
#endif
      return
      end
c
      subroutine westshift(vv,vvwest) 
c
c This routine returns an array of boundary points from the sub array
c that is north of the on-processor array in physical space.
c 
      include 'mpif.h' 
c 
#include "param.inc" 
#include "communicate.inc" 
c 
      real vv(0:idim1,0:jdim1,0:kdim1),vvwest(0:jdim1,0:kdim1) 
      real buf(0:jdim1,0:kdim1)
      integer request(2),status(MPI_STATUS_SIZE,2),ierr,j,k
c
c 
c post receives 
c 
      call MPI_IRECV(vvwest(0,0),jdim2*kdim2,MPI_REAL,nbr_west,
     &               mpitag_wshift,comm,request(2),ierr) 
 
c
c
c pack data before sending 
c 
      do k=0,kdim1 
        do j=0,jdim1 
          buf(j,k) = vv(idim-1,j,k) 
        enddo 
      enddo
c
c post sends
c
      call MPI_ISEND(buf,jdim2*kdim2,MPI_REAL,nbr_east,
     &               mpitag_wshift,comm,request(1),ierr)
c
      call MPI_WAITALL(2,request,status,ierr)
#if defined(OPENBC)
      if(iblkstart(my_pe+1).eq.1) then
        do k=0,kdim1
          do j=0,jdim1
            vvwest(j,k) = vv(1,j,k)
          enddo
        enddo
      endif
#endif
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
#if defined(CHANNEL)
      if(jblkstart(my_pe+1)+jdim.gt.jtotal)then
          do i=0,idim1
            vvnorth(i) = vv(i,jdim)
          enddo
      endif
#endif

      return
      end
c
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
#if defined(CHANNEL)
      if(jblkstart(my_pe+1).eq.1)then
          do i=0,idim1
            vvsouth(i) = vv(i,1)
          enddo
      endif
#endif

      return
      end
