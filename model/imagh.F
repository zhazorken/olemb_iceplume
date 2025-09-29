      subroutine imagh(phi)
c
c    $Id: imag.F,v 1.6 1998/10/21 21:56:38 dwd Exp $
c
c this routine sets the lateral image points
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "moddef.inc"
#include "communicate.inc"
      integer idimh,jdimh,kdimh,idimh1,jdimh1,kdimh1
      parameter(idimh=idim/2,jdimh=jdim/2,kdimh=kdim/2,
     *          idimh1=idimh+1,jdimh1=jdimh+1,kdimh1=kdimh+1)

c
      real phi(0:idimh1,0:jdimh1,0:kdimh1)
      real pyze(1:jdimh,0:kdimh1)
      real pyzein(1:jdimh,0:kdimh1)
      real pyzw(1:jdimh,0:kdimh1)
      real pyzwin(1:jdimh,0:kdimh1)
      real pxzn(0:idimh1,0:kdimh1)
      real pxznin(0:idimh1,0:kdimh1)
      real pxzs(0:idimh1,0:kdimh1)
      real pxzsin(0:idimh1,0:kdimh1)
      integer j,k
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
#else
      integer i
#endif
c
#if defined(USEMPI)
      do k=0,kdimh1
        do j=1,jdimh
c
c first start with the east-west exchange
c 
c left boundary image point
c 
         pyze(j,k) = phi(idimh,j,k)
c 
c right boundary image point
c
         pyzw(j,k) = phi(1,j,k)
c
        enddo
      enddo
      call MPI_IRECV(pyzwin(1,0),kdimh2*jdimh,MPI_REAL,nbr_west,
     &              mpitag_eshift,comm,request(3),ierr)
      call MPI_IRECV(pyzein(1,0),kdimh2*jdimh,MPI_REAL,nbr_east,
     &              mpitag_wshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pyzw(1,0),kdimh2*jdimh,MPI_REAL,nbr_west,
     &               mpitag_wshift,comm,request(1),ierr)
      call MPI_ISEND(pyze(1,0),kdimh2*jdimh,MPI_REAL,nbr_east,
     &               mpitag_eshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
      call MPI_BARRIER(comm,ierr)
c
c fill in image points
c
      do k=0,kdimh1
        do j=1,jdimh
c 
c east boundary image point
c 
         phi(idimh1,j,k)= pyzein(j,k)
c 
c west boundary image point
c
         phi(0,j,k) = pyzwin(j,k)
c
        enddo
      enddo
c
c next do the north-south exchange and include image points
c that were just updated
c
      do k=0,kdimh1
        do i=0,idimh1
c
c north boundary image point
c 
         pxzn(i,k) = phi(i,jdimh,k)
c 
c south boundary image point
c
         pxzs(i,k) = phi(i,1,k)
c
        enddo
      enddo
c
      call MPI_IRECV(pxznin(0,0),kdimh2*idimh2,MPI_REAL,nbr_north,
     &              mpitag_sshift,comm,request(3),ierr)
      call MPI_IRECV(pxzsin(0,0),kdimh2*idimh2,MPI_REAL,nbr_south,
     &              mpitag_nshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pxzn(0,0),kdimh2*idimh2,MPI_REAL,nbr_north,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(pxzs(0,0),kdimh2*idimh2,MPI_REAL,nbr_south,
     &               mpitag_sshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
      call MPI_BARRIER(comm,ierr)
c
c fill in image points
c
      do k=0,kdimh1
        do i=0,idimh1
c 
c east boundary image point
c 
         phi(i,jdimh1,k)= pxznin(i,k)
c 
c west boundary image point
c
         phi(i,0,k) = pxzsin(i,k)
c
        enddo
      enddo
c
#else
c
c the usual serial method
c
      do k=0,kdimh1
        do i=0,idimh1
          phi(i,0,k) = phi(i,jdimh,k)
          phi(i,jdimh1,k) = phi(i,1,k)
        enddo
      enddo
c
#endif
c
      return
      end
c
c
