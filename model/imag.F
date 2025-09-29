      subroutine imag(phi)
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
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real pyze(0:jdim1,0:kdim1)
      real pyzein(0:jdim1,0:kdim1)
      real pyzw(0:jdim1,0:kdim1)
      real pyzwin(0:jdim1,0:kdim1)
      real pxzn(0:idim1,0:kdim1)
      real pxznin(0:idim1,0:kdim1)
      real pxzs(0:idim1,0:kdim1)
      real pxzsin(0:idim1,0:kdim1)
      integer j,k
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
#else
      integer i
#endif
c
#if defined(USEMPI)
      do k=0,kdim1
        do j=0,jdim1
c
c first start with the east-west exchange
c 
c left boundary image point
c 
         pyze(j,k) = phi(idim,j,k)
c 
c right boundary image point
c
         pyzw(j,k) = phi(1,j,k)
c
        enddo
      enddo
      call MPI_IRECV(pyzwin(0,0),kdim2*jdim2,MPI_REAL,nbr_west,
     &              mpitag_eshift,comm,request(3),ierr)
      call MPI_IRECV(pyzein(0,0),kdim2*jdim2,MPI_REAL,nbr_east,
     &              mpitag_wshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pyzw(0,0),kdim2*jdim2,MPI_REAL,nbr_west,
     &               mpitag_wshift,comm,request(1),ierr)
      call MPI_ISEND(pyze(0,0),kdim2*jdim2,MPI_REAL,nbr_east,
     &               mpitag_eshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
      call MPI_BARRIER(comm,ierr)
c
c fill in image points
c

#if defined(OPENBC)
      do k=0,kdim1
        do j=0,jdim1
c 
c east boundary image point
c 
         if(iblkstart(my_pe+1)+idim.ge.itotal) then
           phi(idim1,j,k)= phi(idim,j,k)
         else
           phi(idim1,j,k)= pyzein(j,k)
         endif
c 
c west boundary image point
c
         if(iblkstart(my_pe+1).eq.1) then
           phi(0,j,k) = phi(1,j,k)
         else
           phi(0,j,k) = pyzwin(j,k)
         endif
c
        enddo
      enddo
#else
      do k=0,kdim1
        do j=0,jdim1
c 
c east boundary image point
c 
         phi(idim1,j,k)= pyzein(j,k)
c 
c west boundary image point
c
         phi(0,j,k) = pyzwin(j,k)
c
        enddo
      enddo

#endif
c
c next do the north-south exchange and include image points
c that were just updated
c
      do k=0,kdim1
        do i=0,idim1
c
c north boundary image point
c 
         pxzn(i,k) = phi(i,jdim,k)
c 
c south boundary image point
c
         pxzs(i,k) = phi(i,1,k)
c
        enddo
      enddo
c
      call MPI_IRECV(pxznin(0,0),kdim2*idim2,MPI_REAL,nbr_north,
     &              mpitag_sshift,comm,request(3),ierr)
      call MPI_IRECV(pxzsin(0,0),kdim2*idim2,MPI_REAL,nbr_south,
     &              mpitag_nshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pxzn(0,0),kdim2*idim2,MPI_REAL,nbr_north,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(pxzs(0,0),kdim2*idim2,MPI_REAL,nbr_south,
     &               mpitag_sshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
      call MPI_BARRIER(comm,ierr)
c
c fill in image points
c
#if defined(CHANNEL)
      do k=0,kdim1
        do i=0,idim1
c 
c east boundary image point
c 
         if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
           phi(i,jdim1,k)= phi(i,jdim,k)
         else
           phi(i,jdim1,k)= pxznin(i,k)
         endif
c 
c west boundary image point
c
         if(jblkstart(my_pe+1).eq.1) then
           phi(i,0,k) = phi(i,1,k)
         else
           phi(i,0,k) = pxzsin(i,k)
         endif
c
        enddo
      enddo

#else

      do k=0,kdim1
        do i=0,idim1
c 
c east boundary image point
c 
         phi(i,jdim1,k)= pxznin(i,k)
c 
c west boundary image point
c
         phi(i,0,k) = pxzsin(i,k)
c
        enddo
      enddo
#endif
c
#else
c
c the usual serial method
c
      do k=0,kdim1
        do i=0,idim1
          phi(i,0,k) = phi(i,jdim,k)
          phi(i,jdim1,k) = phi(i,1,k)
        enddo
      enddo
c
#endif
c
      return
      end
c
c
      subroutine imag2d(phi)
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
c
      real phi(0:idim1,0:jdim1)
      real pyze(1:jdim1)
      real pyzein(1:jdim1)
      real pyzw(1:jdim1)
      real pyzwin(1:jdim1)
      real pxzn(0:idim1)
      real pxznin(0:idim1)
      real pxzs(0:idim1)
      real pxzsin(0:idim1)
      integer j,k
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
#else
      integer i
#endif
c
#if defined(USEMPI)
      do j=1,jdim
c
c first start with the east-west exchange
c 
c left boundary image point
c 
         pyze(j) = phi(idim,j)
c 
c right boundary image point
c
         pyzw(j) = phi(1,j)
c
      enddo
      call MPI_IRECV(pyzwin,jdim,MPI_REAL,nbr_west,
     &              mpitag_eshift,comm,request(3),ierr)
      call MPI_IRECV(pyzein,jdim,MPI_REAL,nbr_east,
     &              mpitag_wshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pyzw,jdim,MPI_REAL,nbr_west,
     &               mpitag_wshift,comm,request(1),ierr)
      call MPI_ISEND(pyze,jdim,MPI_REAL,nbr_east,
     &               mpitag_eshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
c
c fill in image points
c
#if defined(OPENBC)
      do j=1,jdim
c 
c east boundary image point
c 
         if(iblkstart(my_pe+1)+idim.ge.itotal) then
           phi(idim1,j)= phi(idim,j)
         else
           phi(idim1,j)= pyzein(j)
         endif
c 
c west boundary image point
c
         if(iblkstart(my_pe+1).eq.1) then
           phi(0,j) = phi(1,j)
         else
           phi(0,j) = pyzwin(j)
         endif
c
      enddo

#else
      do j=1,jdim
c 
c east boundary image point
c 
         phi(idim1,j)= pyzein(j)
c 
c west boundary image point
c
         phi(0,j) = pyzwin(j)
c
      enddo
#endif
c
c next do the north-south exchange and include image points
c that were just updated
c
      do i=0,idim1
c
c north boundary image point
c 
         pxzn(i) = phi(i,jdim)
c 
c south boundary image point
c
         pxzs(i) = phi(i,1)
c
      enddo
c
      call MPI_IRECV(pxznin,idim2,MPI_REAL,nbr_north,
     &              mpitag_sshift,comm,request(3),ierr)
      call MPI_IRECV(pxzsin,idim2,MPI_REAL,nbr_south,
     &              mpitag_nshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(pxzn,idim2,MPI_REAL,nbr_north,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(pxzs,idim2,MPI_REAL,nbr_south,
     &               mpitag_sshift,comm,request(2),ierr)

      call MPI_WAITALL(4,request,status,ierr)
c
c fill in image points
c
#if defined(CHANNEL)
        do i=0,idim1
c 
c east boundary image point
c 
         if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
           phi(i,jdim1)= phi(i,jdim)
         else
           phi(i,jdim1)= pxznin(i)
         endif
c 
c west boundary image point
c
         if(jblkstart(my_pe+1).eq.1) then
           phi(i,0) = phi(i,1)
         else
           phi(i,0) = pxzsin(i)
         endif
c
        enddo

#else

      do i=0,idim1
c 
c east boundary image point
c 
         phi(i,jdim1)= pxznin(i)
c 
c west boundary image point
c
         phi(i,0) = pxzsin(i)
c
      enddo
#endif
c
#else
c
c the usual serial method
c
      do i=0,idim1
          phi(i,0) = phi(i,jdim)
          phi(i,jdim1) = phi(i,1)
      enddo
c
#endif
c
      return
      end
