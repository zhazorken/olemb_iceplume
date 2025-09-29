      subroutine imagd(phi)
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
      real phi(0:idimd1,0:jdimd1,0:kdimd1)
      real phit(0:idimd1,0:jdimd1,0:kdimd1)
      integer j,k,np,leftb,rightb,outb
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
      integer stat(MPI_STATUS_SIZE),req(2)
#endif
      integer i
c
      do k=0,kdimd1
        do j=0,jdimd1
c 
c left boundary image point
c 
         phi(0,j,k) = phi(idimd,j,k)
c 
c right boundary image point
c 
c
         phi(idimd1,j,k) = phi(1,j,k)
        enddo
      enddo

c next do the j direction image points
c
c MPI Code Section
c
#if defined(USEMPI)
#if defined(OPENBC)
c
c This series of receives and sends exchanges the ghost points 
c for the j coordinate boundaries. The communication is buffered
c by using a derived datatype that is defined in setup.F.  This
c datatype pulls out individual blocks of length 0:jdimd1 from the
c local 3-d array, which are separated by (jdimd+2)*(idimd+2)
c points.

c        write(*,*)"mark1 -imag2,relstep,my_pe",relstep,my_pe
c all proc

      do k=0,kdimd1
        do j=0,jdimd1
          do i=0,idimd1
            phit(i,j,k) = phi(i,j,k)
          enddo
        enddo
      enddo
c
c post receives 
c
      call MPI_IRECV(phit(0,jdimd1,0),1,ns_type_real,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phit(0,0,0),1,ns_type_real,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phit(0,1,0),1,ns_type_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phit(0,jdimd,0),1,ns_type_real,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c
c
c if open b.c., then don't set the real 0 and jtotal points
c
      if((blockstart(my_pe+1)+jdimd.lt.jtotal)
     &    .and.(blockstart(my_pe+1).gt.1)) then

c        write(*,*)"mark2 -imag2,relstep,my_pe",relstep,my_pe
c  proc 1 and 2 only

        do k=0,kdimd1
          do i=0,idimd1
             phi(i,jdimd1,k) = phit(i,jdimd1,k)
             phi(i,0,k) = phit(i,0,k)
          enddo
        enddo
      endif
c
c don't do the top and bottom b.c. image points
c
      if((blockstart(my_pe+1)+jdimd).ge.jtotal) then
c        write(*,*)"mark3 -imag2,relstep,my_pe",relstep,my_pe
c proc 3 only
        do k=0,kdimd1
          do i=0,idimd1
             phi(i,0,k) = phit(i,0,k)
             phi(i,jdimd1,k) = phi(i,jdimd,k)
          enddo
        enddo
      endif
c
#if defined(RECIRC)
c
c post receive if on blockstart 1
c
      do np=1,nprocs
        if(blockstart(np).eq.1) then
          leftb = np-1
        else if(blockstart(np).eq.((np_recirc-1)*jdimd+1)) then
c       else if(blockstart(np).eq.( 2           *jdimd+1)) then
          rightb = np-1
        else if(blockstart(np).eq.((np_recirc  )*jdimd+1)) then
c       else if(blockstart(np).eq.( 3           *jdimd+1)) then
          outb = np-1
        endif
      enddo

c         write(*,*)"mark6a-imag2,relstep,np,leftb ",relstep,np,leftb
c         write(*,*)"mark6b-imag2,relstep,np,rightb",relstep,np,rightb
c         write(*,*)"mark6c-imag2,relstep,np,outb  ",relstep,np,outb
c mark csmith - leftb = 0; rightb = 2; outb = 3;

      if(blockstart(my_pe+1).eq.1) then
c        write(*,*)"mark7 -imag2,relstep,my_pe",relstep,my_pe
c proc 0 only
        call MPI_RECV(phi(0,0,0),1,ns_type_real,rightb,
     &              mpitag_sshift,comm,stat,ierr)
c        call MPI_RECV(phi(0,1,0),1,ns_type_real,rightb,
c     &              mpitag_nshift,comm,stat,ierr)


        do k=0,kdimd1
          do i=0,idimd1
           phi(i,jdimd1,k) = phit(i,jdimd1,k)
          enddo
        enddo

      else if(blockstart(my_pe+1).eq.((np_recirc-1)*jdimd+1)) then
c     else if(blockstart(my_pe+1).eq.( 2           *jdimd+1)) then
c       write(*,*)"mark8 -imag2,relstep,my_pe",relstep,my_pe
c proc 2 only
        call MPI_SEND(phi(0,jdimd,0),1,ns_type_real,leftb,
     &               mpitag_sshift,comm,ierr)
c        call MPI_SEND(phi(0,jdimd1,0),1,ns_type_real,leftb,
c     &               mpitag_nshift,comm,ierr)

      endif
c
c send other cyclic point for internal boundary
c
      if(blockstart(my_pe+1).eq.((np_recirc-1)*jdimd+1)) then
c     if(blockstart(my_pe+1).eq.( 2           *jdimd+1)) then
c       write(*,*)"mark9 -imag2,relstep,my_pe",relstep,my_pe
c proc 2 only
        call MPI_RECV(phi(0,jdimd1,0),1,ns_type_real,leftb,
     &               mpitag_nshift,comm,stat,ierr)
      else if(blockstart(my_pe+1).eq. 1) then
c proc 0 only
        call MPI_SEND(phi(0,1,0),1,ns_type_real,rightb,
     &              mpitag_nshift,comm,ierr)

      endif
      call MPI_BARRIER(comm,ierr)
c
c finally, update neighber to internal boundary
c
      if(blockstart(my_pe+1).eq.((np_recirc  )*jdimd+1)) then
c     if(blockstart(my_pe+1).eq.( 3           *jdimd+1)) then
c proc 3 only
c       write(*,*)"mark10-imag2,relstep,my_pe",relstep,my_pe



c j=jdimd  proc 2 --> j=0 proc3
        call MPI_RECV(phi(0,0,0),1,ns_type_real,rightb,
     &               mpitag_nshift,comm,stat,ierr)
c j=jdimd1 proc 2 --> j=1 proc3
        call MPI_RECV(phi(0,1,0),1,ns_type_real,rightb,
     &               mpitag_sshift,comm,stat,ierr)
c
      else if(blockstart(my_pe+1).eq.((np_recirc-1)*jdimd+1)) then
c     else if(blockstart(my_pe+1).eq.( 2           *jdimd+1)) then
c proc 2 only
c       write(*,*)"mark8-imag2,relstep,my_pe",relstep,my_pe


c j=jdimd  proc 2 --> j=0 proc3
        call MPI_SEND(phi(0,jdimd,0),1,ns_type_real,outb,
     &               mpitag_nshift,comm,ierr)
c j=jdimd1 proc 2 --> j=1 proc3
        call MPI_SEND(phi(0,jdimd1,0),1,ns_type_real,outb,
     &               mpitag_sshift,comm,ierr)




      endif
      call MPI_BARRIER(comm,ierr)
c
#else
      if((blockstart(my_pe+1)).eq.1) then
c       write(*,*)"mark11-imag2,relstep,my_pe",relstep,my_pe
c proc not used
        do k=0,kdimd1
          do i=0,idimd1
             phi(i,jdimd1,k) = phit(i,jdimd1,k)
             phi(i,0,k) = phi(i,1,k)
          enddo
        enddo
      endif
c
      call MPI_BARRIER(comm,ierr)
#endif
c
#else
c       write(*,*)"mark12-imag2,relstep,my_pe",relstep,my_pe
c proc not used

c not open bc, do old periodic b.c. transfer
c post receives
c
      call MPI_IRECV(phi(0,jdimd1,0),1,ns_type_real,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phi(0,0,0),1,ns_type_real,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phi(0,1,0),1,ns_type_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phi(0,jdimd,0),1,ns_type_real,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c

      call MPI_BARRIER(comm,ierr)
#endif
c
#else
c       write(*,*)"mark13-imag2,relstep,my_pe",relstep,my_pe
c proc not used

c not mpi
c the usual serial method
c
#if defined(OPENBC)
c       write(*,*)"mark14-imag2,relstep,my_pe",relstep,my_pe
c proc not used
      do k=0,kdimd1
        do i=0,idimd1
          phi(i,0,k) = phi(i,1,k)
          phi(i,jdimd1,k) = phi(i,jdimd,k)
        enddo
      enddo
#else
c       write(*,*)"mark15-imag2,relstep,my_pe",relstep,my_pe
c proc not used
      do k=0,kdimd1
        do i=0,idimd1
          phi(i,0,k) = phi(i,jdimd,k)
          phi(i,jdimd1,k) = phi(i,1,k)
        enddo
      enddo
#endif
c
#endif
c
      return
      end
c
c
      subroutine imagdpd(phi)
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
      real*8 phi(0:idimd1,0:jdimd1,0:kdimd1)
      integer j,k
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
#endif
      integer i
c
      do k=0,kdimd1
        do j=0,jdimd1
c
c left boundary image point
c
         phi(0,j,k) = phi(idimd,j,k)
c
c right boundary image point
c
         phi(idimd1,j,k) = phi(1,j,k)
        enddo
      enddo
c
c next do the j direction image points
c
c MPI Code Section
c
#if defined(USEMPI)
c
c This series of receives and sends exchanges the ghost points 
c for the j coordinate boundaries. The communication is buffered
c by using a derived datatype that is defined in setup.F.  This
c datatype pulls out individual blocks of length 0:jdimd1 from the
c local 3-d array, which are separated by (jdimd+2)*(idimd+2)
c points.
c
c post receives 
c
      call MPI_IRECV(phi(0,jdimd1,0),1,ns_type_double,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phi(0,0,0),1,ns_type_double,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phi(0,1,0),1,ns_type_double,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phi(0,jdimd,0),1,ns_type_double,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c
#else
c
c the usual serial method
c
      do k=0,kdimd1
        do i=0,idimd1
          phi(i,0,k) = phi(i,jdimd,k)
          phi(i,jdimd1,k) = phi(i,1,k)
        enddo
      enddo
c
#endif
c
      return

      end
      subroutine imag2dd(phi)
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
      real phi(0:idimd1,0:jdimd1)
      real phit(0:idimd1,0:jdimd1)
      integer j,rightb,leftb,np,outb
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
      integer stat(MPI_STATUS_SIZE),req(2)

#endif
      integer i
c
      do j=0,jdimd1
c
c left boundary image point
c
        phi(0,j) = phi(idimd,j)
c
c right boundary image point
c
c
        phi(idimd1,j) = phi(1,j)
      enddo
c next do the j direction image points
c
c MPI Code Section
c
#if defined(USEMPI)
#if defined(OPENBC)
      do j=0,jdimd1
        do i=0,idimd1
          phit(i,j) = phi(i,j)
        enddo
      enddo
c This series of receives and sends exchanges the ghost points
c for the j coordinate boundaries. The communication is buffered
c by using a derived datatype that is defined in setup.F.  This
c datatype pulls out individual blocks of length 0:jdimd1 from the
c local 3-d array, which are separated by (jdimd+2)*(idimd+2)
c points.
c
c post receives
c
      call MPI_IRECV(phit(0,jdimd1),idimd+2,mpi_real,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phit(0,0),idimd+2,mpi_real,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phit(0,1),idimd+2,mpi_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phit(0,jdimd),idimd+2,mpi_real,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c
c open b.c., don't set the real 0 and jtotal points
c
      if((blockstart(my_pe+1)+jdimd.lt.jtotal)
     &    .and.(blockstart(my_pe+1).gt.1)) then
        do i=0,idimd1
           phi(i,jdimd1) = phit(i,jdimd1)
           phi(i,0) = phit(i,0)
        enddo
      endif
c
c don't do the top and bottom b.c. image points
c
      if((blockstart(my_pe+1)+jdimd).ge.jtotal) then
          do i=0,idimd1
             phi(i,0) = phit(i,0)
             phi(i,jdimd1) = phi(i,jdimd)
          enddo
      endif

#if defined(RECIRC)
c
c post receive if on blockstart 1
c
      do np=1,nprocs
        if(blockstart(np).eq.1) then
          leftb = np-1
        else if(blockstart(np).eq.2*jdimd+1) then
          rightb = np-1
        else if(blockstart(np).eq.(3*jdimd+1)) then
          outb = np-1
        endif
      enddo
      if(blockstart(my_pe+1).eq.1) then
c
        call MPI_RECV(phi(0,0),idimd+2,mpi_real,rightb,
     &              mpitag_sshift,comm,stat,ierr)
c
        do i=0,idimd1
          phi(i,jdimd1) = phit(i,jdimd1)
        enddo
c
      else if(blockstart(my_pe+1).eq.2*jdimd+1) then
c
        call MPI_SEND(phi(0,jdimd),idimd+2,mpi_real,leftb,
     &               mpitag_sshift,comm,ierr)
      endif
c
c send other cyclic point for internal boundary
c
      if(blockstart(my_pe+1).eq.(2*jdimd+1)) then
        call MPI_RECV(phi(0,jdimd1),idimd+2,mpi_real,leftb,
     &               mpitag_nshift,comm,stat,ierr)
      else if(blockstart(my_pe+1).eq. 1) then
        call MPI_SEND(phi(0,1),idimd+2,mpi_real,rightb,
     &              mpitag_nshift,comm,ierr)

      endif

      call MPI_BARRIER(comm,ierr)
c
c finally, update neighber to internal boundary
c
      if(blockstart(my_pe+1).eq.(3*jdimd+1)) then
        call MPI_RECV(phi(0,0),idimd+2,mpi_real,rightb,
     &               mpitag_nshift,comm,stat,ierr)
        call MPI_RECV(phi(0,1),idimd+2,mpi_real,rightb,
     &               mpitag_sshift,comm,stat,ierr)
c
      else if(blockstart(my_pe+1).eq.(2*jdimd+1)) then
c
        call MPI_SEND(phi(0,jdimd),idimd+2,mpi_real,outb,
     &               mpitag_nshift,comm,ierr)
        call MPI_SEND(phi(0,jdimd1),idimd+2,mpi_real,outb,
     &               mpitag_sshift,comm,ierr)
      endif
      call MPI_BARRIER(comm,ierr)

c
#else
      if((blockstart(my_pe+1)).eq.1) then
        do i=0,idimd1
           phi(i,jdimd1) = phit(i,jdimd1)
           phi(i,0) = phi(i,1)
        enddo
      endif
c
      call MPI_BARRIER(comm,ierr)
#endif

c

#else
c
c This series of receives and sends exchanges the ghost points
c for the j coordinate boundaries. The communication is buffered
c by using a derived datatype that is defined in setup.F.  This
c datatype pulls out individual blocks of length 0:jdimd1 from the
c local 3-d array, which are separated by (jdimd+2)*(idimd+2)
c points.
c
c post receives
c
      call MPI_IRECV(phi(0,jdimd1),idimd+2,mpi_real,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phi(0,0),idimd+2,mpi_real,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phi(0,1),idimd+2,mpi_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phi(0,jdimd),idimd+2,mpi_real,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c
#endif
#else
c
c the usual serial method
c
#if defined(OPENBC)
      do i=0,idimd1
        phi(i,0) = phi(i,1)
        phi(i,jdimd1) = phi(i,jdimd)
      enddo
#else
      do i=0,idimd1
        phi(i,0) = phi(i,jdimd)
        phi(i,jdimd1) = phi(i,1)
      enddo
#endif
c
#endif
c
      return
      end
c
      subroutine imag2did(phi)
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
      real phi(0:idimd1,0:jdimd1)
      real phit(0:idimd1,0:jdimd1)
      integer j,rightb,leftb,np,outb
#if defined(USEMPI)
      integer ierr,request(4),status(MPI_STATUS_SIZE,4)
      integer stat(MPI_STATUS_SIZE),req(2)

#endif
      integer i
c
      do j=0,jdimd1
c
c left boundary image point
c
        phi(0,j) = phi(idimd,j)
c
c right boundary image point
c
c
        phi(idimd1,j) = phi(1,j)
      enddo
c next do the j direction image points
c
c MPI Code Section
c
#if defined(USEMPI)
c
c This series of receives and sends exchanges the ghost points
c for the j coordinate boundaries. The communication is buffered
c by using a derived datatype that is defined in setup.F.  This
c datatype pulls out individual blocks of length 0:jdimd1 from the
c local 3-d array, which are separated by (jdimd+2)*(idimd+2)
c points.
c
c post receives
c
      call MPI_IRECV(phi(0,jdimd1),idimd+2,mpi_real,nbr_north,
     &              mpitag_nshift,comm,request(3),ierr)
      call MPI_IRECV(phi(0,0),idimd+2,mpi_real,nbr_south,
     &              mpitag_sshift,comm,request(4),ierr)
c
c post sends
c
      call MPI_ISEND(phi(0,1),idimd+2,mpi_real,nbr_south,
     &               mpitag_nshift,comm,request(1),ierr)
      call MPI_ISEND(phi(0,jdimd),idimd+2,mpi_real,nbr_north,
     &               mpitag_sshift,comm,request(2),ierr)
c
c wait till everyone is done
c
      call MPI_WAITALL(4,request,status,ierr)
c
#else
c
c the usual serial method
c
#if defined(OPENBC)
      do i=0,idimd1
        phi(i,0) = phi(i,1)
        phi(i,jdimd1) = phi(i,jdimd)
      enddo
#else
      do i=0,idimd1
        phi(i,0) = phi(i,jdimd)
        phi(i,jdimd1) = phi(i,1)
      enddo
#endif
c
#endif
c
      return
      end
c
