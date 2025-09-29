c
#if !defined(SNGLSALT)
       subroutine collect_3ddp(temp,phi)
c
c  $Id: iompi.F,v 1.8 1998/12/08 00:47:40 mikeh Exp $
c

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a full 3d field, including the image points.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req,k
      integer i
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real*8 temp(0:idim1,0:jdim1,0:kdim1)
      real*8 phi(0:idim1,0:jtotal1,0:kdim1)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c
c
      if (iope) then
        do np=1,nprocs
c
c post receives first
c
          call MPI_IRECV(phi(0,blockstart(np),0), 1,mpi_global_double,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
        enddo
      endif
      call MPI_ISEND(temp(0,1,0), 1,mpi_local_double, mastertask,
     &             mpitag_io, comm, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      call MPI_BARRIER(comm,ier)
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c need to do j image points
c
        do k=0,kdim1
          do i=0,idim1
            phi(i,0,k) = phi(i,jtotal,k)
            phi(i,jtotal1,k) = phi(i,1,k)
          enddo
        enddo
      endif

      return
      end

      subroutine distrib_3ddp(temp,phi)

      implicit none

      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"

c-----------------------------------------------------------------------
c
c This routine takes a full 3-d field (phi) and ships it out to each of the
c processors through "temp"
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req, send_req
      integer status(MPI_STATUS_SIZE)

      real*8 temp(0:idim1,0:jdim1,0:kdim1)
      real*8 phi(0:idim1,0:jtotal1,0:kdim1)
c-----------------------------------------------------------------------
c
c post recieves for each of the data chunks
c
c-----------------------------------------------------------------------

      call MPI_IRECV(temp(0,1,0), 1, mpi_local_double, mastertask,
     &              mpitag_io, comm, recv_req, ier)
c
c now send out parts of global variable to corresponding processors
c based on global address space defined in setupmpi routine
c
      if (iope) then
        do np=1,nprocs
c
c post the send for this processors data
c
          call MPI_ISEND(phi(0,blockstart(np),0), 1, 
     &                  mpi_global_double, np-1,
     &                  mpitag_io, comm, send_req, ier)
          call MPI_WAIT(send_req, status, ier)
        enddo
      endif
c
c wait for everything to communicate
c
      call MPI_WAIT(recv_req, status, ier)
      call MPI_BARRIER(comm,ier)
      call imagdp(temp)

      return
      end
#endif
c
       subroutine collect_xy(temp,phi,klev)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a  2d field at level klev. Image points are
c collected
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req
      integer i,j,klev
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1,0:kdim1)
      real phi(idim,jtotal),tphi(0:idim1,jtotal)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c     we don't need a special data structure since the data blocks are
c     continuous for a zlevel
c
c-----------------------------------------------------------------------
c
c
      if (iope) then
        do np=1,nprocs
c
c post receives first, use temporary array to avoid image point
c problems
c
          call MPI_IRECV(tphi(0,blockstart(np)), idim2*jdim, mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
        enddo
      endif
c
c now send data chunks from each node
c
      call MPI_ISEND(temp(0,1,klev),idim2*jdim,mpi_real, mastertask,
     &             mpitag_io, comm, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c
        do j=1,jtotal
          do i=1,idim
            phi(i,j) = tphi(i,j)
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)

      return
      end
c
c
      subroutine collect_yz(temp,phi,islb)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a 2d y-z slab for a constant islb value, 
c including the image points.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req,k
      integer j,islb
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real temp(0:idim1,0:jdim1,0:kdim1),ttemp(kdim,jdim)
      real phi(jtotal,kdim),tphi(kdim,jtotal)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c
c
      if (iope) then
        do np=1,nprocs
c
c post receives first, overwrite the image points to stay with
c a contiguous block of data
c
          call MPI_IRECV(tphi(1,blockstart(np)), kdim*jdim,mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
        enddo
      endif
      do k=1,kdim
        do j=1,jdim
          ttemp(k,j) = temp(islb,j,k)
        enddo
      enddo
      call MPI_ISEND(ttemp(1,1), jdim*kdim,mpi_real, 
     &             mastertask,
     &             mpitag_io, comm, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
        do j=1,jtotal
          do k=1,kdim
            phi(j,k) = tphi(k,j)
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)

      return
      end


       subroutine collect_xz(temp,phi,jslb)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from the processor containing j to the 
c master task. Easiest slab routine.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req,send_req,k
      integer i,jslb,jtop,jbot,jlocal
      integer status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1,0:kdim1)
      real phi(idim,kdim),tphi(0:idim1,0:kdim1)
c
c
      jbot = blockstart(mastertask+1)
      jtop = blockstart(mastertask+1)+jdim
      if(jslb.ge.jbot.and.jslb.le.jtop) then
c
c slab is on the mastertask, so do straight copy
c
        if (iope) then
          write(*,*)"Collect xz pe number ", my_pe,iope,jslb
          do i=1,idim 
            do k=1,kdim 
              phi(i,k) = temp(i,jslb,k) 
            enddo 
          enddo
        endif
      else
c
c need to transfer slab from other processor
c
        if(iope) then
c
c first figure out processor and post receive
c
          do np=1,nprocs
            jbot = blockstart(np)
            jtop = blockstart(np)+jdim
            if(jslb.ge.jbot.and.jslb.le.jtop) then
              call MPI_IRECV(tphi(0,0), idim2*kdim2,
     &                mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req, ier)
            endif
          enddo
        endif
c
c then see if current processor is the sender and post send if it is
c
        jbot = blockstart(my_pe+1)
        jtop = blockstart(my_pe+1)+jdim
        if(jslb.ge.jbot.and.jslb.le.jtop) then
          jlocal = jslb-jbot+1
          call MPI_ISEND(temp(0,jlocal,0), 1,ns_type_real, mastertask,
     &             mpitag_io, comm, send_req, ier)
          call MPI_WAIT(send_req, status, ier)
        endif
c
c make sure everything is communicated before moving on
c
        call MPI_BARRIER(comm,ier)
        if(iope) then
           do k=1,kdim
              do i=1,idim
                 phi(i,k) = tphi(i,k)
              enddo
           enddo
        endif
      endif

      return
      end
c
c
       subroutine collect_point(ix,jx,kx,buf)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from the processor containing j to the
c master task. Easiest slab routine.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req,send_req
      integer jtop,jbot,ix,jx,kx,jlocal
      integer status(MPI_STATUS_SIZE)

      real buf(11)
c
c first see if point in on mastertask processor
c
      jbot = blockstart(mastertask+1)
      jtop = blockstart(mastertask+1)+jdim
      if(jx.ge.jbot.and.jx.le.jtop) then
c
c slab is on the mastertask, so do straight copy
c
        if (iope) then
          jlocal = jx-jbot+1
          buf(1) = um(ix,jlocal,kx)
          buf(2) = vm(ix,jlocal,kx)
          buf(3) = wm(ix,jlocal,kx)
          buf(4) = p(ix,jlocal,kx)
          buf(5) = th(ix,jlocal,kx)
          buf(6) = sal(ix,jlocal,kx)
          buf(7) = tke(ix,jlocal,kx)**2
          buf(8) = kheat(ix,jlocal,kx)
          buf(9) = kmom(ix,jlocal,kx)
          buf(10) = rhop(ix,jlocal,kx)
#if defined(TRACER)
          buf(11) = trc(ix,jlocal,kx,1)
#endif
        endif
      else
c
c need to transfer point from other processor
c
        if(iope) then
c
c first figure out processor and post receive
c we use a buffer to transfer the point values
c
          do np=1,nprocs
            jbot = blockstart(np)
            jtop = blockstart(np)+jdim
            if(jx.ge.jbot.and.jx.le.jtop) then
              call MPI_IRECV(buf, 11,
     &                mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req, ier)
            endif
          enddo
        endif
c
c then see if current processor is the sender and post send if it is
c
        jbot = blockstart(my_pe+1)
        jtop = blockstart(my_pe+1)+jdim
        if(jx.ge.jbot.and.jx.le.jtop) then
          jlocal = jx-jbot+1
          buf(1) = um(ix,jlocal,kx)
          buf(2) = vm(ix,jlocal,kx)
          buf(3) = wm(ix,jlocal,kx)
          buf(4) = p(ix,jlocal,kx)
          buf(5) = th(ix,jlocal,kx)
          buf(6) = sal(ix,jlocal,kx)
          buf(7) = tke(ix,jlocal,kx)**2
          buf(8) = kheat(ix,jlocal,kx)
          buf(9) = kmom(ix,jlocal,kx)
          buf(10) = rhop(ix,jlocal,kx)
#if defined(TRACER)
          buf(11) = trc(ix,jlocal,kx,1)
#endif
          call MPI_ISEND(buf, 11, mpi_real, mastertask,
     &             mpitag_io, comm, send_req, ier)
          call MPI_WAIT(send_req, status, ier)
        endif
      endif
c
c make sure everything is communicated before moving on
c
      call MPI_BARRIER(comm,ier)
c
c copy values from buffer
c

      return
      end

c
       subroutine collect_2d(temp,phi)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a full 2d field, including the image points.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req
      integer i,j
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1)
      real phi(0:idim1,0:jtotal1)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c
c
      if (iope) then
        do np=1,nprocs
c
c post receives first
c
          call MPI_IRECV(phi(0,blockstart(np)), idim2*jdim,mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
        enddo
      endif
      call MPI_ISEND(temp(0,1), idim2*jdim, mpi_real, mastertask,
     &             mpitag_io, comm, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      call MPI_BARRIER(comm,ier)
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c need to do j image points
c
        do i=0,idim1
          phi(i,0) = phi(i,jtotal)
          phi(i,jtotal1) = phi(i,1)
        enddo
        do j=0,jtotal1
	  phi(idim1,j) = phi(1,j)
	  phi(0,j) = phi(idim,j)
	enddo
      endif

      return
      end

      subroutine distrib_2d(temp,phi)

      implicit none

      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
 
c-----------------------------------------------------------------------
c
c This routine takes a full 2-d field (phi) and ships it out to each of
c the processors through "temp"
c
c-----------------------------------------------------------------------
 
      integer np,ier,recv_req, send_req
      integer status(MPI_STATUS_SIZE)
 
      real   temp(0:idim1,0:jdim1)
      real phi(0:idim1,0:jtotal1)
c-----------------------------------------------------------------------
c
c post recieves for each of the data chunks
c
c-----------------------------------------------------------------------
 
      call MPI_IRECV(temp(0,1), idim2*jdim, mpi_real, mastertask,
     &              mpitag_io, comm, recv_req, ier)
c
c now send out parts of global variable to corresponding processors
c based on global address space defined in setupmpi routine
c
      if (iope) then
        do np=1,nprocs
c
c post the send for this processors data
c
          call MPI_ISEND(phi(0,blockstart(np)), idim2*jdim,
     &                  mpi_real, np-1,
     &                  mpitag_io, comm, send_req, ier)
          call MPI_WAIT(send_req, status, ier)
        enddo
      endif
c
c wait for everything to communicate
c
      call MPI_WAIT(recv_req, status, ier)
      call MPI_BARRIER(comm,ier)
 
      return
      end
c
c 
       subroutine collect_2ds(temp,phi)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a  2d field. Image points are
c collected
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req
      integer i,j
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1)
      real phi(idim,jtotal),tphi(0:idim1,jtotal)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c     we don't need a special data structure since the data blocks are
c     continuous for a zlevel
c
c-----------------------------------------------------------------------
c
c
      if (iope) then
        do np=1,nprocs
c
c post receives first, use temporary array to avoid image point
c problems
c
          call MPI_IRECV(tphi(0,blockstart(np)), idim2*jdim, mpi_real,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
        enddo
      endif
c
c now send data chunks from each node
c
      call MPI_ISEND(temp(0,1),idim2*jdim,mpi_real, mastertask,
     &             mpitag_io, comm, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c
        do j=1,jtotal
          do i=1,idim
            phi(i,j) = tphi(i,j)
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)

      return
      end
