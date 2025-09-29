      subroutine collect_yz_ave(temp,phi,jslb)

      implicit none
      include 'mpif.h'
#include "param.inc"
#include "communicate.inc"
c
c-----------------------------------------------------------------------
c
c This routine sends data from each processor to the master task
c and combines it into a 2d y-z slab for a constant islb value, 
c excluding the image points.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req,k
      integer i,j,jslb,jbot,jtop,ibot,itop
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real temp(0:idim1,0:jdim1,0:kdim1)
      real phi(jtotal,kdim1),tphi(jdim,kdim1)
      real phip(jtotal,kdim1)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c           
c      write(*,*)"my_pe jslb jglb ",my_pe,jslb,(jglobal(j),j=1,jdim)
      if (iope) then
        do k=1,kdim1
          do j=1,jtotal
            phip(j,k) = 0.0
          enddo
        enddo
        do np=1,nprocs
c
c post receives for processors that have a section of islb
c          write(*,*)"iblkstart, iope " ,iblkstart(np),np
c
          if(np-1.ne.my_pe) then
c        write(*,*)"rec jslb, np ",jslb, np-1
c
              call MPI_RECV(tphi(1,1), kdim*jdim,MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, status, ier)
c
c copy aveerage section to full slab
c
              do k=1,kdim
                do j=1,jdim
                  phip(jblkstart(np)+j-1,k) = 
     *                   phip(jblkstart(np)+j-1,k)+tphi(j,k)
                enddo
              enddo
          endif
        enddo
      else
c compute average along j axis and send to mastertask
        do k=1,kdim
          do j=1,jdim
            tphi(j,k) = 0.0
          enddo
        enddo
        do k=1,kdim
          do i=1,idim
            do j=1,jdim
              tphi(j,k) = tphi(j,k)+temp(i,j,k)
            enddo
          enddo      
        enddo
        call MPI_SEND(tphi(1,1), jdim*kdim,MPI_REAL, 
     &             mastertask,
     &             mpitag_io, comm, ier)
c        call MPI_WAIT(send_req, status, ier)
      endif
c
c make sure everything is communicated before moving on
c
      if(iope) then
c        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c if a section of the slab is on iope, copy it into the output array
c
          do k=1,kdim
            do i=1,idim
              do j=1,jdim
                phip(jblkstart(my_pe+1)+j-1,k) =
     *                phip(jblkstart(my_pe+1)+j-1,k)+temp(i,j,k)
              enddo
            enddo
          enddo
      endif
c
c compute average along the y axis
c
      do k=1,kdim
        do j=1,jtotal
          phip(j,k) = phip(j,k)/itotal
        enddo
      enddo
      call MPI_BARRIER(comm,ier)
c
c next we need to send out the x-z phi array to all of the processors
c
c first send out parts of global variable to corresponding processors
c based on global address space defined in setupmpi routine
c
      if (iope) then
        do np=1,nprocs
c
c get limits for processor np-1
c
c           ibot = iblkstart(np)
c           itop = iblkstart(np)+idim-1

c           write(*,*)"ibot,itop,jbot,jtop ",ibot,itop,jbot,jtop 
c           do k=1,kdim
c             do i=ibot,itop
c               tphi(i-ibot+1,k) = phi(i,k)
c             enddo
c           enddo
           if(np-1.ne.my_pe) then
c
c post the send for this processors data
c

             call MPI_SEND(phip, jtotal*kdim1,
     &                  MPI_REAL, np-1,
     &                  mpitag_io, comm, ier)
c             call MPI_WAIT(send_req, status, ier)
            endif

        enddo
      else
c
c if not iope, then collect chunk sent from iope
c
        call MPI_RECV(phi(1,1), jtotal*kdim1, MPI_REAL, mastertask,
     &              mpitag_io, comm, status, ier)
c        call MPI_WAIT(recv_req, status, ier)
         
      endif
c
c make copy of array for mastertask
c
      if(iope) then
        do k=1,kdim1
          do j=1,jtotal
            phi(j,k) = phip(j,k)
          enddo
        enddo
      endif
c
c wait for everything to communicate
c
      call MPI_BARRIER(comm,ier)




      return
      end
c
c
