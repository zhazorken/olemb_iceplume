      subroutine collect_xz_ave(temp,phi,jslb)

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
      integer i,j,jslb,jbot,jtop
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real temp(0:idim1,0:jdim1,0:kdim1)
      real phi(itotal,kdim1),tphi(idim,kdim)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c           
c      write(*,*)"my_pe jslb jglb ",my_pe,jslb,(jglobal(j),j=1,jdim)
      if (iope) then
        do k=1,kdim1
          do i=1,itotal
            phi(i,k) = 0.0
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
              call MPI_RECV(tphi(1,1), kdim*idim,MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
c
c copy aveerage section to full slab
c
              do k=1,kdim
                do i=1,idim
                  phi(iblkstart(np)+i-1,k) = 
     *                   phi(iblkstart(np)+i-1,k)+tphi(i,k)
                enddo
              enddo
          endif
        enddo
      else
c compute average along j axis and send to mastertask
        do k=1,kdim
          do i=1,idim
            tphi(i,k) = 0.0
          enddo
        enddo
        do k=1,kdim
          do i=1,idim
            do j=1,jdim
              tphi(i,k) = tphi(i,k)+temp(i,j,k)
            enddo
          enddo      
        enddo
        call MPI_SEND(tphi(1,1), idim*kdim,MPI_REAL, 
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
                phi(iblkstart(my_pe+1)+i-1,k) =
     *                phi(iblkstart(my_pe+1)+i-1,k)+temp(i,j,k)
              enddo
            enddo
          enddo
      endif
c
c compute average along the y axis
c
      do k=1,kdim
        do i=1,itotal
          phi(i,k) = phi(i,k)/jtotal
        enddo
      enddo
      call MPI_BARRIER(comm,ier)

      return
      end
