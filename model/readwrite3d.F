      subroutine write3d(cdfout,cdfid,rstart,rcount,fld,phi)
c
c Routine to write out netcdf file in chunks to avoid large array
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
c

      real fld(0:idim1,0:jdim1,0:kdim1)
      real phi(0:idim1,0:jdim1,0:kdim1)
      integer cdfout,cdfid,start(4),count(4)
      integer rstart(4),rcount(4),rcode
      integer np,ier,recv_req(nprocs),send_req,k,nprec,nps
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)


c
      count(1) = idim+2
      count(2) = jdim+2
      count(3) = kdim+2
      count(4) = 1
      if(iope) then
c        write(*,*)"Processing data chunks on iope"
        do np=1,nprocs
          start(1) = 1
          if(np.eq.1) then
            start(2) = 1
          else
            start(2) = jdim*(np-1)+1
          endif
          start(3) = 1
          start(4) = 1
          do nps=1,nprocs
            if(blockcord(nps).eq.np-1) then
               nprec = nps-1
            endif
          enddo
c          write(*,*)"Post Recv ",nprec
          if(nprec.ne.0)then
            call MPI_RECV(phi(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $         nprec,mpitag_io,comm,status,ier)
            call ncvpt(cdfout, cdfid, start, count, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite", nprec
          else
c            write(*,*)"writing iope chunk"
            call ncvpt(cdfout, cdfid, start, count, fld, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite on iope"
          endif            
c          write(*,*)"Output phi ",phi(1,1,1),nprec,start(2)
        enddo
      else
c        write(*,*)"Post Sends ",my_pe
        call MPI_SEND(fld(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $               mastertask, mpitag_io,comm,ier)
      endif
      call MPI_BARRIER(comm,ier)
      return
      end
c
c
      subroutine read3d(cdfin,cdfid,rstart,rcount,fld,phi)
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      real fld(0:idim1,0:jdim1,0:kdim1)
      real phi(0:idim1,0:jdim1,0:kdim1)
      integer cdfin,cdfid,start(4),count(4),rcode
      integer rstart(4),rcount(4)
      integer np,ier,recv_req, send_req,i,k,nprec,nps
      integer status(MPI_STATUS_SIZE)


c
      count(1) = idim+2
      count(2) = jdim+2
      count(3) = kdim+2
      count(4) = 1
      if(iope) then
        do np=1,nprocs
          start(1) = 1
          if(np.eq.1) then
            start(2) = 1
          else
            start(2) = jdim*(np-1)+1
          endif
          start(3) = 1
          start(4) = 1
          do nps=1,nprocs
            if(blockcord(nps).eq.np-1) then
               nprec = nps-1
            endif
          enddo
          if(nprec.ne.0) then
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            call MPI_SEND(phi(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $         nprec,
     $         mpitag_io,comm,ier)
          else
            call ncvgt(cdfin,cdfid,start,count,fld,rcode)
          endif
        enddo
      else
        call MPI_RECV(fld(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $                mastertask,mpitag_io,comm,status,ier)
      endif
      call MPI_BARRIER(comm,ier)
      write(*,*)"my pe ",my_pe,phi(1,1,1)
      return
      end



