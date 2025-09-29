      subroutine read3d(cdfin,cdfid,rstart,rcount,fld,phi)
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      integer idimh,jdimh,kdimh,idimh1,jdimh1,kdimh1
      integer idimh2,jdimh2,kdimh2
      parameter(idimh=idim/2,jdimh=jdim/2,kdimh=kdim/2,
     *          idimh1=idimh+1,jdimh1=jdimh+1,kdimh1=kdimh+1,
     *          idimh2=idimh+2,jdimh2=jdimh+2,kdimh2=kdimh+2)
      real fld(0:idimh1,0:jdimh1,0:kdimh1)
      real phi(0:idimh1,0:jdimh1,0:kdimh1)
      integer cdfin,cdfid,start(4),count(4),rcode
      integer rstart(4),rcount(4)
      integer np,ier,recv_req, send_req,i,k,nprec,nps
      integer status(MPI_STATUS_SIZE)


c
      count(1) = idimh+2
      count(2) = jdimh+2
      count(3) = kdimh+2
      count(4) = 1
      if(iope) then
        do np=1,nprocs
          start(1) = 1
          if(np.eq.1) then
            start(2) = 1
          else
            start(2) = (jdimh)*(np-1)+1
          endif
          start(3) = 1
          start(4) = 1
          do nps=1,nprocs
            if(blockcord(nps).eq.np-1) then
               nprec = nps-1
            endif
          enddo
          if(nprec.ne.0) then
c            write(*,*)"getting phi ",nprec,count(1),count(2),count(3),
c     *          start(2)
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            call MPI_SEND(phi(0,0,0),idimh2*jdimh2*kdimh2,MPI_REAL,
     $         nprec,
     $         mpitag_io,comm,ier)
          else
c            write(*,*)"getting fld ",nprec,count(1),count(2),count(3),
c     *          start(2)
            call ncvgt(cdfin,cdfid,start,count,fld,rcode)
          endif
        enddo
      else
        call MPI_RECV(fld(0,0,0),idimh2*jdimh2*kdimh2,MPI_REAL,
     $                mastertask,mpitag_io,comm,status,ier)
      endif
      call MPI_BARRIER(comm,ier)
c      write(*,*)"my pe ",my_pe,phi(1,1,1)
      return
      end



