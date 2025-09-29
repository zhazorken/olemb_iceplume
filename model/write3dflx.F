      subroutine write3dflx(cdfout,cdfid,fld)
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
      real phi(idim,jdim,kdim)
      integer cdfout,cdfid,startz(4),countz(4)
      integer rcode
#if defined(USEMPI)
      integer np,ier,recv_req(nprocs),send_req,k,nprec,nps
      integer npx,npy
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)
#endif


c
      countz(1) = idim
      countz(2) = jdim
      countz(3) = kdim
      countz(4) = 1
      startz(1) = 1
      startz(2) = 1
      startz(3) = 1
      startz(4) = 1
c
      startz(1) = iblkstart(my_pe+1)
      startz(2) = jblkstart(my_pe+1)
c
c limit x axis to itotal/4.
c
c      if(startz(1).lt.itotal/3.) then
      write(*,*)"flx write ",my_pe
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            phi(i,j,k) = fld(i,j,k)
          enddo
        enddo
      enddo
c      if(my_pe.eq.1) then
c        write(*,*)"write3d,varid,phi(3,3) ",cdfid,phi(3,3),startz(1),
c     *          startz(2)
c      endif
c      call ncvpt(cdfout,cdfid,start,count,phi,rcode)
c      write(*,*)"my_pe,startz,phi(2,2)",my_pe,cdfid,startz(1),startz(2),
c     *           phi(2,2,4)
c      rcode = nf_var_par_access(cdfout,cdfid,nf_collective)
c      if(rcode.ne.0) write(*,*)"Error in par_access in write3d"
      rcode = nf_put_vara_real(cdfout,cdfid,startz,countz,phi)
      if(rcode.ne.0) write(*,*)"Error in write3dflux",rcode,my_pe
c      call MPI_BARRIER(comm,ier)
c      endif
c
      return
      end
c
