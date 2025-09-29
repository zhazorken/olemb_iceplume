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

      integer np,npx,npy,ier,recv_req(nprocs),send_req
      integer i,j,k,klev
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1,0:kdim1)
      real phi(itotal,jtotal),tphi(idim,jdim)
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
          if(np-1.ne.mastertask) then
            call MPI_RECV(tphi(1,1), 
     &                idim*jdim, MPI_REAL, np-1,
     &                mpitag_io, comm, recv_req(np), ier)
            do i=0,idim-1
              do j=0,jdim-1
                phi(iblkstart(np)+i,jblkstart(np)+j) = tphi(i+1,j+1)
              enddo
            enddo
          endif
        enddo
      else
c
c now send data chunks from each node
c
        do j=1,jdim
          do i=1,idim
c            tphi(i,j) = 0.0
c            do k=1,12
c            tphi(i,j) = temp(i,j,k)+tphi(i,j)
c            enddo
c            tphi(i,j) = tphi(i,j)/12.0
            tphi(i,j) = temp(i,j,klev)
          enddo
        enddo
        call MPI_SEND(tphi(1,1),idim*jdim,MPI_REAL, mastertask,
     &             mpitag_io, comm, ier)
c        call MPI_WAIT(send_req, status, ier)
      endif
c
c make sure everything is communicated before moving on
c
      if(iope) then
c        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c finally, copy local data to global array
c
         do i=0,idim-1
           do j=0,jdim-1
              phi(iblkstart(my_pe+1)+i,jblkstart(my_pe+1)+j) = 
     *          temp(i+1,j+1,klev)
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
c excluding the image points.
c
c-----------------------------------------------------------------------

      integer np,ier,recv_req(nprocs),send_req,k
      integer j,islb,ibot,itop
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real temp(0:idim1,0:jdim1,0:kdim1)
      real phi(jtotal,kdim1),tphi(jdim,kdim)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c-----------------------------------------------------------------------
c           
      ibot = iblkstart(my_pe+1)
      itop = iblkstart(my_pe+1)+idim-1
      if (iope) then
        do np=1,nprocs
c
c post receives for processors that have a section of islb
c
          if(np-1.ne.my_pe) then
            if(islb.ge.iblkstart(np).and.islb.lt.iblkstart(np)+idim)
     &       then
              call MPI_RECV(tphi(1,1), kdim*jdim,MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
c
c copy subsection to full slab
c
              do k=1,kdim
                do j=1,jdim
                  phi(jblkstart(np)+j-1,k) = tphi(j,k)
                enddo
              enddo
            endif
          endif
        enddo

      else if(islb.ge.ibot.and.islb.le.itop)then

        do k=1,kdim
          do j=1,jdim
            tphi(j,k) = temp(islb-ibot+1,j,k)
          enddo      
        enddo
        call MPI_SEND(tphi(1,1),jdim*kdim,MPI_REAL, 
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
        if(islb.ge.ibot.and.islb.le.itop) then
          do k=1,kdim
            do j=1,jdim
              phi(jblkstart(my_pe+1)+j-1,k) = temp(islb-ibot+1,j,k)
            enddo
          enddo
        endif
      endif
      call MPI_BARRIER(comm,ier)

      return
      end
c
c now routine for xz slab
c

      subroutine collect_xz(temp,phi,jslb)

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
        do np=1,nprocs
c
c post receives for processors that have a section of islb
c          write(*,*)"iblkstart, iope " ,iblkstart(np),np
c
          if(np-1.ne.my_pe) then
            if(jslb.ge.jblkstart(np).and.jslb.lt.jblkstart(np)+jdim)
     &        then
c        write(*,*)"rec jslb, np ",jslb, np-1
c
              call MPI_RECV(tphi(1,1), kdim*idim,MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req(np), ier)
c
c copy subsection to full slab
c
              do k=1,kdim
                do i=1,idim
                  phi(iblkstart(np)+i-1,k) = tphi(i,k)
                enddo
              enddo
            endif
          endif
        enddo
      else if(jslb.ge.jblkstart(my_pe+1).and.
     $        jslb.lt.jblkstart(my_pe+1)+jdim)then

c        write(*,*)"send jslb, my_pe ",jslb,my_pe
        do k=1,kdim
          do i=1,idim
            tphi(i,k) = temp(i,jslb-jblkstart(my_pe+1)+1,k)
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
        jbot = jblkstart(my_pe+1)
        jtop = jblkstart(my_pe+1)-1+jdim
        if(jslb.ge.jbot.and.jslb.le.jtop) then
          do k=1,kdim
            do i=1,idim
              phi(iblkstart(my_pe+1)+i-1,k) = temp(i,jslb-jbot+1,k)
            enddo
          enddo
        endif
      endif
      call MPI_BARRIER(comm,ier)

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

      integer np,ier,send_req
      integer itop,ibot,jtop,jbot,ix,jx,kx,jlocal,ilocal
      integer status(MPI_STATUS_SIZE)
      integer recv_req(MPI_STATUS_SIZE)

      real buf(11)
c
c first see if point in on mastertask processor
c
      jbot = jblkstart(mastertask+1)
      ibot = iblkstart(mastertask+1)
      jtop = jblkstart(mastertask+1)+jdim
      itop = iblkstart(mastertask+1)+idim
      if(jx.ge.jbot.and.jx.lt.jtop.and.ix.ge.ibot.and.ix.lt.itop) then
c
c slab is on the mastertask, so do straight copy
c
        if (iope) then
          jlocal = jx-jbot+1
          ilocal = ix-ibot+1
          buf(1) = um(ilocal,jlocal,kx)
          buf(2) = vm(ilocal,jlocal,kx)
          buf(3) = wm(ilocal,jlocal,kx)
          buf(4) = p(ilocal,jlocal,kx)
          buf(5) = th(ilocal,jlocal,kx)
          buf(6) = sal(ilocal,jlocal,kx)
          buf(7) = tke(ilocal,jlocal,kx)**2
          buf(8) = kheat(ilocal,jlocal,kx)
          buf(9) = kmom(ilocal,jlocal,kx)
          buf(10) = rhop(ilocal,jlocal,kx)
#if defined(TRACER)
          buf(11) = trc(ilocal,jlocal,kx,1)
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
            jbot = jblkstart(np)
            ibot = iblkstart(np)
            jtop = jblkstart(np)+jdim
            itop = iblkstart(np)+idim
            if(jx.ge.jbot.and.jx.lt.jtop.and.
     &         ix.ge.ibot.and.ix.lt.itop) then
              call MPI_RECV(buf, 11,
     &                MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req, ier)
            endif
          enddo
        endif
c
c then see if current processor is the sender and post send if it is
c
        jbot = jblkstart(my_pe+1)
        ibot = iblkstart(my_pe+1)
        jtop = jblkstart(my_pe+1)+jdim
        itop = iblkstart(my_pe+1)+idim

        if(jx.ge.jbot.and.jx.lt.jtop.and.
     &         ix.ge.ibot.and.ix.lt.itop) then

          jlocal = jx-jbot+1
          ilocal = ix-ibot+1
          buf(1) = um(ilocal,jlocal,kx)
          buf(2) = vm(ilocal,jlocal,kx)
          buf(3) = wm(ilocal,jlocal,kx)
          buf(4) = p(ilocal,jlocal,kx)
          buf(5) = th(ilocal,jlocal,kx)
          buf(6) = sal(ilocal,jlocal,kx)
          buf(7) = tke(ilocal,jlocal,kx)**2
          buf(8) = kheat(ilocal,jlocal,kx)
          buf(9) = kmom(ilocal,jlocal,kx)
          buf(10) = rhop(ilocal,jlocal,kx)
#if defined(TRACER)
          buf(11) = trc(ilocal,jlocal,kx,1)
#endif
          call MPI_SEND(buf, 11, MPI_REAL, mastertask,
     &             mpitag_io, comm, ier)
c          call MPI_WAIT(send_req, status, ier)
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
c and combines it into a full 2d field. Image points are not
c collected
c
c-----------------------------------------------------------------------

      integer np,npx,npy,ier,recv_req(nprocs),send_req
      integer i,j
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)

      real   temp(0:idim1,0:jdim1)
      real phi(itotal,jtotal),tphi(idim,jdim)
c-----------------------------------------------------------------------
c
c     receive chunks from processors who own them
c
c     we don't need a special data structure since the data blocks are
c     continuous for a zlevel
c
c-----------------------------------------------------------------------
c
c first copy to buffer area
c
      do j=1,jdim
        do i=1,idim
          tphi(i,j) = temp(i,j)
        enddo
      enddo
c
      if (iope) then
        do np=1,nprocs
c
c post receives first, use temporary array to avoid image point
c problems
c
          if(my_pe.ne.np-1) then
            call MPI_RECV(tphi(1,1), 
     &                idim*jdim, MPI_REAL, np-1,
     &                mpitag_io, comm, recv_req(np), ier)
            do i=0,idim-1
              do j=0,jdim-1
                phi(iblkstart(np)+i,jblkstart(np)+j) = tphi(i+1,j+1)
              enddo
            enddo
          endif
        enddo
      else
c
c now send data chunks from each node
c
        call MPI_SEND(tphi(1,1),idim*jdim,MPI_REAL, mastertask,
     &             mpitag_io, comm, ier)
c        call MPI_WAIT(send_req, status, ier)
      endif
c
c make sure everything is communicated before moving on
c
      if(iope) then
c        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
c
c finally, copy local data to global array
c
        do j=0,jdim-1
          do i=0,idim-1
            phi(iblkstart(my_pe+1)+i,jblkstart(my_pe+1)+j) = 
     *         temp(i+1,j+1)
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)

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
      integer i,j,ibot,jbot,itop,jtop
      integer status(MPI_STATUS_SIZE)
 
      real temp(1:idim,1:jdim)
      real phi(1:itotal,1:jtotal)
c
c-----------------------------------------------------------------------
c first send out parts of global variable to corresponding processors
c based on global address space defined in setupmpi routine
c
      if (iope) then
        do np=1,nprocs
c
c get limits for processor np-1
c
           jbot = jblkstart(np)
           ibot = iblkstart(np)
           jtop = jblkstart(np)+jdim-1
           itop = iblkstart(np)+idim-1

           do j=jbot,jtop
             do i=ibot,itop
               temp(i-ibot+1,j-jbot+1) = phi(i,j)
             enddo
           enddo
           if(np-1.ne.my_pe) then
c
c post the send for this processors data
c
          
             call MPI_SEND(temp, idim*jdim,
     &                  MPI_REAL, np-1,
     &                  mpitag_io, comm, ier)
c             call MPI_WAIT(send_req, status, ier)
           endif
            
        enddo
      else
c
c if not iope, then collect chunk sent from iope
c
        call MPI_RECV(temp(1,1), idim*jdim, MPI_REAL, mastertask,
     &              mpitag_io, comm, recv_req, ier)
c        call MPI_WAIT(recv_req, status, ier)
      endif
c
c wait for everything to communicate
c
      call MPI_BARRIER(comm,ier)
 
      return
      end
c
c read in 3-d files
c
      subroutine read3d(cdfin,cdfid,fld)
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
      real phi(idimh,jdimh,kdimh)
      integer cdfin,cdfid,start(4),count(4),rcode
      integer rstart(4),rcount(4)
      integer np,ier,recv_req, send_req,i,k,nprec,nps
      integer status(MPI_STATUS_SIZE)


c
      count(1) = idimh
      count(2) = jdimh
      count(3) = kdimh
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1

c
c #if defined(USEMPI)
c #if defined(USEMPI_OLDNETCDF)
      if(iope) then
        write(*,*)"input count ",idimh,jdimh,kdimh
        do npx=1,nprocx
        do npy=1,nprocy
          start(1) = idimh*(npx-1)+1
          start(2) = jdimh*(npy-1)+1
          start(3) = 1
          start(4) = 1
          do np=1,nprocs
            if(jblkcord(np).eq.npy-1.and.iblkcord(np).eq.npx-1)
     $       then
               nprec = np-1
            endif
          enddo
          if(nprec.ne.my_pe) then
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            call MPI_SEND(phi(1,1,1),idimh*jdimh*kdimh,MPI_REAL,
     $         nprec,
     $         mpitag_io,comm,ier)
          else
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            do k=1,kdimh
              do j=1,jdimh
                do i=1,idimh
                  fld(i,j,k) = phi(i,j,k)
                enddo
              enddo
            enddo
          endif
        enddo
        enddo
      else
        call MPI_RECV(phi(1,1,1),idimh*jdimh*kdimh,MPI_REAL,
     $                mastertask,mpitag_io,comm,status,ier)
        do k=1,kdimh
          do j=1,jdimh
            do i=1,idimh
              fld(i,j,k) = phi(i,j,k)
            enddo
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)
c      write(*,*)"my pe ",my_pe,phi(1,1,1)
c #if defined(NEWUSEMPI)
#if defined(NEWUSEMPI)
      start(1) = iblkstart(my_pe+1)
      start(2) = jblkstart(my_pe+1)
c      write(*,*)"reading input file on ",my_pe,start(1),start(2)
      call ncvgt(cdfin,cdfid,start,count,phi,rcode)
c      write(*,*)"finished reading on ",my_pe
      do k=1,kdimh
        do j=1,jdimh
          do i=1,idimh
            fld(i,j,k) = phi(i,j,k)
          enddo
        enddo
      enddo
      call MPI_BARRIER(comm,ier)
#endif
      return
      end



      subroutine write3d(cdfout,cdfid,fld)
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
#if defined(USEMPI_OLDNETCDF)
      if(iope) then
c        write(*,*)"Processing data chunks on iope"
        do npx=1,nprocx
        do npy=1,nprocy
c
          startz(1) = idim*(npx-1)+1
          startz(2) = jdim*(npy-1)+1
          startz(3) = 1
          startz(4) = 1
          do nps=1,nprocs
            if(jblkcord(nps).eq.npy-1.and.iblkcord(nps).eq.npx-1)
     *       then
               nprec = nps-1
            endif
          enddo
c          write(*,*)"Post Recv ",nprec
          if(nprec.ne.mastertask)then
            call MPI_RECV(phi,idim*jdim*kdim,MPI_REAL,
     $         nprec,mpitag_io,comm,status,ier)
            call ncvpt(cdfout, cdfid, startz, countz, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite", nprec
          else
c            write(*,*)"writing iope chunk"
            do k=1,kdim
              do j=1,jdim
                do i=1,idim
                  phi(i,j,k) = fld(i,j,k)
                enddo
              enddo
            enddo
            call ncvpt(cdfout, cdfid, startz, countz, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite on iope"
          endif            
c          write(*,*)"Output phi ",phi(1,1,1),nprec,startz(2)
        enddo
        enddo
      else
c        write(*,*)"write3d Post Sends ",my_pe
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              phi(i,j,k) = fld(i,j,k)
            enddo
          enddo
        enddo
        call MPI_SEND(phi,idim*jdim*kdim,MPI_REAL,
     $               mastertask, mpitag_io,comm,ier)
      endif
      call MPI_BARRIER(comm,ier)
c #if defined(NEWUSEMPI)
#else
      startz(1) = iblkstart(my_pe+1)
      startz(2) = jblkstart(my_pe+1)
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
      if(rcode.ne.0) write(*,*)"Error in ncwrite in write3d",rcode,my_pe
c      call MPI_BARRIER(comm,ier)
#endif
      return
      end
c
c
      subroutine write2d(cdfout,cdfid,fld)
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

      real fld(0:idim1,0:jdim1)
      real phi(idim,jdim)
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
      countz(3) = 1
      countz(4) = 1
      startz(1) = 1
      startz(2) = 1
      startz(3) = 1
      startz(4) = 1
c
#if defined(USEMPI_OLDNETCDF)
      if(iope) then
c        write(*,*)"Processing data chunks on iope"
        do npx=1,nprocx
        do npy=1,nprocy
c
          startz(1) = idim*(npx-1)+1
          startz(2) = jdim*(npy-1)+1
          startz(3) = 1
          startz(4) = 1
          do nps=1,nprocs
            if(jblkcord(nps).eq.npy-1.and.iblkcord(nps).eq.npx-1)
     *       then
               nprec = nps-1
            endif
          enddo
c          write(*,*)"Post Recv ",nprec
          if(nprec.ne.mastertask)then
            call MPI_RECV(phi,idim*jdim,MPI_REAL,
     $         nprec,mpitag_io,comm,status,ier)
            call ncvpt(cdfout, cdfid, startz, countz, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite", nprec
          else
            do j=1,jdim
              do i=1,idim
                phi(i,j) = fld(i,j)
              enddo
            enddo
            call ncvpt(cdfout, cdfid, startz, countz, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite on iope"
          endif            
c          write(*,*)"Output phi ",phi(1,1),nprec,startz(2)
        enddo
        enddo
      else
c        write(*,*)"write2d Post Sends ",my_pe
        do j=1,jdim
          do i=1,idim
            phi(i,j) = fld(i,j)
          enddo
        enddo
        call MPI_SEND(phi,idim*jdim,MPI_REAL,
     $               mastertask, mpitag_io,comm,ier)
      endif
      call MPI_BARRIER(comm,ier)
c#if defined(NEWUSEMPI)
#else
      startz(1) = iblkstart(my_pe+1)
      startz(2) = jblkstart(my_pe+1)
        do j=1,jdim
          do i=1,idim
            phi(i,j) = fld(i,j)
          enddo
        enddo
c      if(my_pe.eq.1) then
c        write(*,*)"write2d,varid,phi(3,3) ",cdfid,phi(3,3),startz(1),
c     *          startz(2)
c      endif
c      call ncvpt(cdfout,cdfid,start,count,phi,rcode)
c      write(*,*)"my_pe,startz,phi(2,2)",my_pe,cdfid,startz(1),startz(2),
c     *           phi(2,2,4)
c      rcode = nf_var_par_access(cdfout,cdfid,nf_collective)
c      if(rcode.ne.0) write(*,*)"Error in par_access in write2d"
      rcode = nf_put_vara_real(cdfout,cdfid,startz,countz,phi)
      if(rcode.ne.0) write(*,*)"Error in ncwrite in write2d",rcode,my_pe
c      call MPI_BARRIER(comm,ier)
#endif
      return
      end

      subroutine collect_sfc(buf)

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

      integer np,ier,send_req,m
      integer itop,ibot,jtop,jbot,ix,jx,kx,jlocal,ilocal
      integer recv_req(MPI_STATUS_SIZE)

      real buf(nprocs,4),tbuf(4)
c
c
      ilocal = 2
      jlocal = 2
c
c slab is on the mastertask, so do straight copy
c
      if (iope) then
          buf(my_pe+1,1) = um(ilocal,jlocal,1)
          buf(my_pe+1,2) = vm(ilocal,jlocal,1)
          buf(my_pe+1,3) = p(ilocal,jlocal,1)
          buf(my_pe+1,4) = th(ilocal,jlocal,1)
c          buf(my_pe+1,5) = qv(ilocal,jlocal,1)
c          buf(my_pe+1,6) = botheat(ilocal,jlocal)
c          buf(my_pe+1,7) = botlhf(ilocal,jlocal)
c          buf(my_pe+1,8) = uistar(ilocal,jlocal)
c          buf(my_pe+1,9) = rhop(ilocal,jlocal,1)+rhob(1)
c
c need to transfer point from other processor
c
c first figure out processor and post receive
c we use a buffer to transfer the point values
c
          do np=1,nprocs
            if(np-1.ne.my_pe) then
              call MPI_RECV(tbuf, 4,
     &                MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req, ier)
              do m=1,4
                buf(np,m) = tbuf(m)
              enddo
            endif
          enddo
      else
c
c then see if current processor is the sender and post send if it is
c
          tbuf(1) = um(ilocal,jlocal,1)
          tbuf(2) = vm(ilocal,jlocal,1)
          tbuf(3) = p(ilocal,jlocal,1)
          tbuf(4) = th(ilocal,jlocal,1)
c          tbuf(5) = qv(ilocal,jlocal,1)
c          tbuf(6) = botheat(ilocal,jlocal)
c          tbuf(7) = botlhf(ilocal,jlocal)
c          tbuf(8) = uistar(ilocal,jlocal)
c          tbuf(9) = rhop(ilocal,jlocal,1)+rhob(1)
          call MPI_SEND(tbuf, 4, MPI_REAL, mastertask,
     &             mpitag_io, comm, ier)
c          call MPI_WAIT(send_req, status, ier)
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
c
c
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
     *          idimh1=idimh+1,jdimh1=jdimh+1,kdimh1=kdimh+1,
     *          kdimh2=kdimh+2,idimh2=idimh+2)

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
