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
      integer i,j,klev
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

      integer np,ier,recv_req,send_req
      integer itop,ibot,jtop,jbot,ix,jx,kx,jlocal,ilocal
      integer status(MPI_STATUS_SIZE)

      real buf(11)
c
c first see if point in on mastertask processor
c
      jbot = jblkstart(mastertask+1)
      ibot = iblkstart(mastertask+1)
      jtop = jblkstart(mastertask+1)+jdim
      itop = iblkstart(mastertask+1)+idim
c
c make sure valid point
c
      if(jx.gt.jtotal.or.ix.gt.itotal.or.kx.gt.kdim) return
c
c      write(*,*)"start collecting"
      if(jx.ge.jbot.and.jx.lt.jtop.and.ix.ge.ibot.and.ix.lt.itop) then
c
c slab is on the mastertask, so do straight copy
c
c        write(*,*)"point on iope"
        if (iope) then
c          write(*,*)"point on iope",ix,jx
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
c              write(*,*)"collect recv ",np,jbot,jtop,ibot,itop,ix,jx
              call MPI_RECV(buf, 11,
     &                MPI_REAL,
     &                np-1,
     &                mpitag_io, comm, recv_req, ier)
              write(*,*)"collect got it, get out"
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

c          write(*,*)"collect send ",my_pe,ibot,jbot,itop,jtop,ix,jx
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
c      if(iope) write(*,*)"done collecting",my_pe
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
      real fld(0:idim1,0:jdim1,0:kdim1)
      real phi(idim,jdim,kdim)
      integer cdfin,cdfid,start(4),count(4),rcode
#if defined(USEMPI)
      integer np,ier,recv_req, send_req,i,k,nprec,nps
      integer npx,npy
      integer status(MPI_STATUS_SIZE)
#endif


c
      count(1) = idim
      count(2) = jdim
      count(3) = kdim
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1

c
c #if defined(USEMPI)
c #if defined(USEMPI_OLDNETCDF)
      if(iope) then
        do npx=1,nprocx
        do npy=1,nprocy
          start(1) = idim*(npx-1)+1
          start(2) = jdim*(npy-1)+1
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
            call MPI_SEND(phi(1,1,1),idim*jdim*kdim,MPI_REAL,
     $         nprec,
     $         mpitag_io,comm,ier)
          else
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            do k=1,kdim
              do j=1,jdim
                do i=1,idim
                  fld(i,j,k) = phi(i,j,k)
                enddo
              enddo
            enddo
          endif
        enddo
        enddo
      else
        call MPI_RECV(phi(1,1,1),idim*jdim*kdim,MPI_REAL,
     $                mastertask,mpitag_io,comm,status,ier)
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              fld(i,j,k) = phi(i,j,k)
            enddo
          enddo
        enddo
      endif
      call MPI_BARRIER(comm,ier)
c      write(*,*)"my pe ",my_pe,phi(1,1,1)
#if defined(NEWUSEMPI)
#if defined(USEMPI)
      start(1) = iblkstart(my_pe+1)
      start(2) = jblkstart(my_pe+1)
#endif
c      write(*,*)"reading input file on ",my_pe,start(1),start(2)
      call ncvgt(cdfin,cdfid,start,count,phi,rcode)
c      write(*,*)"finished reading on ",my_pe
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
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
c#if defined(USEMPI_OLDNETCDF)
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
          write(*,*)"Post Recv ",nprec
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
          write(*,*)"Output phi ",phi(1,1,1),nprec,startz(2)
        enddo
        enddo
      else
        write(*,*)"Post Sends ",my_pe
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
#if defined(NEWUSEMPI)
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
