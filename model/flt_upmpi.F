      subroutine flt_update(n)
      integer n
c
c     $Id: flt_update.F,v 1.3 1998/07/24 23:15:28 dwd Exp $
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "float.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "communicate.inc"
c
      integer l,i,j,k,nn
      real tmpxu,tmpxw,tmpyu,tmpyv,tmpzu,tmpzw,xi,yi,zi
      real ui,u1,u2,u3,u4,u5,u6,u7,u8
      real vi,v1,v2,v3,v4,v5,v6,v7,v8
      real wi,w1,w2,w3,w4,w5,w6,w7,w8
      real xsize,ysize,zsize
      real ufl(0:idim1,0:jtotal1),vfl(0:idim1,0:jtotal1)
      real uff(0:idim1,jtotal),vff(0:idim1,jtotal)
c
c first gather all horizontal velocity data from k=kdim
c and fix boundaries
c
      write(*,*)"Collecting u,v ",my_pe
      call collect_xyf(u,uff,kdim)
      call collect_xyf(v,vff,kdim)
c
      if(iope) then
      do j=1,jtotal
        do i=0,idim1
          ufl(i,j) = uff(i,j)
          vfl(i,j) = vff(i,j)
        enddo
      enddo
      do i=0,idim1
         ufl(i,0) = ufl(i,jtotal)
         vfl(i,0) = vfl(i,jtotal)
         ufl(i,jtotal1) = ufl(i,1)
         vfl(i,jtotal1) = vfl(i,1)
      enddo
c
c
      write(*,*)"Moving floats"
      xsize=dx*idim
      ysize=dy*jtotal
c
      do 10 l=1,ndesc
         if(n.ge.fltstrt(l) .and. n.le.fltstop(l))then
            do 11 nn=1,fltdescsize(l)

            tmpxu = fltx(l,nn)/dx+1.0
            tmpxw = fltx(l,nn)/dx+0.5
            tmpyu = flty(l,nn)/dy+0.5
            tmpyv = flty(l,nn)/dy+1.0
c
c u velocity
c
            i = ifix(tmpxu)
            j = ifix(tmpyu)
            xi = tmpxu - i
            yi = tmpyu - j
c
               u1=ufl(i,j)
               u2=ufl(i+1,j)
               u5=ufl(i,j+1)
               u6=ufl(i+1,j+1)
               ui=u1+(u2-u1)*xi+(u5-u1)*yi+(u6-u5+u1-u2)*xi*yi
c
c v velocity
c
            i = ifix(tmpxw)
            j = ifix(tmpyv)
            xi = tmpxw - i
            yi = tmpyv - j
c
               v1=vfl(i,j)
               v2=vfl(i+1,j)
               v5=vfl(i,j+1)
               v6=vfl(i+1,j+1)
               vi=v1+(v2-v1)*xi+(v5-v1)*yi+(v6-v5+v1-v2)*xi*yi
c
c
            if(nn.eq.1) write(*,*)"ui, vi ",ui,vi
            fltx(l,nn)=fltx(l,nn)+ui*delt
            flty(l,nn)=flty(l,nn)+vi*delt
            fltz(l,nn) = 1.0
c
c test for wrap and bounds
c
            if(fltx(l,nn).lt.0) fltx(l,nn) = xsize+fltx(l,nn)
            if(fltx(l,nn).gt.xsize) fltx(l,nn) = fltx(l,nn)-xsize
            if(flty(l,nn).lt.0) flty(l,nn) = ysize+flty(l,nn)
            if(flty(l,nn).gt.ysize) flty(l,nn) = flty(l,nn)-ysize
 11         continue
         endif
 10   continue
c
c iope endif
c
      endif
      return
      end
c
c collection routine
c
      subroutine collect_xyf(temp,phi,klev)

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
      real phi(0:idim1,jtotal)
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
          call MPI_IRECV(phi(0,blockstart(np)), idim2*jdim, mpi_real,
     &                np-1,
     &                mpitag_io, MPI_COMM_WORLD, recv_req(np), ier)
        enddo
      endif
c
c now send data chunks from each node
c
      call MPI_ISEND(temp(0,1,klev),idim2*jdim,mpi_real, mastertask,
     &             mpitag_io, MPI_COMM_WORLD, send_req, ier)
      call MPI_WAIT(send_req, status, ier)
c
c make sure everything is communicated before moving on
c
      if(iope) then
        call MPI_WAITALL(nprocs, recv_req, recv_status, ier)
      endif
c
c
      call MPI_BARRIER(MPI_COMM_WORLD,ier)

      return
      end
