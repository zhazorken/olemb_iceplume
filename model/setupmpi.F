      subroutine setupmpi
c
c  $Id: setupmpi.F,v 1.7 1998/10/21 21:56:41 dwd Exp $
c

      include 'mpif.h'
#include "param.inc"
#include "olemp.inc"
#include "communicate.inc"
#include "fftw_f77_param.h"

c
c setup a 2-d cartesian layout
c      
      integer dims(2),coords(2)

      integer j,ierr,np,ik2dim,ncord
      integer j_global(jtotal)
      integer i_global(itotal)
      integer numblock,blocklen,stride,ndims
      logical periods(2), reorder
c
c Initiate MPI
c
c      call MPI_INIT(ierr)
c
#if defined(PRINTDEBUG)
      write(*,*)"Set up MPI"
#endif
      ndims = 2
      dims(1) = nprocx
      dims(2) = nprocy
      periods(1) = .true.
      periods(2) = .true.
      reorder = .true.
      mpitag_nshift = 1
      mpitag_sshift = 2
      mpitag_eshift = 3
      mpitag_wshift = 4
      mastertask = 0

      call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods,
     &                     reorder, comm, ierr)
c
c compute shifts for north and south
c
      call MPI_CART_SHIFT(comm,1,1,nbr_south,nbr_north,ierr)
      call MPI_CART_SHIFT(comm,0,1,nbr_west,nbr_east,ierr)
c
c get the process number
c
      call MPI_COMM_RANK(comm,my_pe,ierr)
c
c get the local cartesian group relative to my_pe (probably = my_pe)
c
      call MPI_CART_COORDS(comm,my_pe,2,coords,ierr)
c #if defined(PRINTDEBUG)
      write(*,*)"Coords created ",my_pe, coords(1),coords(2)
c #endif
      write(*,*)"r, e,w,s,n ",my_pe,nbr_east,nbr_west,
     *   nbr_south,nbr_north
c
c set flag for master task part of program
c
      if(my_pe.eq.mastertask) then
        iope = .true.
      else
        iope = .false.
      endif
c      write(*,*)"iope ",iope
c
c setup plan for fft used in poislk
c
c      call fftw_f77_create_plan(iplan,idim,FFTW_BACKWARD,
c     &    FFTW_ESTIMATE)
c      call fftw_f77_create_plan(rfplan,idim,FFTW_FORWARD,
c     &    FFTW_ESTIMATE)
c      call fftw_f77_create_plan(jplan,jtotal,FFTW_BACKWARD,
c     &    FFTW_ESTIMATE)
c
c      call fftw2d_f77_mpi_create_plan
c     &   (comm,iplan,idim,jtotal,FFTW_BACKWARD,
c     &    FFTW_ESTIMATE)
c
c set up plans for radiative b.c. routine
c
c      call fftw2d_f77_mpi_create_plan
c     &   (comm,rfplan,idim,jtotal,FFTW_FORWARD,
c     &    FFTW_ESTIMATE)
c
c      call fftw2d_f77_mpi_create_plan
c     &   (comm,rbplan,idim,jtotal,FFTW_BACKWARD,
c     &    FFTW_ESTIMATE)
#if defined(PRINTDEBUG)
c      write(*,*)"Done with fft plan setup"
#endif
c
c
c set up data blocks, assume fortran ordering with row changing fastest
c for example, for 32 x 32 array, divide into 32 x 8 chunks.  jglobal
c has the beginning, zero-based cartesian coordinate for each processor.
c jdim is the local maximum coordinate defined as jtotal/NPROC. jglobal is
c used for i/o functions and does not consider the ghost points 
c
c      do j=1,jdim
c        jglobal(j) = jdim*coords(2) + j 
c      enddo
c      do i=1,idim
c        iglobal(i) = idim*coords(1) + i
c      enddo
c
c figure out block starting address (basically an offset)
c and broadcast to each node
c
c i-axis
      do j=1,jdim
        j_global(j) = jdim*coords(2)+j
      enddo
      do i=1,idim
        i_global(i) = idim*coords(1)+i
      enddo
c
c      write(*,*)"my_pe, jglobal ",my_pe,(jglobal(j),j=1,jdim)
      do np=1,nprocs
        if(np-1.eq.my_pe) then
           i = i_global(1)
           j = j_global(1)
           ncordx = coords(1)
           ncordy = coords(2)
c #if defined(PRINTDEBUG)
c           write(*,*)"processor # ",my_pe," j_str,i_str ",j,i,
c     * "cordx,cordy ",ncordx,ncordy," jglb,iglb ",jglobal(jdim),
c     *  iglobal(idim)
c #endif
        endif
        call broadcast_scalar_int(i,np-1)
        call broadcast_scalar_int(j,np-1)
        call broadcast_scalar_int(ncordx,np-1)
        call broadcast_scalar_int(ncordy,np-1)
        call MPI_barrier(comm, ierr)
        iblkstart(np) = i
        iblkcord(np) = ncordx
        jblkstart(np) = j
        jblkcord(np) = ncordy
      enddo

      return
      end
