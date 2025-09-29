      subroutine poislk(y1,gzb,gzt,im1,jm1,lk,k2,dx,dz)
c
c    $Id: poislk.F,v 1.3 1998/07/24 23:15:35 dwd Exp $
c
c k2=2*lk,  I think im1 = idim-1, jm1 = jdim-1, gzb and gzt are the top
c and bottom values of p (del2p = y1), and y1 is the right hand side of
c the poisson equation.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "fftw_f77_param.h"
#include "param.inc"
#include "communicate.inc"
c      implicit real*8 (a-h,o-z)
c
      integer im1,jm1,lk,k2,k1,k12,ierr,jstart
      integer hmany,istride,ostride,idist,odist
      complex y1(idim,jdim,kdim)
      complex yt(idim,jdim,kdim)
      real
     * a(idim,jdim,kdimt2)
     *,gzb(idim,jdim)
     *,gzt(idim,jdim)
      real dx,dz
c
      real p1,p2,dy,dzx2,dzy2,al
      integer n2,m2,nebot,netop,n2d,nmx,i,j,k,nnmx,nm4
c
      ostride = 1
      idist = 0
      odist = 0
      p1= 3.14159265358979d0
      p2=2.*p1
c  ********** begin specification of paremeters **********
c  ** n2,m2,lk,dx,dy,dz,nebot,netop,n2d,y1,gzb,gzt ***
      n2=im1
      m2=jm1
      dy=dx
#if defined(RADIATE)
c
c set boundary conditionto Dirichlet for radiation boundary condition
c on bottom
c
      nebot=0
#else
c
c otherwise, set boundary to Neumann
c
      nebot=1
#endif
#if defined(RADTOP)
      netop = 0
#else
      netop=1
#endif
      n2d=1
c
c   ***** get nmx=log2(lk)******
c
      nmx=0
      i=1
      do j=1,30
         i=i*2
         if(i.eq.lk) nmx=j
      enddo
c
c
c
      nnmx=.5000001*lk
      nm4=2*lk
      dzx2=(dz/dx)*(dz/dx)
      dzy2=(dz/dy)*(dz/dy)
c
c
c Set the bottom boundary
c 
      do j=1,jm1
        do i=1,im1
          y1(i,j,1)=y1(i,j,1)-gzb(i,j)
        enddo
      enddo
c
c set the top boundary
c
      do j=1,jm1
        do i=1,im1 
          y1(i,j,lk)=y1(i,j,lk)-gzt(i,j)
        enddo    
      enddo
c
c start scheme by computing the 2d fft in place using the mpi fftw
c routine
c
#if defined(USEMPI)
c setup plan for fft used in poislk
c
      hmany = jdim*kdim
      idist = idim
      istride = 1
      call fftw_f77(iplan,hmany,y1,istride,idist,
     *              yt,istride,idist)
c      call fft1d(y1,yt,-1,idim)
      call MPI_BARRIER(comm,ierr)
c      write(*,*)"yt after 1st fft ",yt(1,1,kdim-1),iplan,my_pe
      call transpose(yt,y1,comm,my_pe)
c      write(*,*)"y1 after trans ",y1(1,1,kdim-1),iplan,hmany,idist,my_pe
c      call fft1d(y1,yt,-1,idim)
      call fftw_f77(iplan,hmany,y1,istride,idist,
     *              yt,istride,idist)
c      write(*,*)"yt after 2nd fft ",yt(10,3,30),iplan,my_pe
#else
      call fftwnd_f77_mpi(iplan,many,y2,stride,idist,y2,ostride,odist)
#endif
c
c 
c
c      write(*,*)"poislk yt ",yt(1,1,kdim),my_pe
      k1=nnmx+1
      k2=nnmx+lk
      do k=1,nm4
        do j=1,m2
          do i=1,n2
            a(i,j,k) = 1.
          enddo
       enddo
      enddo
      do k12=k1,k2
c  **** al=.5 for Neumann bdry cond on top or bot, 1. otherwise ****
         al=1.
         if(k12.eq.k1) al=(2-nebot)*.5
         if(k12.eq.k2) al=(2-netop)*.5
         do j=0,m2-1
           do i=0,n2-1
#if defined(USEMPI)
             jstart = j+blockstart(my_pe+1)-1
#else
             jstart = j
#endif
             a(i+1,j+1,k12) =
     *        ((-2.+2.*cos(i*p2/n2))*dzx2 
     *        +(-2.+2.*cos(jstart*p2/jtotal))*dzy2-2.*al)
           enddo
         enddo
      enddo
c
c **** a(i1:i2,0,:) and a(i1:i2,:,0) are 1, not matter ,bdry pts ******
c
      call cycelim(a,yt,lk,nmx,nm4,n2,m2)
c
c *********
c
      do k=1,lk
        do j=1,m2
          do i=1,n2
#if defined(USEMPI)
            yt(i,j,k)=conjg(yt(i,j,k))/(n2*jtotal)
#else
            y1(i,j,k)=conjg(y1(i,j,k))/(n2*m2)
#endif
          enddo
        enddo
      enddo
c
#if defined(USEMPI)
c setup plan for fft used in poislk
c
c
c      call fft1d(yt,y1,-1,idim)
c      write(*,*)"y1 after round 2 fft ",y1(10,3,30),iplan,my_pe
      call fftw_f77(iplan,hmany,yt,istride,idist,
     *              y1,istride,idist)
      call MPI_BARRIER(comm,ierr)
      call transpose(y1,yt,comm,my_pe)
c      write(*,*)"yt after round 2 tran ",yt(10,3,30),iplan,my_pe
c      call fft1d(yt,y1,-1,idim)
c      write(*,*)"y1 after round 2 2nd fft ",y1(10,3,30),iplan,my_pe
      call fftw_f77(iplan,hmany,yt,istride,idist,
     *              y1,istride,idist)
#else
        call fftwnd_f77_mpi(iplan,many,y2,stride,idist,y2,ostride,odist)
#endif
c
c      call fft341(y1,lk,n2,m2)
c      call fft342(y1,lk,n2,m2)
c
c all done
c
      return
      end
c
c define cyclic reduction routine
c
c for this to work, must have a power of 2 vertical dimension
c
