c #ifdef RADIATE
#if defined(USEMPI)
      subroutine radbct(tstep)
c
c   $Id: radbc.F,v 1.5 1998/07/24 23:50:22 dwd Exp $
c
c This subroutine computes the radiation boundary condition for pressure
c along the model bottom
c
      include 'mpif.h'
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "communicate.inc"
c
      integer i,j,tstep,ierr
c
      integer hmany,istride,idist
      complex ww(idim,jdim),wo(idim,jdim)
c
c
c not needed if rigid lower boundary is used
c
      if(tstep.gt.0) then
c
c
        do j=1,jdim
          do i=1,idim
            ww(i,j) = cmplx(wm(i,j,kdim1))
          enddo
        enddo
c
c perform forward transform
c
c we use the normal transpose version of fftwnd (_n) to avoid mapping
c problems for the output array.  A better option might be to output
c the different array with dimensions jdim,idim from fftwnd
c
c
        call MPI_BARRIER(comm,ierr)
      hmany = jdim
      idist = idim
      istride = 1
      call fftw_f77(rfplan,hmany,ww,istride,idist,
     *              wo,istride,idist)
	call transpose1(wo,ww,comm,my_pe)
      call fftw_f77(rfplan,hmany,ww,istride,idist,
     *              wo,istride,idist)
c
c multiple amplitudes by N/k
c
         do j=1,jdim
           do i=1,idim
             wo(i,j) = bvft*wo(i,j)/rk(i,j)
           enddo
         enddo
c
c zero out mean
c
c         do j=1,jdim
c           wo(1,j)=cmplx(0.0)
c         enddo
c         if(blockstart(my_pe+1).eq.1) then
c           do i=1,idim
c             wo(i,1)=cmplx(0.0)
c           enddo
c         endif
c
c
c perform reverse transform
c

      call MPI_BARRIER(comm,ierr)
      call fftw_f77(iplan,hmany,wo,istride,idist,
     *              ww,istride,idist)
      call transpose1(ww,wo,comm,my_pe)
      call fftw_f77(iplan,hmany,wo,istride,idist,
     *              ww,istride,idist)
c
c
      else
c
         do j=1,jdim
           do i=1,idim
             ww(i,j) = 0.0
           enddo
         enddo
c
      endif
c
c
c update the pressure field
c divide by total number of points for backwards fft
c
      do j=1,jdim
        do i=1,idim
          p(i,j,kdim1) =  real(ww(i,j))/float(idim*jtotal)
        enddo
      enddo
c
c all done
c
      return
      end
c
#else
      subroutine radbct(tstep)
c
c   $Id: radbc.F,v 1.5 1998/07/24 23:50:22 dwd Exp $
c
c This subroutine computes the radiation boundary condition for pressure
c along the model bottom
c
#include "param.inc"
#include "olemf.inc"
#include "cyclic.inc"
#include "moddef.inc"
c
      integer i,j,tstep,nsign
c
      real ww(idim,jdim),wi(idim,jdim),wf(idim,jdim),
     *   wfi(idim,jdim)
      real pp(idim,jdim),pi(idim,jdim),pf(idim,jdim),
     *   pfi(idim,jdim)

c
c not needed if rigid lower boundary is used
c
      if(tstep.ne.0) then
         do 10 j=1,jdim
         do 10 i=1,idim
            ww(i,j) = wm(i,j,kdim1)
            wi(i,j) = 0.0
 10      continue
         nsign = 1
         call fft2d(ww,wi,wf,wfi,nstage,nsign,idim,jdim)
c         wf(idim/2+1,jdim/2+1) = 0.0
c         wfi(idim/2+1,jdim/2+1) = 0.0
         do 20 i=2,idim
         do 20 j=2,jdim
           pp(i,j) =  bvft*wf(i,j)/rk(i,j)
           pi(i,j) =  bvft*wfi(i,j)/rk(i,j)
 20      continue
c         do 21 i=1,idim
c             pp(i,1) = 0.0
c             pi(i,1) = 0.0
c 21      continue
c         do 22 j=1,jdim
c             pp(1,j) = 0.0
c             pi(1,j) = 0.0
c 22      continue
         nsign = -1
         call fft2d(pp,pi,pf,pfi,nstage,nsign,idim,jdim)
c
      else
c
         do 40 j=1,jdim
         do 40 i=1,idim
            pf(i,j) = 0.0
 40      continue
c
      endif
c
c
c update the pressure field
c
      do 60 i=1,idim
      do 60 j=1,jdim
         p(i,j,kdim1) = pf(i,j)
60    continue
c
c all done
c
      return
      end
#endif
