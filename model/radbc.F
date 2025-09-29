c #ifdef RADIATE
#if defined(USEMPI)
      subroutine radbc(tstep)
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
            ww(i,j) = cmplx(wm(i,j,1),0.0)
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
             wo(i,j) = cmplx(-1.*bvf*real(wo(i,j))/rk(i,j),
     $                -1.*bvf*imag(wo(i,j))/rk(i,j))
c             wo(i,j) = -1.*bvf*wo(i,j)/rk(i,j)
           enddo
         enddo
c
c zero out mean
c
c #if defined(ZERO)
         do j=1,jdim
           wo(1,j)=cmplx(0.0,0.0)
         enddo
         if(blockstart(my_pe+1).eq.1) then
           do i=1,idim
             wo(i,1)=cmplx(0.0,0.0)
           enddo
         endif
c #endif
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
          p(i,j,0) =  real(ww(i,j))/float(idim*jtotal)
        enddo
      enddo
c
c all done
c
      return
      end
c
#else
      subroutine radbc(tstep)
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
            ww(i,j) = wm(i,j,1)
            wi(i,j) = 0.0
 10      continue
         nsign = 1
         call fft2d(ww,wi,wf,wfi,nstage,nsign,idim,jdim)
         wf(idim/2+1,jdim/2+1) = 0.0
         wfi(idim/2+1,jdim/2+1) = 0.0
         do 20 i=2,idim
         do 20 j=2,jdim
           pp(i,j) =  -bvf*wf(i,j)/rk(i,j)
           pi(i,j) =  -bvf*wfi(i,j)/rk(i,j)
 20      continue
         do 21 i=1,idim
             pp(i,1) = 0.0
             pi(i,1) = 0.0
 21      continue
         do 22 j=1,jdim
             pp(1,j) = 0.0
             pi(1,j) = 0.0
 22      continue
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
         p(i,j,0) = pf(i,j)
60    continue
c
c all done
c
      return
      end
c
c the following subroutine computes the two-dimensional fast
c fourier transform
c
      subroutine fft2d(ww1,wi1,wf1,wfi1,nstg,nsn,iimaxx,jjmaxx)
c
c
#include "param.inc"
c
      real data(ijdimt2)
      real ww1(idim,jdim),wf1(idim,jdim)
      real wi1(idim,jdim),wfi1(idim,jdim)
      integer ii,i,j,nstg,nsn,iimaxx,jjmaxx,nij(2)

      ii=1
      nij(1) = iimaxx
      nij(2) = jjmaxx

      do 10 j=1,jjmaxx
      do 10 i=1,iimaxx
          data(ii) = ww1(i,j)
          data(ii+1) = wi1(i,j)
          ii = ii+2
10    continue
      call fourn(data,nij,2,nsn)
      if(nsn.lt.0) then
          ii = 1
          do 15 j=1,jjmaxx
          do 15 i=1,iimaxx
#ifdef cray
             wf1(i,j) = data(ii)/float(iimaxx*jjmaxx)
             wfi1(i,j) = data(ii+1)/float(iimaxx*jjmaxx)
#else
             wf1(i,j) = data(ii)/dble(float(iimaxx*jjmaxx))
             wfi1(i,j) = data(ii+1)/dble(float(iimaxx*jjmaxx))
#endif
             ii = ii+2
15        continue
      else
          ii = 1
          do 16 j=1,jjmaxx
          do 16 i=1,iimaxx
             wf1(i,j) = data(ii)
             wfi1(i,j) = data(ii+1)
             ii = ii+2
16        continue
      endif

      return
      end

#ifdef cray
        subroutine fourn(data,nn,ndim,isign)
        integer ndim
        real wr,wi,wpr,wpi,wtemp,theta
        real tempr,tempi
        integer nn(ndim)
        real data(ndim)
        integer ntot,idim,nprev,ip1,ip2,ip3,n,nrem,i2rev,i3rev
        integer i1,i2,i3,ibit,ifp1,ifp2,isign,k1,k2
        ntot=1
        do 9 idim=1,ndim
                ntot=ntot*nn(idim)
 9      continue
        nprev=1
        do 16 idim=1,ndim
           n=nn(idim)
           nrem=ntot/(n*nprev)
           ip1=2*nprev
           ip2=ip1*n
           ip3=ip2*nrem
           i2rev=1
           do 12 i2=1,ip2,ip1
                if(i2.lt.i2rev)then
                    do 11 i1=i2,i2+ip1-2,2
                        do 10 i3=i1,ip3,ip2
                                i3rev=i2rev+i3-i2
                                tempr=data(i3)
                                tempi=data(i3+1)
                                data(i3)=data(i3rev)
                                data(i3+1)=data(i3rev+1)
                                data(i3rev)=tempr
                                data(i3rev+1)=tempi
 10                            continue
 11                         continue
                endif
                ibit=ip2/2
1               if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
                        i2rev=i2rev-ibit
                        ibit=ibit/2
                go to 1
                endif
                i2rev=i2rev+ibit
 12             continue
                ifp1=ip1
2               if(ifp1.lt.ip2)then
                        ifp2=2*ifp1
                        theta=isign*6.28318530717959/(ifp2/ip1)
                        wpr=-2.0*sin(0.5*theta)**2
                        wpi=sin(theta)
                        wr=1.0
                        wi=0.0
                        do 15 i3=1,ifp1,ip1
                        do 14 i1=i3,i3+ip1-2,2
                           do 13 i2=i1,ip3,ifp2
                                k1=i2
                                k2=k1+ifp1
                                tempr=
     *                  wr*data(k2)-wi*data(k2+1)
                                tempi=
     *                  wr*data(k2+1)+wi*data(k2)
                                data(k2)=data(k1)-tempr
                                data(k2+1)=data(k1+1)-tempi
                                data(k1)=data(k1)+tempr
                                data(k1+1)=data(k1+1)+tempi
 13                          continue
 14                     continue
                        wtemp=wr
                        wr=wr*wpr-wi*wpi+wr
                        wi=wi*wpr+wtemp*wpi+wi
 15             continue
                ifp1=ifp2
                go to 2
                endif
                nprev=n*nprev
 16     continue
        return
        end
#else

        subroutine fourn(data,nn,ndim,isign)
        integer ndim
        real*8 wr,wi,wpr,wpi,wtemp,theta
        real tempr,tempi
        integer nn(ndim)
        real data(ndim)
        integer ntot,idim,nprev,ip1,ip2,ip3,n,nrem,i2rev,i3rev
        integer i1,i2,i3,ibit,ifp1,ifp2,isign,k1,k2
        ntot=1
        do 9 idim=1,ndim
                ntot=ntot*nn(idim)
 9      continue
        nprev=1
        do 16 idim=1,ndim
           n=nn(idim)
           nrem=ntot/(n*nprev)
           ip1=2*nprev
           ip2=ip1*n
           ip3=ip2*nrem
           i2rev=1
           do 12 i2=1,ip2,ip1
                if(i2.lt.i2rev)then
                    do 11 i1=i2,i2+ip1-2,2
                        do 10 i3=i1,ip3,ip2
                                i3rev=i2rev+i3-i2
                                tempr=data(i3)
                                tempi=data(i3+1)
                                data(i3)=data(i3rev)
                                data(i3+1)=data(i3rev+1)
                                data(i3rev)=tempr
                                data(i3rev+1)=tempi
 10                            continue
 11                         continue
                endif
                ibit=ip2/2
1               if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
                        i2rev=i2rev-ibit
                        ibit=ibit/2
                go to 1
                endif
                i2rev=i2rev+ibit
 12             continue
                ifp1=ip1
2               if(ifp1.lt.ip2)then
                        ifp2=2*ifp1
                        theta=isign*6.28318530717959d0/(ifp2/ip1)
                        wpr=-2.d0*dsin(0.5d0*theta)**2
                        wpi=dsin(theta)
                        wr=1.d0
                        wi=0.d0
                        do 15 i3=1,ifp1,ip1
                        do 14 i1=i3,i3+ip1-2,2
                           do 13 i2=i1,ip3,ifp2
                                k1=i2
                                k2=k1+ifp1
                                tempr=
     *                  sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                                tempi=
     *                  sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                                data(k2)=data(k1)-tempr
                                data(k2+1)=data(k1+1)-tempi
                                data(k1)=data(k1)+tempr
                                data(k1+1)=data(k1+1)+tempi
 13                          continue
 14                     continue
                        wtemp=wr
                        wr=wr*wpr-wi*wpi+wr
                        wi=wi*wpr+wtemp*wpi+wi
 15             continue
                ifp1=ifp2
                go to 2
                endif
                nprev=n*nprev
 16     continue
        return
        end
#endif /* cray */
c #endif /* RADIATE */
#endif /* cmf */
