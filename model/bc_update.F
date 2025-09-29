      subroutine bc_update (tstep)
c
c    $Id: bc_update.F,v 1.4 1998/08/04 17:25:32 eds Exp $
c
c routine to update the conditions for boundary forcing
c 
c tstep  =  time step number
c ivar   =  variable selected, currently only temperature
c 
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "moddef.inc"
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "cyclic.inc"
#include "olemtke.inc"
#include "refdate.inc"
c
      integer tstep,i,j,n1,n2
      integer htime,nn,dtime,datediff
      real trat,ran2,denom
      real phi2(itotal,jtotal)
      real randt(idim,jdim)
      logical ranon
      save ranon
#if defined(STOKES)
      integer k
      real stokabs
#endif
c
      data ranon/.true./
      data n1 /1/
      data n2 /2/
c 
c update heat flux term first
c
c first, make sure b.c. file has been entered
c
      if(bcon) then

         htime = int(float(tstep)*delt)
         if(datemode) then
            dtime = datediff(refdate, htime, bcref, tt(n2))
         else
            dtime = htime - tt(n2)
         endif
         if(dtime.gt.0) then
            do 100 nn=1,nsize
               if(n2.gt.timsiz) then
c
c set fields to last value model exceeds bc file time period
c
                  if(heaton) heatf = heat(n1)
                  if(swheaton) swheatf = swheat(n1)
c                  if(lheaton) latent = lheat(n1)
                  if(rainon) rainf = 0.001*rain(n1)/3600.
                  if(uston) ustr = ust(n1) 
                  if(vston) vstr = vst(n1) 
#if defined(STOKES)
                  if(wavelon) wavelen = wavel(n1)
                  if(wavehon) waveamp = waveh(n1)
                  if(angleon) waveang = angle(n1)
#endif
                  n2=n1
                  goto 110
               endif
               if(datemode) then
                  dtime = datediff(refdate, htime, bcref, tt(n2))
               else
                  dtime = htime - tt(n2)
               endif
               if(dtime.gt.0) then
                  n1 = n2
                  n2 = n2+1
               else
                  goto 105
               endif
 100        continue
         endif
 105     continue
c
c interpolate the b.c. value
c
         if(datemode) then
            dtime = datediff(refdate, htime, bcref, tt(n1))
         else
            dtime = htime - tt(n1)
         endif
         if(dtime.lt.0) then
            if(swheaton) swheatf = swheat(n1)
            if(heaton) heatf = heat(n1)
            if(lheaton) latent = lheat(n1)
            if(rainon) rainf = 0.001*rain(n1)/3600.
            if(uston) ustr = ust(n1)
            if(vston) vstr = vst(n1)
#if defined(STOKES)
            if(wavelon) wavelen = wavel(n1)
            if(wavehon) waveamp = waveh(n1)
            if(angleon) waveang = angle(n1)
#endif
         else
            if(datemode) then
               trat = float(datediff(refdate, htime, bcref, tt(n1)))/
     *              float(datediff(bcref, tt(n2), bcref, tt(n1)))
            else
               trat = float(htime-tt(n1))/     
     *              float(tt(n2)-tt(n1)) 
            endif
            if(heaton) heatf = heat(n1) + (heat(n2)-heat(n1))*trat
            if(swheaton) swheatf = swheat(n1) + 
     *           (swheat(n2)-swheat(n1))*trat
            if(lheaton) latent = lheat(n1) +
     *           (lheat(n2)-lheat(n1))*trat
c     
c rainfall rate is calculated for upcoming period
c
            if(rainon) rainf = 0.001*rain(n1)/(tt(n2)-tt(n1))
            if(uston) ustr = ust(n1) + (ust(n2)-ust(n1))*trat
            if(vston) vstr = vst(n1) + (vst(n2)-vst(n1))*trat
#if defined(STOKES)
            if(wavelon) wavelen = wavel(n1) + (wavel(n2)-wavel(n1))*trat
            if(wavehon) waveamp = waveh(n1) + (waveh(n2)-waveh(n1))*trat
            if(angleon) waveang = angle(n1) + (angle(n2)-angle(n1))*trat
#endif
         endif
c
 110     continue
         if(wavelen.le.5) wavelen = 5.0
c
c
c construct exponential stokes drift profiles from input surface value
c
c SET WIND and WAVES to zero
c
c         waveamp = 0.0
c         ustr = 0.0
c         vstr = 0.0
c
#if defined(STOKES)
         if(ustr.ne.0.0.or.vstr.ne.0.0) then
           denom = sqrt(ustr*ustr+vstr*vstr)
	 else
	   denom = 1.0
         endif
         do 120 k=1,kdim
            stokabs = ((3.14159*waveamp/wavelen)**2)*
     *           sqrt(1.561*wavelen)*
     *           exp(-12.57*(k*dz-0.5*dz)/wavelen)
#if defined(MONO) 
            us(kdim-k+1) =  sin(0.017453293*waveang)*stokabs
            vs(kdim-k+1)  = cos(0.017453293*waveang)*stokabs 
#else
            us(kdim-k+1) =  stokabs*ustr/denom
            vs(kdim-k+1)  = stokabs*vstr/denom
#endif
c     us(kdim-k+1) = waveamp*exp(-0.102*(k*dz-0.5*dz)) 
c     vs(kdim-k+1) = 0.0
 120     continue
         stokabs = ((3.14159*waveamp/wavelen)**2)*sqrt(1.561*wavelen)
#if defined(MONO)
         us(kdim1) =  sin(0.017453293*waveang)*stokabs
         vs(kdim1)  = cos(0.017453293*waveang)*stokabs
#else
         us(kdim1) =  stokabs*ustr/denom
         vs(kdim1)  = stokabs*vstr/denom
#endif
c     us(kdim1) = waveamp
c     vs(kdim1) = 0.0
#endif
#if defined(DEBUG) && defined(STOKES)
         write(*,*)'time ',tstep*delt
         write(*,*)'wind stress terms ',ustr,vstr
         write(*,*)'heat flux ',heatf
         write(*,*)'surface stokes drift components',us(kdim1),vs(kdim1)
#endif
c 
c end of time step
c change random cooling every 300 time steps (keep it random!)
c
      endif
      if(.not.bcuseran) then
         ranon = .false.
         do i=0,idim1
           do j=0,jdim1
             randn(i,j) = 1.0
           enddo
         enddo
      endif

      if(ranon)then
         if (((mod(tstep,300).eq.0) .and. (tstep .le. nrndom)) .and.
     *        bcuseran) then
            do 18 j=1,jtotal
               do 10 i=1,itotal
                  phi2(i,j)=(1.0-ranamp)+ranamp*ran2(nf)
 10            continue
c     
c     filter in x-direction
c     
               do 11 i=2,itotal-1
                  phi2(i,j)=(2.*phi2(i,j)+phi2(i+1,j)+
     *                 phi2(i-1,j))/4.
 11            continue
 18         continue
c
c
c filter in y-direction
c
            do 30 i=1,itotal
               do 31 j=2,jtotal-1
                  phi2(i,j)=(2.*phi2(i,j)+phi2(i,j+1)+
     *                 phi2(i,j-1))/4.
 31            continue
 30         continue
#if defined(USEMPI)
            call distrib_2d(randt,phi2)
#else
            do i=1,idim
              do j=1,jdim
                randt(i,j) = phi2(i,j)
              enddo
            enddo
#endif
c
            do i=1,idim
              do j=1,jdim
                randn(i,j) = randt(i,j)
              enddo
            enddo
            call imag2d(randn)
         else if(tstep.gt.nrndom) then
            ranon=.false.
            do i=0,idim1
              do j=0,jdim1
                randn(i,j) = 1.0
              enddo
            enddo
         endif
      endif
c
      return
      end
