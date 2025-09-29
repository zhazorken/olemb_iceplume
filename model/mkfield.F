c---------------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine mkfield(fld,idm,jdm,kdm,im,jm,km,dd,f,dz,icnt)

#include "param.inc"
      integer psize
      parameter (psize=500)

      integer idm,jdm,im,jm,icnt,i,k,j,kk,km,kdm
      real fld(idm,kdm,jdm),dd(icnt),f(icnt)
      logical flip
      real dz,tmp,ddz,slope,z,d(psize)


      do 10 i=1,icnt
         d(i)=dd(i)
 10   continue

C-----Bubble sort data
20    flip=.false.
      do 30 i=2,icnt
         if(d(i).lt.d(i-1))then
            flip=.true.
            tmp=d(i)
            d(i-1)=d(i)
            d(i)=tmp
            tmp=f(i)
            f(i-1)=f(i)
            f(i)=tmp
         endif
30    continue
      if(flip)goto 20

C-----Interpolate data to grid points
      if(dz*0.5.lt.d(1) .or. (km-0.5)*dz.gt.d(icnt))then
         write(*,*)' Field not made: Profile file does not cover ',
     &             'depths.',(km-0.5)*dz,d(icnt),dz*0.5,d(1)
         stop
      endif
c
c if using RADTOP then the profiles should be referenced from
c the bottom, not the top:
c
c    z   t
c    0   5
c   20   6
c
c  where z is now distance from the bottom and not depth.
c
c
#if defined(RADTOP)
      do 40 k=1,km
         z=(k-0.5)*dz
#else
      do 40 k=km,1,-1
         z=(km-k+0.5)*dz
#endif
         do 45 kk=1,icnt-1
            if(z.ge.d(kk) .and. z.lt.d(kk+1))goto 46
45       continue

C--------Interpolate temperature and remove background
46       ddz=z-d(kk)
         slope=(f(kk)-f(kk+1))/(d(kk)-d(kk+1))
         tmp=f(kk)+slope*ddz
         do 50 i=1,im
           do 60 j=1,jm
             fld(i,k,j)=tmp
60         continue
50       continue
40    continue
      return
      end
