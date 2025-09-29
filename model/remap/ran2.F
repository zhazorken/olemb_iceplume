      real function ran2 (idum)
c
c   $Id: ran2.F,v 1.1.1.1 1998/07/20 22:03:48 eds Exp $
c
c     returns a uniform random deviate between 0.0 and 1.0.  Set
c     IDUM to any negative value to initialize or reinitialize
c     the sequence.
c
      save
      integer m,ia,ic
      real rm 
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      integer ir(97),idum,iff,j,mod,iy
      data iff/0/
      if (idum.lt.0.or.iff.eq.0) then
          iff=1
          idum=mod(ic-idum,m)
          do 10 j=1,97
              idum=mod(ia*idum+ic,m)
              ir(j)=idum
 10       continue
          idum=mod(ia*idum+ic,m)
          iy=idum
      end if
      j=1+(97*iy)/m
      if (j.gt.97.or.j.lt.1)then
         write(*,*)'Problem with random number'
         return
      endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
