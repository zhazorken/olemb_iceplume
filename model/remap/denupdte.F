      subroutine denupdte
c
c   $Id: denupdte.F,v 1.3 1998/10/21 21:56:36 dwd Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "dentable.inc"
c
      real s,t
      integer k
c
      integer i,kz,index,ii,kk,jj
      real y1,y2,y3,y4,ti,si
      real tmpt,tmps,ies
c
      do 175 jj=1,jdim
         do 180 ii=1,idim
            do 185 kk=1,kdim
               s=sal(ii,jj,kk)
               t=th(ii,jj,kk)
               kz=2*kk
c
               tmpt=(t-tmin)/dtemp
               tmps=(s-smin)/dsal
c
               i=ifix(tmpt)+1
               k=ifix(tmps)+1
c
               index=i+(k-1)*tsize+(kz-1)*tssize
               y1=rhotbl(index)
               y2=rhotbl(index+1)
               y3=rhotbl(index+1+tsize)
               y4=rhotbl(index+tsize)
               ti=tmpt-(i-1)
               si=tmps-(k-1)
               ies=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
               rhop(ii,jj,kk)=ies-rhobm1(kk)
 185        continue
c
c top and bottom B.C.
c
            rhop(ii,jj,kdim1) = rhop(ii,jj,kdim)
            rhop(ii,jj,0) = rhop(ii,jj,1)
 180     continue
c
c
c do ii=0 and ii=idim1 cyclic B.C.
c
c         do 190 kk=1,kdim
c            rhop(0,jj,kk)=rhop(idim,jj,kk)
c            rhop(idim1,jj,kk)=rhop(1,jj,kk)
c 190     continue
c
c corners
c
c         rhop(0,jj,kdim1) = rhop(0,jj,kdim)
c         rhop(0,jj,0) = rhop(0,jj,1)
c         rhop(idim1,jj,kdim1) = rhop(idim1,jj,kdim)
c         rhop(idim1,jj,0) = rhop(idim1,jj,1)
 175  continue
      call imag(rhop)
c
      return
      end
