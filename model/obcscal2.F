      subroutine obcscal(phi,phip)
c
c    $Id: obcscal.F,v 1.6 1998/10/21 21:56:38 eds Exp $
c
c this routine sets the lateral image points
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "moddef.inc"
#include "communicate.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real phip(0:idim1,0:jdim1,0:kdim1)
      real vbl,vbr
      integer i,j,k
c
#if !defined(RECIRC)
      if(blockstart(my_pe+1).eq.1) then
      do k=0,kdim1
        do i=1,idim
          j=1
c 
c compute the Doppler shifted phase speeds at the left and right 
c boundaries
c 
         vbl = v(i,1,k)-maxv/1.2
c
c         if(vbl.gt.vgeo) then
c           vbl = vgeo
c         endif
c
c
         j = 1
         if(vbl.gt.0.0) then
            phip(i,j,k) = phi(i,j,k)
         else
           phip(i,j,k) = phi(i,j,k) - 
     *            vbl*(phi(i,j+1,k)-phi(i,j,k))*dtdx
         endif
         phip(i,0,k) = phip(i,1,k)
       enddo
      enddo
      endif
#endif
      if(blockstart(my_pe+1)+jdim.ge.jtotal) then
      do k=0,kdim1
        do i=1,idim
          j = jdim1
c
c open on right side
c
c #if defined(OPENR)

         vbr = maxv/1.2 + v(i,jdim,k)
c         if(v(i,j-1,k).gt.0.0) then
c           vbr = v(i,j-1,k)
c	 endif
         if(vbr.lt.0.0) then
            phip(i,j,k) = phi(i,j,k)
         else
            phip(i,j,k) = phi(i,j,k) - 
     *            delt*vbr*(phi(i,j,k)-phi(i,j-1,k))*ody
c            phip(i,j,k) = phip(i,j-1,k)
         endif
c 
c #else
c
c zero gradient
c
c         phip(i,j,k) = phip(i,j-1,k)
c #endif
        enddo
      enddo
      endif
c
      return
      end
