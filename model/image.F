      subroutine image(phi)
c
c    $Id: image.F,v 1.2 1998/07/24 16:59:42 dwd Exp $
c
c this routine sets the lateral image points
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "moddef.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      integer j,k
c
      do j=1,jdim
         do k=0,kdim1
c 
c left boundary image point
c 
            phi(0,j,k) = phi(idim,j,k)
c 
c right boundary image point
c 
            phi(idim1,j,k) = phi(1,j,k)
c
         enddo
      enddo
c
      return
      end
c
#if !defined(SNGLSALT)
      subroutine imagedp(phi)
c
c    $Id: image.F,v 1.2 1998/07/24 16:59:42 dwd Exp $
c
c this routine sets the lateral image points
c
#include "param.inc"
#include "olemp.inc"
#include "moddef.inc"
c
      real*8 phi(0:idim1,0:jdim1,0:kdim1)
      integer j,k
c
      do j=1,jdim
         do k=0,kdim1
c 
c left boundary image point
c 
            phi(0,j,k) = phi(idim,j,k)
c 
c right boundary image point
c 
            phi(idim1,j,k) = phi(1,j,k)
c
         enddo
      enddo
c
      return
      end
#endif
