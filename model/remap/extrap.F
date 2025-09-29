      subroutine extrap(vc,vf)
#include "param.inc"
      real*4 vc(0:idim1,0:jdim1,0:kdim1)
      real*4 vf(0:idimd1,0:jdimd1,0:kdimd1)
      integer i,j,k


      do k=1,kdim1
        do j=1,jdim1
          do i=1,idim1
            vf(i*2-1,j*2-1,k*2-1) = vc(i,j,k)
          enddo
        enddo
      enddo
      do k=1,kdimd,2
c
c first do the i-edge points
c
        do j=1,jdimd,2
          do i=2,idimd1,2
            vf(i,j,k) = 0.5*(vf(i+1,j,k)+vf(i-1,j,k))
          enddo
        enddo
c
c next do the j-edge points
c
        do j=2,jdimd1,2
          do i=1,idimd,2
            vf(i,j,k) = 0.5*(vf(i,j+1,k)+vf(i,j-1,k))
          enddo
        enddo
c
c finish with the middle point
c
        do j=2,jdimd1,2
          do i=2,idimd1,2
            vf(i,j,k) = 0.25*(vf(i-1,j-1,k)+vf(i+1,j-1,k)+
     *                       vf(i-1,j+1,k)+vf(i+1,j+1,k))
          enddo
        enddo

      enddo
      do k=2,kdimd,2
        do i=1,idimd
          do j=1,jdimd
            vf(i,j,k)=0.5*(vf(i,j,k+1)+vf(i,j,k-1))
          enddo
        enddo
      enddo
c
      call imagd(vf)
      return
      end
            
            
