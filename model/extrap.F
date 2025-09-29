      subroutine extrap(vc,vf)
#include "param.inc"
      integer idimh,jdimh,kdimh,idimh1,jdimh1,kdimh1
      parameter(idimh=idim/2,jdimh=jdim/2,kdimh=kdim/2,
     *          idimh1=idimh+1,jdimh1=jdimh+1,kdimh1=kdimh+1)

      real*4 vc(0:idimh1,0:jdimh1,0:kdimh1)
      real*4 vf(0:idim1,0:jdim1,0:kdim1)
      integer i,j,k

      do j=0,jdimh1
        do i=0,idimh1
          vc(i,j,kdimh1) = vc(i,j,kdimh)
        enddo
      enddo

      do k=1,kdimh1
        do j=1,jdimh1
          do i=1,idimh1
            vf(i*2-1,j*2-1,k*2-1) = vc(i,j,k)
          enddo
        enddo
      enddo
c
c update edge points
c
      call imag(vf)
c
      do k=1,kdim1,2
c
c first do the i-edge points
c
        do j=1,jdim,2
          do i=2,idim1,2
            vf(i,j,k) = 0.5*(vf(i+1,j,k)+vf(i-1,j,k))
          enddo
        enddo
c
c next do the j-edge points
c
        do j=2,jdim1,2
          do i=1,idim,2
            vf(i,j,k) = 0.5*(vf(i,j+1,k)+vf(i,j-1,k))
          enddo
        enddo
c
c finish with the middle point
c
        do j=2,jdim1,2
          do i=2,idim1,2
            vf(i,j,k) = 0.25*(vf(i-1,j-1,k)+vf(i+1,j-1,k)+
     *                       vf(i-1,j+1,k)+vf(i+1,j+1,k))
          enddo
        enddo

      enddo
      do k=2,kdim,2
        do i=1,idim
          do j=1,jdim
            vf(i,j,k)=0.5*(vf(i,j,k+1)+vf(i,j,k-1))
          enddo
        enddo
      enddo
c
      call imag(vf)
      return
      end
            
            
