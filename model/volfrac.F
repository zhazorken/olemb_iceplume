      subroutine voltop(hght)
#include "param.inc"
#include "finitv.inc"
#include "olemp.inc"
c
c This routine computes the finite element volume and cell boundary area
c given a 2-D terrain array
c
      real hght(0:idim1,0:jdim1)
      real tmpp(0:idim1,0:jdim1),fraction
      real dhx,dhy,ddx,ddy,dys,dxs
c      real frac(0:idim1,0:jdim1)
c      integer kdist(0:idim1,0:jdim1)
      integer i,j,k,kmin,junk
c
c first calculate the grid level for the terrain at each point and the
c fractional distance (or volume) left in the grid cell
c
c set ktop to kdim1, this will need adjustment for bottom terrain
c
c
c
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            ax(i,j,k) = 1.0
            ay(i,j,k) = 1.0
            frac(i,j) = 1.0
            kdist(i,j) = 0
          enddo
        enddo
      enddo
c
      dys = dy*dy
      dxs = dx*dx
      kmin = kdim1
      do i=1,idim
        do j=1,jdim
          kdist(i,j) = int(hght(i,j)/dz)
c          frac(i,j) = hght(i,j)/dz-kdist(i,j)
          fraction = hght(i,j)/dz-kdist(i,j)
          frac(i,j) = 1.0 - fraction
          if(fraction.lt.0.5) then
            frac(i,j) = 1.0
          else if (fraction.ge.0.5) then
            frac(i,j) = 1.0
            kdist(i,j) = kdist(i,j)+1
          endif
          kdist(i,j) = kdim-kdist(i,j)
          if(kdist(i,j).eq. kdim1) then
             kdist(i,j) = kdim
          endif
          kmin = min(kmin,kdist(i,j))
          dhx = 0.5*(hght(i+1,j)-hght(i-1,j))
          dhy = 0.5*(hght(i,j+1)-hght(i,j-1))
          ddx = sqrt(dhx**2 + dxs)
          ddy = sqrt(dhy**2 + dys)
          harea(i,j) = ddx*ddy/(dx*dy)
c          if(i.eq.10) then
c             write(*,*)"harea,k,h ",harea(i,j),kdist(i,j),hght(i,j)
c          endif

        enddo
c        write(*,*)"kdist,frac sample ",kdist(i,1),frac(i,1)
      enddo
      call imag2d(harea)
c
c do boundaries for kdist using real imag2d
c
      do j=1,jdim
        do i=1,idim
          tmpp(i,j) = kdist(i,j)
        enddo
      enddo
      call imag2d(tmpp)
      do j=0,jdim1
        do i=0,idim1
          kdist(i,j) = nint(tmpp(i,j))
        enddo
      enddo
c
c
      call imag2d(frac)
      do i=1,idim
      do j=1,jdim
c
c first do ax
c
c if k levels are different, mask off vertical wall
c
        if(kdist(i,j).gt.kdist(i-1,j)) then
          do k = kdist(i-1,j)+1,kdist(i,j)
            ax(i,j,k) = 0.0
          enddo
          ax(i,j,kdist(i-1,j)) = frac(i-1,j)
c          ax(i,j,kdist(i,j)) = (frac(i,j))
c
c if k levels are the same, compute minimum area
c
        else if(kdist(i,j).eq.kdist(i-1,j)) then
          if(frac(i,j).lt.frac(i-1,j)) then
            ax(i,j,kdist(i,j)) = frac(i,j)
          else
            ax(i,j,kdist(i,j)) = frac(i-1,j)
          endif
        else
c
c do other vertical wall
c
          do k=kdist(i,j)+1,kdist(i-1,j)
            ax(i,j,k) = 0.0
          enddo
          ax(i,j,kdist(i,j)) = frac(i,j)
c          ax(i,j,kdist(i-1,j)) = frac(i-1,j)
c
        endif
c
c next do ay
c
c if k levels are different, mask off vertical wall
c
        if(kdist(i,j).gt.kdist(i,j-1)) then
          do k = kdist(i,j-1)+1,kdist(i,j)
            ay(i,j,k) = 0.0
          enddo
          ay(i,j,kdist(i,j-1)) = frac(i,j-1)
c
c if k levels are the same, compute minimum area
c
        else if(kdist(i,j).eq.kdist(i,j-1)) then
          if(frac(i,j).lt.frac(i,j-1)) then
            ay(i,j,kdist(i,j)) = frac(i,j)
          else
            ay(i,j,kdist(i,j)) = frac(i,j-1)
          endif
        else
c
c do other vertical wall
c
          do k=kdist(i,j)+1,kdist(i,j-1)
            ay(i,j,k) = 0.0
          enddo
          ay(i,j,kdist(i,j)) = frac(i,j)
        endif

      enddo
      enddo
c
c update boundaries
c
      call imag(ax)
      call imag(ay)
c
c
      do i=1,idim
        do j=1,jdim
          do k=kdist(i,j)+1,kdim1
            ax(i,j,k) = 0.0
            ay(i,j,k) = 0.0
            ax(i+1,j,k) = 0.0
            ay(i,j+1,k) = 0.0
          enddo
          k = kdist(i,j)
          if(j.eq.jdim/2) then
c            write(*,*)"kdist, ax ", kdist(i,j),
c     *                ax(i,j,k)
          endif
        enddo
      enddo
      call imag(ax)
      call imag(ay)
      return
      end
