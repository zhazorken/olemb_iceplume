      subroutine volbot(hght)
#include "param.inc"
#include "finitv.inc"
#include "olemp.inc"
c
c This routine computes the finite element volume and cell boundary area
c given a 2-D terrain array
c
      real hght(0:idim1,0:jdim1)
      real tmpp(0:idim1,0:jdim1)
c      real way(0:idim1,0:jdim1,0:kdim1)
c      real frac(0:idim1,0:jdim1)
c      integer kdist(0:idim1,0:jdim1)
      integer i,j,k,kmin,junk,olev
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
            vay(i,j,k) = 1.0
            bay(i,j,k) = 1.0
            way(i,j,k) = 1.0
            frac(i,j) = 0.0
            kdist(i,j) = 0
          enddo
        enddo
      enddo
c
      do i=1,idim
        do j=1,jdim
          kdist(i,j) = int(hght(i,j)/dz)
          frac(i,j) = hght(i,j)/dz-kdist(i,j)
          if(frac(i,j).lt.0.05) then
            frac(i,j) = 0.0
          else if (frac(i,j).gt.0.95) then
            frac(i,j) = 0.0
            kdist(i,j) = kdist(i,j)+1
          endif
          kdist(i,j) = kdist(i,j)+1
          way(i,j,kdist(i,j)) = 1.-frac(i,j)
#if defined(HBOTOLD)
          ay(i,j,kdist(i,j)) = 1.-frac(i,j)
#endif
        enddo
c        write(*,*)"kdist,frac sample ",kdist(i,1),frac(i,1)
      enddo
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
c
c compute volume for v component as an average
c
      call imag(way)
c
      do i=1,idim
        do j=1,jdim
          olev = kdist(i,j)
          if(kdist(i,j).eq.kdist(i,j-1)) then
            if(way(i,j,olev).lt.way(i,j-1,olev)) then
              vay(i,j,olev) = way(i,j,olev)
              ay(i,j,olev) = way(i,j,olev)
            else
              vay(i,j,olev) = way(i,j-1,olev)
              ay(i,j,olev) = way(i,j-1,olev)
            endif
c          else if(kdist(i,j).gt.kdist(i,j-1)) then
c            vay(i,j,olev) = 1.0
c          else
c            vay(i,j,olev) = 0.0
          endif
            
c          vay(i,j,olev) = 
c     *         0.5*(way(i,j-1,olev)+way(i,j,olev))
c        if(kdist(i,j).gt.kdist(i,j-1)) then
c          do k=kdist(i,j-1),kdist(i,j)
c            vay(i,j,k) = 0.0
c          enddo
c        else if(kdist(i,j).gt.kdist(i,j+1)) then
c          do k=kdist(i,j+1),kdist(i,j)
c            vay(i,j+1,k) = 0.0
c          enddo
c        endif
        enddo
      enddo
c      do i=1,idim
c        do j=1,jdim
c          olev = kdist(i,j)
c          way(i,j,olev)= 0.5+way(i,j,olev)*0.5
c        enddo
c      enddo
c
      call imag(ax)
      call imag(ay)
      call imag(bay)
      call imag(vay)
      call imag(way)
      do i=0,idim
        do j=0,jdim
          if(kdist(i,j).gt.1) then
          do k=0,kdist(i,j)-1
            ax(i,j,k) = 0.0
            ay(i,j,k) = 0.0
            bay(i,j,k) = 0.0
c            vay(i,j,k) = 0.0
            ax(i+1,j,k) = 0.0
            ay(i,j+1,k) = 0.0
            bay(i,j+1,k) = 0.0
c            vay(i,j+1,k) = 0.0
          enddo
          endif
          bay(i,j,0) = 0.0
          bay(i,j+1,0) = 0.0
c          ax(i,j,0) = 0.0
c          ay(i,j,0) = 0.0
c          if(j.eq.jdim/2.and.my_pe.eq.3) then
c            write(*,*)"kdist, ax ", kdist(i,j),
c     *                ax(i,j,k),i
c          endif
        enddo
      enddo
c      read(*,*)junk
c      call imag(ax)
c      call imag(ay)
c      call imag(bay)
c      call imag(way)
c      call imag(vay)
      return
      end
