      subroutine voltop(hght)
#include "param.inc"
#include "finitv.inc"
#include "olemp.inc"
#include "communicate.inc"
c
c This routine computes the finite element volume and cell boundary area
c given a 2-D terrain array
c
      real hght(0:idim1,0:jdim1)
      real tmpp(0:idim1,0:jdim1),fraction
c      real way(0:idim1,0:jdim1,0:kdim1)
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
            vay(i,j,k) = 1.0
            bay(i,j,k) = 1.0
            bax(i,j,k) = 1.0
            way(i,j,k) = 1.0
            frac(i,j) = 1.0
            kdist(i,j) = 0
            dhdx(i,j) = 0.0
            dhdy(i,j) = 0.0
          enddo
        enddo
      enddo
c
      dys = dy*dy
      kmin = kdim1
      do i=1,idim
        dxs = dx(i)*dx(i)
        do j=1,jdim
          kdist(i,j) = int(hght(i,j)/dz)
          fraction = hght(i,j)/dz-kdist(i,j)
          if(fraction.lt.0.5) then
            frac(i,j) = 1.0
          else if (fraction.ge.0.5) then
            frac(i,j) = 1.0
            kdist(i,j) = kdist(i,j)+1
          endif
c          frac(i,j) = 1.0 - fraction
          kdist(i,j) = kdim-kdist(i,j)
          if(kdist(i,j).eq. kdim1) then
             kdist(i,j) = kdim
          endif
          kmin = min(kmin,kdist(i,j))
          if(kdist(i,j).eq.1) kdist(i,j) = 0
c
c minus sign for upside down topography
c
          if(hght(i+1,j).ne.hght(i,j).and.hght(i,j).ne.hght(i-1,j)) then
            dhdx(i,j) = -0.5*(hght(i+1,j)-hght(i-1,j))/dx(i)
          endif
          if(hght(i,j+1).ne.hght(i,j).and.hght(i,j).ne.hght(i,j-1)) then
            dhdy(i,j) = -0.5*(hght(i,j+1)-hght(i,j-1))/dy
          endif
          ddx = sqrt(dhx**2 + dxs)
          ddy = sqrt(dhy**2 + dys)
          harea(i,j) = ddx*ddy/(dx(i)*dy)
c
c set face to 1.0 for lead case
c
c          harea(i,j) = 1.0
c          if(i.eq.10) then
c             write(*,*)"harea,k,h ",harea(i,j),kdist(i,j),hght(i,j)
c          endif

        enddo
c        write(*,*)"kdist,frac sample ",kdist(i,1),frac(i,1)
      enddo
      call imag2d(harea)
      call imag2d(dhdy)
      call imag2d(dhdx)
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
c
      call imag(ax)
      call imag(ay)
      call imag(bay)
      call imag(bax)
#if defined(CHANNEL)
      do i=0,idim
        do k=0,kdim
          jj = jblkstart(my_pe+1)+jdim
          if(jj.ge.jtotal) then
             ay(i,jdim1,k) = 0.0
             bay(i,jdim1,k) = 0.0
          endif
          if(jblkstart(my_pe+1).eq.1) then
             ay(i,1,k) = 0.0
             bay(i,1,k) = 0.0
          endif
        enddo
      enddo
#endif
c
      do i=0,idim
        do j=0,jdim
          if(kdist(i,j).lt.kdim) then
          do k=kdist(i,j)+1,kdim1
            ax(i,j,k) = 0.0
            ay(i,j,k) = 0.0
            bay(i,j,k) = 0.0
            bax(i,j,k) = 0.0
            ax(i+1,j,k) = 0.0
            bax(i+1,j,k) = 0.0
            ay(i,j+1,k) = 0.0
            bay(i,j+1,k) = 0.0
            way(i,j,k) = 0.0
          enddo
          endif
c          if(j.eq.jdim/2) then
c            write(*,*)"kdist, ax ", kdist(i,j),
c     *                ax(i,j,k)
c          endif
        enddo
      enddo
      call imag(ax)
      call imag(bax)
c      call imag(ay)
c      call imag(bay)
      call imag(way)
      return
      end
