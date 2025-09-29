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
      real tmpp(0:idim1,0:jdim1)
      real ang,yy
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
            ax(i,j,k)  = 1.0
            uax(i,j,k)  = 1.0
            ay(i,j,k)  = 1.0
            vay(i,j,k) = 1.0
            bay(i,j,k) = 1.0
            bax(i,j,k) = 1.0
            way(i,j,k) = 1.0
            frac(i,j)  = 0.0
            kdist(i,j) = 0
            dhdy(i,j) = 0.0
            dhdx(i,j) = 0.0
            calpha(i,j) = 1.0
            salpha(i,j) = 0.0
          enddo
        enddo
      enddo
c
      do i=1,idim
        do j=1,jdim
          kdist(i,j) = int(hght(i,j)/dz)+1
          dhdy(i,j) = (hght(i,j)-hght(i,j-1))/dy
c           yy = (blockstart(my_pe+1)+j-1)*dy - jtotal*dy*0.6
c
c           dhdy(i,j) = -2.*yy*h_mtn*a_mtn**2
c     *                      /(yy**2+a_mtn**2)**2
           if(i.eq.5) then
             write(*,*)"dhdy,j ",dhdy(i,j),blockstart(my_pe+1)+j-1
           endif
          dhdx(i,j) = 0.5*(hght(i+1,j)-hght(i-1,j))/dx
c           xx = i*dx - idim*dx*0.5

c           dhdx(i,j) = -2.*yy*h_mtn*a_mtn**2
c     *                      /(xx**2+yy**2+a_mtn**2)**2


          if(abs(dhdy(i,j)).lt.0.003) then
            dhdy(i,j) = 0.0
          endif
          if(abs(dhdx(i,j)).lt.0.003) then
            dhdx(i,j) = 0.0
          endif
          ang = atan(dhdy(i,j)+0.00001)
          calpha(i,j) = cos(ang)
          salpha(i,j) = sin(ang)
c          if(i.eq.10) then
c             write(*,*)"j ",blockstart(my_pe+1)+j,calpha(i,j),salpha(i,j),
c     *                 dhdy(i,j),hght(i,j),dy,dx
c          endif

c         this is wrong
c           frac(i,j) = hght(i,j)/dz-kdist(i,j)-1
c         this is right
          frac(i,j) = hght(i,j)/dz-kdist(i,j)+1

c           if (i.eq.4) then
c 
c             write(lunit,919) j,kdist(4,j),hght(4,j),frac(4,j)
c 919         format('mark 1 - j = ',i3,' kdist = ',i3,' hght = ',f5.2,
c      *       ' frac = ',f5.2)
c 
c           endif

          if(frac(i,j).le.0.5) then
            frac(i,j) = 0.0
          else if (frac(i,j).gt.0.5) then
            if(kdist(i,j-1).gt.kdist(i,j))then
c
c transition cell
c


c   mark i have commented this out for runs12
c
c              way(i,j,kdist(i,j)+1) = 0.5

              frac(i,j) = 0.0
              kdist(i,j) = kdist(i,j)+1
            else if(kdist(i,j-1).eq.kdist(i,j)) then
              frac(i,j) = 0.0
              kdist(i,j) = kdist(i,j)+1
            else
              frac(i,j) = 0.0
              kdist(i,j)=kdist(i,j)+1
            endif
            
          endif
           frac(i,j) = 1.-frac(i,j)

c           if (i.eq.4) then
c 
c             write(lunit,920) j,kdist(4,j),hght(4,j),frac(4,j)
c 920         format('mark 2 - j = ',i3,' kdist = ',i3,' hght = ',f5.2,
c      *       ' frac = ',f5.2)
c 
c           endif

        enddo
      enddo
c
c Next reverse vertical to get topography on the top
c
     do j=1,jdim
       do i=1,idim
         kdist(i,j) = kdim-kdist(i,j)
         dhdy(i,j) = -dhdy(i,j)
         dhdx(i,j) = -dhdx(i,j)
       enddo
     enddo
 
c
c do boundaries for kdist using real imag2d
c
      do j=1,jdim
        do i=1,idim
          tmpp(i,j) = kdist(i,j)
        enddo
      enddo
      call imag2di(tmpp)
      do j=0,jdim1
        do i=0,idim1
          kdist(i,j) = nint(tmpp(i,j))
        enddo
      enddo
c
c
      call imag2di(frac)
      call imag2di(dhdy)
      call imag2di(dhdx)
      call imag2di(calpha)
      call imag2di(salpha)
c
c compute volume for v component as an average
c
      call imag(way)
c
      do i=1,idim
        do j=1,jdim
          olev = kdist(i,j)
          if(kdist(i,j).eq.kdist(i,j-1)) then
            if(frac(i,j).lt.frac(i,j-1)) then
              vay(i,j,olev) = frac(i,j)
              ay(i,j,olev) = frac(i,j)
            else
              vay(i,j,olev) = frac(i,j-1)
              ay(i,j,olev) = frac(i,j-1)
            endif
          endif
          if(kdist(i,j).eq.kdist(i-1,j)) then
            if(frac(i,j).lt.frac(i-1,j)) then
              uax(i,j,olev) = frac(i,j)
              ax(i,j,olev) = frac(i,j)
            else
              uax(i,j,olev) = frac(i-1,j)
              ax(i,j,olev) = frac(i-1,j)
            endif
          endif
        enddo
      enddo

      call imag(ax)
      call imag(ay)

c      do i=1,idim
c        do j=1,jdim
c          olev = kdist(i,j)
c          way(i,j,olev) = (ax(i+1,j,olev)+ax(i,j,olev)+ay(i,j+1,olev)+
c     *                     ay(i,j,olev))*0.25
c        enddo
c      enddo
c
      call imag(bay)
      call imag(bax)
      call imag(vay)
      call imag(uax)
      call imag(way)
      do i=0,idim
        do j=0,jdim
          bax(i,j,0) = 0.0
          bay(i,j,0) = 0.0
          if(kdist(i,j).gt.1) then
          do k=0,kdist(i,j)-1
            ax(i,j,k) = 0.0
            bax(i,j,k) = 0.0
c            uax(i,j,k) = 0.0
            ay(i,j,k) = 0.0
            bay(i,j,k) = 0.0
c            vay(i,j,k) = 0.0
            ax(i+1,j,k) = 0.0
            bax(i+1,j,k) = 0.0
c            uax(i+1,j,k) = 0.0
            ay(i,j+1,k) = 0.0
            bay(i,j+1,k) = 0.0
c            vay(i,j+1,k) = 0.0
c            way(i,j,k) = 0.0
          enddo
          endif
        enddo
      enddo
      call imag(ax)
      call imag(bax)
      call imag(ay)
      call imag(bay)
      call imag(way)
      call imag(vay)
      call imag(uax)
      return
      end
