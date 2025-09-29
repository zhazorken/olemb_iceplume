      subroutine volbot(hght)
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
      real ax2(0:idim1,0:kdim1)
      real hd2(0:idim1,0:kdim1)
      real dx2,dz2,frc1,mslp1,mslp2
c      real way(0:idim1,0:jdim1,0:kdim1)
c      real frac(0:idim1,0:jdim1)
c      integer kdist(0:idim1,0:jdim1)
      integer ind1,kmin,junk,olev
      integer ibeg,iend,i,j,k,flip
c
c first calculate the grid level for the terrain at each point and the
c fractional distance (or volume) left in the grid cell
c
c set ktop to kdim1, this will need adjustment for bottom terrain
c
c
c slope greater that one is steep, do horizontal
c
      do k=0,kdim1
        do i=0,idim1
          ax2(i,k) = 1.0
        enddo
      enddo 
c
c start with 45 degree slope to see if ax is correct
c
c undercut positive
c
c      mslope = 1.0
c      mslope = 1.0
      mslope = 8.0
      mslp1 = 5
      mslp2 = 3
      mslope = 5./3.
c
      mslope = 2.0
c      mslope = 16.0
      yinter = 0.0
c
c all overcut is negative
c
c      mslope = -8.0
c      mslope = -2.0
c      mslope = -0.5
c      yinter = kdim/(-mslope)

      if(abs(mslope).gt.0.0) then
        dz2 = dz
        dx2 = dxx
        ibeg = iblkstart(my_pe+1)
        iend = ibeg+idim-1
        flip = 1
        do k=0,kdim1
c
c terrain coordinate in absolute units
cc
c          yinter = 1.0
          if(k.lt.(kdim-5))then
c
          hx(k) =  yinter+(k-0.5)/mslope
c          hx(k) =  yinter+(k-0.5)*mslp2/mslp1
          else
          hx(k) = yinter+((kdim-5)-0.5)/mslope
c          hx(k) = yinter+((kdim-5)-0.5)*mslp2/mslp1
          endif
c
c overcut bottom fix miniumum of 8 vertical levels
c
          nvert = nint(6/abs(mslope))+1
          if(iope) write(*,*)"number of levels nvert ",nvert
          if(mslope.lt.0.0) then
            if(hx(k).gt.yinter-nvert) hx(k) = yinter-nvert
          endif
c sin shape
c
cc         hx(k) = kdim+5.0 + 0.5*kdim*sin((kdim+k)*1.57/kdim)
c         hx(k) = 205.-200*sin(k*1.57/kdim)
         
c
c constant no slope
c
          if(int(hx(k)).le.3.) then
             hx(k) = 3.
          endif
c
c          hx(k) = 5.
c          hx(k) = hx(k)+4.
c          if(mod(k,4).eq.0) then
c             flip = -flip
c          endif
c          hx(k) = 8.+(flip+1)*0.5
c          angwgt = 1.0

 
c
          ind1 = int(hx(k))
          if(iope)write(*,*)"k,ind,ibeg,iend,hx ",k,ind1,ibeg,iend,hx(k)
c
          if(ind1.ge.ibeg.and.ind1.le.iend) then
            i = ind1-ibeg+1
            do ii=0,i
              ax2(ii,k) = 1.0
              hd2(ii,k) = hx(k)-ii-ibeg
c              hd2(ii,k) = hx(k)-ii
            enddo
c            ax2(i+1,k) = hx(k)-ind1
            ax2(i,k) = hx(k)-ind1
            if(iope)write(*,*)"i,k,ind,ibeg,iend,hx,ax2 ",
     *       i,k,ind1,ibeg,iend,hx(k),ax2(i,k)

          else if(ind1.lt.ibeg) then
              do i=0,idim1
                ax2(i,k) = 0.0
                hd2(i,k) = 0.0
              enddo
          else if(ind1.gt.iend) then
            do i=0,idim1
              ax2(i,k) = 1.0
              hd2(i,k) = hx(k)-iend+(idim-i)
            enddo
          endif
        enddo
      endif
c
c do vertical
c
c
        ibeg = iblkstart(my_pe+1)
c        ibeg = iblkstart(my_pe+1)*2-1
        iend = ibeg+idim-1
c        iend = ibeg+2*idim-1
        

      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            ax(i,j,k) = 1.0
            bax(i,j,k) = 0.0
            ay(i,j,k) = 1.0
            vay(i,j,k) = 1.0
            bay(i,j,k) = 0.0
            way(i,j,k) = 1.0
            kdist(i,j) = 1
          enddo
        enddo
      enddo
c
c      call imag(ax)
c      call imag(ay)
c      call imag(way)
      do i=0,idim1
        do j=0,jdim1
          do k=0,kdim1
c            ax(i,j,k) = ax2(i,k)
            if(ax2(i,k).lt.1.0.and.ax2(i,k).gt.0.0) then
              ay(i,j,k) = 1.0
              ax(i,j,k) = 0.0
              ax(i+1,j,k) = 1.0
              way(i,j,k) = 0.0
            else if(ax2(i,k).eq.0.0) then
              ay(i,j,k) = 1.0
              ay(i,j+1,k) = 1.0
              ax(i,j,k) = 1.0
              ax(i+1,j,k) = 1.0
              way(i,j,k) = 1.0
              way(i,j,k+1) = 1.0
            else
              ay(i,j,k) = 0.0
              ax(i+1,j,k) = 0.0
              ax(i,j,k) = 0.0
              way(i,j,k) = 0.0
            endif

              
            hdx(i,j,k) = hd2(i,k)
c
c            way(i-1,j,k-1) = 0.25*(ax2(ii,kk-1)+ax2(ii+1,kk-1)+
c     *                   ax2(ii,kk)+ax2(ii+1,kk))
c            ay(i-1,j,k-1) = 0.25*(ax2(ii,kk)+ax2(ii+1,kk)+
c     *                  ax2(ii,kk+1)+ax2(ii+1,kk+1))
cc
c            hdx(i-1,j,k-1) = 0.25*(hd2(ii-1,kk)+hd2(ii,kk)+
c     *                  +hd2(ii-1,kk+1)+hd2(ii,kk+1))

c
c test for 45 degree slope
c
c          if(kdist(i,j).gt.0) then
c          if(kdist(i,j).le.kdim) then
c            way(i,j,kdist(i,j))=0.5
c            vay(i,j,kdist(i,j))=0.0
c            ax(i,j,kdist(i,j))=0.5
c            bax(i,j,kdist(i,j))=0.0
c            ay(i,j,kdist(i,j))=0.875
c            ay(i,j,kdist(i,j)-1)=0.125
c          else
c            way(i,j,kdist(i,j))=0.0
c            way(i,j,kdist(i,j)+1)=0.0
c            ax(i,j,kdist(i,j))=0.0
c            ax(i+1,j,kdist(i,j))=0.0
c            ay(i,j,kdist(i,j))=0.0
c          endif
c          else
c           way(i,j,1) = 0.5
c          endif
c
c        way(i,j,kdim1) = 0.0
c        ax(i,j,kdim1) = 0.0
c        ay(i,j,kdim1) = 0.0
c        way(i,j,1) = 0.5
        enddo
c        write(*,*)"kdist,frac sample ",kdist(i,1),frac(i,1)
      enddo
      enddo
c
      call imag(hdx)
c
c partial blocks
c
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
c            if(ax(i,j,k).lt.1.0)ax(i,j,k) = 0.0
c            if(ay(i,j,k).lt.1.0)ay(i,j,k) = 0.0
c            if(way(i,j,k).lt.1.0)way(i,j,k) = 0.0
c             way(i,j,k) = ay(i,j,k)
c             if(ay(i,j,k).eq.0.0) bay(i,j,k) = 0.0
c
c compute slope location
c            if(k.lt.(kdim-5))then
c               hx(k) =  0.5*yinter+(k-1.0)/mslope
c            else
c               hx(k) = 0.5*yinter+((kdim-5)-1.0)/mslope
c            endif
c
c         
c             if(ay(i,j,k).lt.1.0.and.ay(i,j,k).ne.0.) then
c               ay(i,j,k) = 0.0
c             endif
           if(ax2(i,k).gt.0.0.and.ax2(i,k).lt.1.0) then
             ax(i,j,k) = 0.0
             if(mslope.gt.0.0) then
              if(ax2(i,k).gt.ax2(i,k-1)) then
                way(i,j,k) = 1.-ax2(i,k)
                ay(i,j,k) = 1.0
              else if(ax2(i,k-1).eq.1.0) then
                way(i,j,k) = 1.0-ax2(i,k)
                ay(i,j,k) = 1.0
              endif
             else
              if(ax2(i,k).lt.ax2(i,k-1)) then
                way(i,j,k) = 1.-ax2(i,k-1)
                ay(i,j,k) = 1.0
                ax(i+1,j,k) = 1.0
c              else if(hdx(i,j,k+1).eq.1.0) then
              endif
              if(ax2(i,k+1).eq.1.0) then
                way(i,j,k+1) = 1.0-ax2(i,k)
                ay(i,j,k) = 1.0
              endif
             endif

           endif
          enddo
       enddo
      enddo
c
c do boundaries for kdist using real imag2d
c
c      do j=1,jdim
c        do i=1,idim
c          tmpp(i,j) = kdist(i,j)
c        enddo
c      enddo
c      call imag2d(tmpp)
c set kdist to kdim1 to avoid using limits
c      do k=0,kdim1     
c      do j=0,jdim1
c        if(iblkstart(my_pe+1).eq.1) then
c          if(ax(1,j,k).gt.0.0) ax(1,j,k) = 1.0
c          if(ay(1,j,k).gt.0.0) ay(1,j,k) = 1.0
c          if(way(1,j,k).gt.0.0) way(1,j,k) = 1.0
c        endif
c        do i=0,idim1
cc          kdist(i,j) = nint(tmpp(i,j))
c          kdist(i,j) = 1
c        enddo
c      enddo
c      enddo
cc
c
c      call imag2d(frac)
c
c compute volume for v component as an average
c
c      call imag(way)
c
      call imag(ax)
      call imag(ay)
c      call imag(bay)
c      call imag(bax)
c      call imag(vay)
      call imag(way)
      do i=0,idim1
        do j=0,jdim1
          way(i,j,kdim1) = 0.0
          way(i,j,0) = 0.0
          do k=1,kdim
           if(way(i,j,k).ne.0.0.and.ax(i+1,j,k).eq.0.0) then
             way(i,j,k) = 0.0
           endif
           bay(i,j,k) = ay(i,j,k)
           bax(i,j,k) = ax(i,j,k)
           bway(i,j,k) = ay(i,j,k)
c           way(i,j,k) = ay(i,j,k)
           vay(i,j,k) = ay(i,j,k)
c
c           if(ax(i,j,k).gt.0.0.and.ax(i,j,k).lt.1.0) then
c             ax(i,j,k) = 1.0
c             bax(i,j,k) = 0.0
c           endif
c           if(ay(i,j,k).gt.0.0.and.ay(i,j,k).lt.1.0) then
c             ay(i,j,k) = 1.0
c             bay(i,j,k) = 0.0
c           endif
c           if(way(i,j,k).gt.0.0.and.way(i,j,k).lt.1.0) then
c             way(i,j,k) = 1.0
c             bway(i,j,k) = 0.0
c           endif

          enddo
        enddo
      enddo
c      do i=0,idim
c        do j=0,jdim1
c           do k=0,kdim
c           if(mod(k,3).eq.0.and.bay(i,j,k).eq.0.0
c     *            .and.bay(i+1,j,k).ne.0.0) then
c            ax(i+1,j,k)= 0.0
c           endif
c          enddo
c        enddo
c      enddo
c      read(*,*)junk
c      call imag(bay)
c      call imag(bax)
c
      return
      end
