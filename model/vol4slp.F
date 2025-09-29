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
      real ax2(0:idim1,0:jdim1,0:kdim1)
      real hd2(0:idim1,0:kdim1)
      real dx2,dz2,frc1
      real yunder
c      real way(0:idim1,0:jdim1,0:kdim1)
c      real frac(0:idim1,0:jdim1)
c      integer kdist(0:idim1,0:jdim1)
      integer ind1,kmin,junk,olev
      integer ibeg,iend,i,j,k,flip
      integer mslope1,mslope2
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
        do j=0,jdim1
          ax2(i,j,k) = 1.0
          if(k.eq.0) ax2(i,j,k) = 0.0
        enddo
      enddo 
      enddo
c
c start with 45 degree slope to see if ax is correct (mslope1)
c
c undercut positive
c
c      mslope1 = 1.0
c      mslope2 = 4.0
cc      mslope2 = 6.0
c      mslope = 4.0
cc      mslope = 6.0
c      yinter = 0.0
c      yinter2 = 0.0
c      startk = 10.0

c
c all overcut is negative
c
      mslope1 = -6.
      mslope2 = -6.
      mslope = -6.
      yinter = kdim/(-mslope1)
      yinter2 = yinter
c
      do j=1,jdim1
      jj = jblkstart(my_pe+1)+j-1

      if(abs(mslope).gt.0.0) then
        dz2 = dz
        dx2 = dxx
        ibeg = iblkstart(my_pe+1)
        iend = ibeg+idim-1
        if(jbeg.lt.495.and.(jbeg+jdim).gt.495) then
          write(*,*),my_pe,ibeg,iend
        endif
        flip = 1
        ax2(i,j,kdim1) = 0.0
        do k=1,kdim
c
c terrain coordinate in absolute units
c
          if(k.le.(kdim))then
c
            if(mslope.gt.0.0) then
              if(k.lt.startk) then
                hx(k) = 0.0
              else
                h1 =  yinter+(k-0.5-startk)/mslope1
                h2 =  yinter2+(k-0.5-startk)/mslope2
                if(h1.lt.(h2+5.)) then
                   hx(k) = h1
                else
                   hx(k) = h2+5.
                endif
              endif
            else
             if(k.lt.(kdim-5))then
               hx(k) =  yinter+(k-0.5)/mslope
             else
               hx(k) = yinter+((kdim-5)-0.5)/mslope
             endif
            endif
          endif
        enddo
c
        do k=1,kdim
c
c overcut bottom fix miniumum of 15 vertical levels
c
          nvert = nint(10/abs(mslope))+1
          if(iope) write(*,*)"number of levels nvert ",nvert
          if(mslope.lt.0.0) then
            if(hx(k).gt.yinter-nvert) hx(k) = yinter-nvert+1
          endif
c sin shape
c
cc         hx(k) = kdim+5.0 + 0.5*kdim*sin((kdim+k)*1.57/kdim)
c         hx(k) = 205.-200*sin(k*1.57/kdim)
         
c
c
          if(hx(k).le.2.1) then
             hx(k) = 2.1
          endif
c          yunder = kdim/mslope2
c          if(hx(k).gt.yunder-2) hx(k) = yunder-2
c
c constant no slope, make sure mslope1 is 1 for undercut notch
c
cc          if(hx(k).ge.4.1) then
c          if(hx(k).ge.6.1) then
cc             hx(k) = 4.1
c             hx(k) = 6.1
c          endif
          hdx(i,j,k) = hx(k)
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
              ax2(ii,j,k) = 0.0
            enddo
              ax2(i,j,k) = 1.- (hx(k)-ind1)
c             jj = jblkstart(my_pe+1)+j-1
c            if(jj.eq.5)write(*,*)"i,k,ind,ibeg,iend,hx,ax2 ",
c     *       i,k,ind1,ibeg,iend,hx(k),ax2(i,j,k)

          else if(ind1.lt.ibeg) then
              do i=0,idim1
                ax2(i,j,k) = 1.0
c                hd2(i,k) = 0.0
              enddo
          else if(ind1.gt.iend) then
            do i=0,idim1
              ax2(i,j,k) = 0.0
c              hd2(i,k) = hx(k)-iend+(idim-i)
            enddo
          endif
        enddo
      endif
      enddo
c
c do vertical
c
      do i=1,idim
        do j=1,jdim
          kdist(i,j) = kdim
          way(i,j,kdim1) = 0.0
          ax(i,j,kdim1) = 0.0
          ay(i,j,kdim1) = 0.0
          do k=0,kdim
            if(ax2(i,j,k).lt.1.0.and.ax2(i,j,k).gt.0.0) then
              ay(i,j,k) = 1.0
              ax(i,j,k) = 0.0
              ax(i+1,j,k) = 1.0
              way(i,j,k) = 1.0
            else if(ax2(i,j,k).eq.1.0) then
              ay(i,j,k) = 1.0
              ay(i,j+1,k) = 1.0
              ax(i,j,k) = 1.0
              ax(i+1,j,k) = 1.0
              way(i,j,k) = 1.0
            else if(ax2(i,j,k).eq.0.0) then
              ay(i,j,k) = 0.0
              ay(i,j+1,k) = 0.0
              ax(i+1,j,k) = 0.0
              ax(i,j,k) = 0.0
              way(i,j,k) = 0.0
              way(i,j,k+1) = 0.0
            endif

              
        enddo
      enddo
      enddo
c
c reset way bc
c
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
            if(ax2(i,j,k).eq.0.0) then
                    way(i,j,k) =0.0
                    way(i,j,k+1) = 0.0
            endif
          enddo
        enddo
      enddo
c

c
      call imag(hdx)
c
c partial blocks
c
      do i=1,idim
        do j=1,jdim
          do k=1,kdim1
c
           if(ax2(i,j,k).gt.0.0.and.ax2(i,j,k).lt.1.0) then
             if(mslope.gt.0.0) then
                if(ax2(i,j,k).lt.ax2(i,j,k-1)) then
                  way(i,j,k) = ax2(i,j,k)
                else if(ax2(i,j,k-1).eq.1.0) then
                  way(i,j,k) = ax2(i,j,k)
                endif
                if(ax2(i,j,k+1).eq.0.0) then
                  way(i,j,k+1)= 0.0
                endif
             else
                if(ax2(i,j,k).gt.ax2(i,j,k-1)) then
                   way(i,j,k) = ax2(i,j,k-1)
                else
                   way(i,j,k) = ax2(i,j,k)
                endif
                if(ax2(i,j,k-1).eq.0.0) way(i,j,k) = 0.0
             endif
             bway(i,j,k) = way(i,j,k)
           endif
           if(way(i,j,k).ne.0.0.and.ax(i+1,j,k).eq.0.0) then
             way(i,j,k) = 0.0
           endif

          enddo
       enddo
      enddo
c
c
c
      call imag(ax)
      call imag(ay)
      call imag(way)
      call imag(bway)
      do i=0,idim1
        do j=0,jdim1
          way(i,j,kdim1) = 0.0
          way(i,j,0) = 0.0
          bax(i,j,kdim1) = 0.0
          bay(i,j,kdim1) = 0.0
          vay(i,j,kdim1) = 0.0
          do k=1,kdim
           vay(i,j,k) = ay(i,j,k)
           bay(i,j,k) = ay(i,j,k)
           bax(i,j,k) = ax(i,j,k)

           if(way(i,j,k).ne.0.0.and.ax(i+1,j,k).eq.0.0) then
             way(i,j,k) = 0.0
           endif
         enddo
        enddo
      enddo
c
      return
      end
