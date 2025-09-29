      subroutine filter(zz,fzz)
c
c   $Id: filter.F,v 1.4 1998/07/24 23:15:27 dwd Exp $
c
c --- zz = field to be filtered
c --- nord = order of the filter
c
c --- this subroutine filters a field of data. it works in combination
c --- with subroutine shpiro. filter picks the column or row of data,
c --- calls shpiro, and then adds a factor to the data (as calculated
c --- by shpiro.
c
c returns the filter coefficients of zz in fzz, should be added to zz!
c
#include "param.inc"
#include "olemp.inc"
c
      real zz(0:idim1,0:jdim1,0:kdim1),fzz(0:idim1,0:jdim1,0:kdim1)
      real horzx(0:idim1),vert(0:kdim1)
#if defined(USEMPI)
      real horzy(0:idim1,0:jdim1,0:kdim1)
#else
      real horzy(0:jdim1)
#endif
      real fac,fsign,facv
      integer nord,iodev,i,j,k
c     integer kord
c
      nord=filtord
c
      iodev=(nord+1)/2 -nord/2
      fac=-1.+2.*float(iodev)
      fac=fac/2.**(2*nord)
c
      fac=filtfac*fac
c
c filter horizontally in x direction
c
      do 90 k=1,kdim
         do 100 j=1,jdim
            do 200 i=1,idim
               horzx(i)=zz(i,j,k)
 200        continue
c
            call shpiro_x(horzx,nord,j,k)
c
            do 300 i=1,idim
               fzz(i,j,k)=fac*horzx(i)
 300        continue
 100     continue
 90   continue
c
c filter horizontally in y direction
c
c first make sure its not really 2-D
c
      if(jdim.eq.1) goto 401
c
c
#if defined(USEMPI)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            horzy(i,j,k) = zz(i,j,k)+fzz(i,j,k)
          enddo
        enddo
      enddo
      call shpiro_y_mpi(horzy,nord)
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            fzz(i,j,k)=fzz(i,j,k)+fac*horzy(i,j,k)
          enddo
        enddo
      enddo
      
#else
      do 400 k=1,kdim
         do 500 i=1,idim
            do 600 j=1,jdim
               horzy(j)=zz(i,j,k) + fzz(i,j,k)
 600        continue
c
            call shpiro_y(horzy,nord,i,k)
c     
            do 700 j=1,jdim
               fzz(i,j,k)=fzz(i,j,k) + fac*horzy(j)
 700        continue    
 500     continue  
 400  continue
#endif
c
 401  if(.not.filtvert) return
c
c filter vertically
c
      do 1400 j=1,jdim
         do 1500 i=1,idim
            do 1600 k=1,kdim
               vert(k)=zz(i,j,k) + fzz(i,j,k)
 1600       continue
c
c implement constant gradient b.c. 
            vert(kdim1) = 2.*vert(kdim)-vert(kdim-1)
c
c use the no taper vertical filter, assumes zero gradient for b.c.'s
c
            call shpiro_vnt(vert,nord,j,k)
c
c do not use the del2 (nord=1) filterd points
c ie, skip vert(2) and vert(kdim-1).  vert(1) and vert(kdim)
c are zero (unfiltered).
c
            fsign=-1.
            facv=fsign*filtfac/16.0
c
c the filter increases in order by 1 for each point into the domain
c comment out the taper
c
c           do 1700 kord=2,nord-1
c              fzz(i,kord+1,j)=fzz(i,kord+1,j)+facv*vert(kord+1)
c              fzz(i,kdim-kord,j)=fzz(i,kdim-kord,j)+
c    *              facv*vert(kdim-kord)
c
c              fsign=-fsign
c              facv=fsign*filtfac/2.**(2*(kord+1))
c1700       continue
c
c apply nord filter to the rest of the points
c
            do 1800 k=nord+1,kdim-nord
               fzz(i,j,k)=fzz(i,j,k)+fac*vert(k)
 1800       continue
 1500    continue
 1400 continue
      return
      end
c
      subroutine shpiro_x(g,nord,j,k)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the x direction
c
#include "param.inc"
#include "olemp.inc"
#include "finitv.inc"
c
      real g(0:idim1)
      integer nord,i,kord,j,k
c
      do 100 kord=1,nord
         g(0)=g(idim)
         g(idim1)=g(1)
c
         do 200 i=0,idim
c            g(i)=(g(i+1)-g(i))
            g(i)=(g(i+1)-g(i))*(ax(i+1,j,k)+ax(i,j,k))*0.5
c            if(k.gt.kdist(i,j))then
c              g(i) = 0.0
c            endif
 200     continue
c
         do 300 i=idim,1,-1
c            g(i)=g(i)-g(i-1)
            g(i)=(g(i)-g(i-1))*(ax(i-1,j,k)+ax(i,j,k))*0.5
c            if(k.gt.kdist(i,j))then
c              g(i) = 0.0
c            endif
 300     continue
 100  continue
c
      return
      end
      subroutine shpiro_y(g,nord,ij,k)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the y direction
c
#include "param.inc"
#include "olemp.inc"
#include "finitv.inc"
c
      real g(0:jdim1)
      integer nord,i,ij,k,kord
c
      do 100 kord=1,nord
         g(0)=g(jdim)
         g(jdim1)=g(1)
c
         do 200 i=0,jdim
c            g(i)=g(i+1)-g(i)
            g(i)=(g(i+1)-g(i))*(ay(ij,i+1,k)+ay(ij,i,k))*0.5
 200     continue
c
         do 300 i=jdim,1,-1
            g(i)=g(i)-g(i-1)
 300     continue
 100  continue
c
      return
      end
c
#if defined(USEMPI)
      subroutine shpiro_y_mpi(g,nord)
c
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the y direction
c
#include "param.inc"
#include "olemp.inc"
#include "finitv.inc"
c
      real g(0:idim1,0:jdim1,0:kdim1)
      integer nord,i,j,k,kord
c
      do 100 kord=1,nord
         call imag(g)
c
         do k=1,kdim
           do  j=0,jdim
             do i=1,idim
               g(i,j,k)=0.5*(ay(i,j,k)+ay(i,j+1,k))*
     *                       (g(i,j+1,k)-g(i,j,k))
             enddo
           enddo
         enddo
c
	 call imag(g)
c
         do k=1,kdim
           do j=jdim,1,-1
             do i=1,idim
               g(i,j,k)=g(i,j,k)-g(i,j-1,k)
             enddo
           enddo
         enddo
 100  continue
      call imag(g)
c
      return
      end
#endif
c
      subroutine shpiro_vnt(g,nord,ii,j)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the z direction assuming
c --- a zero gradient condition except on the first pass at the
c --- top where input b.c.'s are used
c
#include "param.inc"
#include "olemp.inc"
#include "finitv.inc"
c
      real g(0:kdim1)
      integer nord,i,kord,ii,j
c
         g(0)=g(1)
      do 100 kord=1,nord
c
         do 200 i=0,kdim
            g(i)=ax(ii+1,j,i)*ax(ii+1,j,i+1)*(g(i+1)-g(i))
c            if(i.gt.kdist(ii,j)) then
c              g(i) = 0.0
c            endif
 200     continue
c
         do 300 i=kdim,1,-1
c            g(i)=g(i)-g(i-1)
            g(i)=ax(ii+1,j,i)*ax(ii+1,j,i-1)*(g(i)-g(i-1))
c            if(i.gt.kdist(ii,j)) then
c              g(i) = 0.0
c            endif
 300     continue
         g(0)=g(1)
         g(kdim1)=g(kdim)
 100  continue
c
      return
      end
c
      subroutine shpiro_v(g,nord)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
c shpiro_v does not use boundary conditions in the vertical
c it increases the order of the filter by 1 each point into
c the domain.  Thus, g(1),g(kdim)= 0., g(2),g(kdim-1) is order 1
c g(3),g(kdim-2) is order 2, etc...
c
#include "param.inc"
#include "olemp.inc"
c
      real g(0:kdim1),gsave
      integer nord,i,kord
      integer mm,mm1,ms,ms1
c
      mm=kdim
      mm1=kdim-1
      ms=1
      ms1=2
c
      do 100 kord=1,nord
c
c save first point, first pass overwrites
c last point is OK
c
         gsave=g(ms)
         do 10 i=ms,mm1
            g(i)=g(i+1)-g(i)
 10      continue
         do 20 i=mm1,ms1,-1
            g(i)=g(i)-g(i-1)
 20      continue
c
c restore g(ms)
c
         g(ms)=gsave
c
c redefine start and stop points for next application
c of the del2 filter
c
         mm = mm-1
         mm1= mm-1
         ms = ms+1
         ms1= ms+1
 100  continue
      return
      end


