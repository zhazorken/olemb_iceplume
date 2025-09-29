      subroutine flt_write(n)
      integer n
c
c   $Id: flt_write.F,v 1.4 1998/07/27 20:12:34 dwd Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "float.inc"
#include "communicate.inc"
#include "netcdf.inc"
c
      integer corner(4),edges(4),l,iret,tval
      integer i,j,k,nn
      real tmpx,tmpy,tmpz,xi,yi,zi
      real ti(fltdim),t1,t2,t3,t4,t5,t6,t7,t8
      real si(fltdim),s1,s2,s3,s4,s5,s6,s7,s8
      real tmpxyz(fltdim)
      real thf(0:idim1,0:jtotal1),salf(0:idim1,0:jtotal1)
      real thff(0:idim1,jtotal),salff(0:idim1,jtotal)
c
c
c collect the surface th and sal into thf, salf
c
      write(*,*)"Collecting th,sal ",my_pe
      call collect_xyf(th,thff,kdim)
      call collect_xyf(sal,salff,kdim)
      if(iope)then
      do j=1,jtotal
        do i=0,idim1
          thf(i,j) = thff(i,j)
          salf(i,j) = salff(i,j)
        enddo
      enddo
      do i=0,idim1
         thf(i,0) = thf(i,jtotal)
         salf(i,0) = salf(i,jtotal)
         thf(i,jtotal1) = thf(i,1)
         salf(i,jtotal1) = salf(i,1)
      enddo

      write(*,*)"Ready to write sal,th floats"
c
      do 10 i=1,4
         corner(i)= 1
         edges(i) = 1
 10   continue
      flttimecnt=flttimecnt+1
      corner(4) = flttimecnt
c
c first update time stamp
c
      tval=n*delt
      call ncvpt(fltid, flttime_id, flttimecnt, 1, tval, iret)
c
      do 20 l=1,ndesc
         if(n.ge.fltstrt(l) .and. n.le.fltstop(l))then
            edges(1) = fltdescsize(l)
            do 21 i=1,fltdescsize(l)
                 tmpxyz(i) = fltx(l,i)
 21         continue
            call ncvpt(fltid, fltx_id(l), corner, edges,
     *           tmpxyz(1), iret)
            do 22 i=1,fltdescsize(l) 
                 tmpxyz(i) = flty(l,i) 
 22         continue
            call ncvpt(fltid, flty_id(l), corner, edges,
     *           tmpxyz(1), iret)
            do 23 i=1,fltdescsize(l) 
                 tmpxyz(i) = fltz(l,i) 
 23         continue
            call ncvpt(fltid, fltz_id(l), corner, edges,
     *           tmpxyz(1), iret)
c
c compute temp and salinity
c
            do 30 nn=1,fltdescsize(l)
            tmpx = fltx(l,nn)/dx+0.5
            tmpy = flty(l,nn)/dy+0.5
c
            i = ifix(tmpx)
            j = ifix(tmpy)
            xi = tmpx - i
            yi = tmpy - j
c
c
               t1=thf(i,j)
               t2=thf(i+1,j)
               t5=thf(i,j+1)
               t6=thf(i+1,j+1)
               ti(nn)=t1+(t2-t1)*xi+(t5-t1)*yi+(t6-t5+t1-t2)*xi*yi
c
               s1=salf(i,j)
               s2=salf(i+1,j)
               s5=salf(i,j+1)
               s6=salf(i+1,j+1)
               si(nn)=s1+(s2-s1)*xi+(s5-s1)*yi+(s6-s5+s1-s2)*xi*yi
 30         continue
c
            call ncvpt(fltid, fltt_id(l), corner, edges,
     *           ti(1), iret)
            call ncvpt(fltid, flts_id(l), corner, edges,
     *           si(1), iret)
         endif
 20   continue
      call ncsnc(fltid, iret)
c
c iope endif
c
      endif
      return
      end
