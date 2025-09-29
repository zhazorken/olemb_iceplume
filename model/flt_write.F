      subroutine flt_write(n)
      integer n
c
c   $Id: flt_write.F,v 1.4 1998/07/27 20:12:34 dwd Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "float.inc"
#include "netcdf.inc"
c
      integer corner(4),edges(4),l,iret,tval
      integer i,j,k,nn
      real tmpx,tmpy,tmpz,xi,yi,zi
      real ti(fltdim),t1,t2,t3,t4,t5,t6,t7,t8
      real si(fltdim),s1,s2,s3,s4,s5,s6,s7,s8
      real tmpxyz(fltdim)
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
            tmpz = kdim-fltz(l,nn)/dz+0.5
c
            i = ifix(tmpx)
            j = ifix(tmpy)
            k = ifix(tmpz)
            xi = tmpx - i
            yi = tmpy - j
            zi = tmpz - k
c
            if(k.le.0)then
               t1=th(i,j,1)
               t2=th(i+1,j,1)
               t5=th(i,j+1,1)
               t6=th(i+1,j+1,1)
               ti(nn)=t1+(t2-t1)*xi+(t5-t1)*yi+(t6-t5+t1-t2)*xi*yi
            else if(k.ge.kdim)then
               t1=th(i,j,kdim)
               t2=th(i+1,j,kdim)
               t5=th(i,j+1,kdim)
               t6=th(i+1,j+1,kdim)
               ti(nn)=t1+(t2-t1)*xi+(t5-t1)*yi+(t6-t5+t1-t2)*xi*yi
            else
               t1=th(i,j,k)
               t2=th(i+1,j,k)
               t3=th(i+1,j,k+1)
               t4=th(i,j,k+1)
               t5=th(i,j+1,k)
               t6=th(i+1,j+1,k)
               t7=th(i+1,j+1,k+1)
               t8=th(i,j+1,k+1)
               ti(nn)=t1+(t4-t1)*zi+(t2-t1)*xi+(t1-t2+t3-t4)*zi*xi+
     *              (t5-t1)*yi+(t6-t5+t1-t2)*xi*yi+(t8-t5+t1-t4)*zi*yi+
     *              (t7-t8+t5-t6+t4-t3+t2-t1)*xi*yi*zi
            endif
c
            if(k.le.0)then
               s1=sal(i,j,1)
               s2=sal(i+1,j,1)
               s5=sal(i,j+1,1)
               s6=sal(i+1,j+1,1)
               si(nn)=s1+(s2-s1)*xi+(s5-s1)*yi+(s6-s5+s1-s2)*xi*yi
            else if(k.ge.kdim)then
               s1=sal(i,j,kdim)
               s2=sal(i+1,j,kdim)
               s5=sal(i,j+1,kdim)
               s6=sal(i+1,j+1,kdim)
               si(nn)=s1+(s2-s1)*xi+(s5-s1)*yi+(s6-s5+s1-s2)*xi*yi
            else
               s1=sal(i,j,k)
               s2=sal(i+1,j,k)
               s3=sal(i+1,j,k+1)
               s4=sal(i,j,k+1)
               s5=sal(i,j+1,k)
               s6=sal(i+1,j+1,k)
               s7=sal(i+1,j+1,k+1)
               s8=sal(i,j+1,k+1)
               si(nn)=s1+(s4-s1)*zi+(s2-s1)*xi+(s1-s2+s3-s4)*zi*xi+
     *              (s5-s1)*yi+(s6-s5+s1-s2)*xi*yi+(s8-s5+s1-s4)*zi*yi+
     *              (s7-s8+s5-s6+s4-s3+s2-s1)*xi*yi*zi
            endif
 30         continue
c
            call ncvpt(fltid, fltt_id(l), corner, edges,
     *           ti(1), iret)
            call ncvpt(fltid, flts_id(l), corner, edges,
     *           si(1), iret)
         endif
 20   continue
      call ncsnc(fltid, iret)
      return
      end
