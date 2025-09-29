      subroutine parsedate(str, val)
c
c read a string of the format YYYY-MM-DD hh:mm:ss and
c convert to the GeoTime (EPIC) format
c
c   $Id: date.F,v 1.1.1.1 1998/07/20 22:03:47 eds Exp $
c
      character*(*) str
      integer val(2)
c
      integer year, mon, day, hour, min, sec
      integer jm,jy,ja,jul
c
      read(str,999)year,mon,day,hour,min,sec
 999  format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)
c
      if(year .lt. 0) year = year + 1
      if(mon .gt. 2) then
         jy = year
         jm = mon + 1
      else
         jy = year - 1
         jm = mon + 13
      endif
      jul = int(365.25*jy)+int(30.6001*jm)+day+1720995
      if(day + 31*(mon+12*year) .ge. 588829) then
         ja = int(0.01*jy)
         jul = jul + 2 - ja + int(0.25*ja)
      endif
      val(1) = jul
      val(2) = (hour*3600 + min*60)*1000 + sec*1000
      return
      end
c
      subroutine frmtdate(val, str)
c
c convert a GeoTime value to a string with the format
c YYYY-MM-DD hh:mm:ss for use in NetCDF files
c
      character*(*) str
      integer val(2)
c
      real*8 jalpha
      integer ja, jb, jc, jd, je
      integer year, mon, day, hour, min, sec
c
      if(val(1) .gt. 2299161) then
         jalpha = (dble(val(1) - 1867216) - 0.25d0)/36524.25d0
         ja = val(1) + 1 + int(jalpha) - int(0.25*jalpha)
      else
         ja = val(1)
      endif
      jb = ja + 1524
      jc = int(6680.0d0 + (dble(jb - 2439870) - 122.1d0)/365.25d0)
      jd = int(365*jc + (0.25*jc))
      je = int(dble(jb - jd)/30.6001d0)

      day = int(jb - jd) - int(30.6001d0*dble(je))
      mon = int(je - 1)
      if(mon .gt. 12) mon = mon - 12
      year = int(jc - 4715)
      if(mon .gt. 2) year = year - 1
      if(year .le. 0) year = year - 1
      ja = val(2)/1000
      hour = int(ja/3600)
      min = int((ja - hour*3600)/60)
      sec = int((val(2) - (hour*3600 + min*60)*1000)/1000)

      write(str,999)year,mon,day,hour,min,sec
 999  format(i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)
      return
      end
c
      integer function datediff(ref1, sec1, ref2, sec2)
c
c compute the difference in seconds between two times using
c the reference times in GeoTime format
c
      integer ref1(2), sec1, ref2(2), sec2
c
      integer ddiff, mdiff, sdiff
c
      ddiff = ref1(1) - ref2(1)
      mdiff = ref1(2) - ref2(2)
      sdiff = sec1 - sec2
c
c convert all to seconds and add up the differences
c
      datediff = ddiff*86400 + mdiff/1000 + sdiff
      return
      end
