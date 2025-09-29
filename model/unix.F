      character*30 function modid()
c
c   $Id: unix.F,v 1.3 1999/07/15 19:00:17 dwd Exp $
c
c      unix version of modid
c 
#if defined(Win32)
      integer*2 ihr,imin,isec,i100th,iyr,imon,iday
      character str*30
      call gettim(ihr,imin,isec,i100th)
      call getdat(iyr,imon,iday)
      write(str,999)ihr,imin,isec,imon,iday,iyr
  999 format(i2.2,':',i2.2,':',i2.2,' ',i2,'/',i2.2,'/',i4)
      modid=str
#endif
#if defined(sun) || defined(hp) || defined(sgi)
      character dte*24
      call fdate(dte)
      modid=dte
#endif
#if defined(cray)
      character str*30,dte*8,clk*8
      call date(dte)
      call clock(clk)
      str=dte
      str(10:)=clk
      modid=str
#endif
      return
      end
