      subroutine deltime(sec)
      real sec
c   $Id: deltime.F,v 1.3 1999/07/15 19:00:17 dwd Exp $
c
#if defined(cray)
      real oldtime,time
      data oldtime/0.0/
      save oldtime
#endif
#if defined(sun) || defined(sgi)
      real time(2),dtime
#endif
#if defined(hp)
      real time(2),dtime
#endif
#if defined(sp2)
      
      real dtime_
      TYPE TB_TYPE
      SEQUENCE
        real (4) USRTIME
        real (4) SYSTIME
      END TYPE
      type (TB_TYPE) dtime_struct
c
#endif
#if defined(Win32)
      real oldtime,time
      save oldtime
      integer*2 ihr,imin,isec,i100th
      data oldtime/0.0/
#endif
c
c returns the number of seconds of cpu time since last call
c
#if defined(sun) || defined(hp) || defined(sgi)
      sec = dtime(time)
#endif
#if defined(sp2)
      sec = dtime_(dtime_struct)
#endif
#if defined(cray)
      call second(time)
      sec=time-oldtime
      oldtime=time
#endif
#if defined(Win32)
      call gettim(ihr,imin,isec,i100th)
      time = ihr*3600.0+imin*60.0+isec+i100th*0.01
      sec=time-oldtime
      if(sec.lt.-0.1)then
         write(*,*)'deltime: clock overflow reset by 86400'
         write(*,*)'deltime: oldtime = ',oldtime
         write(*,*)'deltime:    time = ',time
         sec=sec+86400.0
      endif
      oldtime=time
#endif
      return
      end
