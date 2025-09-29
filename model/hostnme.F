      subroutine hostnme(name)
      character name*(*)
c
c   $Id: hostnme.F,v 1.2 1999/07/15 19:00:17 dwd Exp $
c
c returns the name of the machine that the code is running on
c
#if defined(sun) || defined(hp)
      call hostnm(name)
#endif
#if defined(cray)
      call gethost(name)
#endif
#if defined(Win32)
#if defined(NDEBUG)
      name = 'Windows NT host (Release)'
#else
      name = 'Windows NT host (DEBUG)'
#endif
#endif
#if defined(sgi)
c      character gethostname
c      name = gethostname()
      name = 'SGI Unix host'
#endif
      return
      end
