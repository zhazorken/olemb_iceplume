#if defined(sun)
      subroutine sig_handler(sig)
c
c  $Id: sig_handler.F,v 1.3 1998/10/16 21:21:25 dwd Exp $
c
c
      integer sig
c
#include "param.inc"
#include "olemp.inc"
#include "signal.inc"
#include "moddef.inc"
c
      if(sig.eq.15) then
         sig_quit=.true.
         write(lunit,1)'TERM',timstp,'TERMINATING'
 1       format(/'Caught ',a,' signal at time step ',i6,' -- ',a)
      else if(sig.eq.24) then
         sig_restart=.true.
         sig_quick=.true.
         write(lunit,1)'XCPU',timstp,'WRITING RESTART AND TERMINATING'
      else if(sig.eq.30) then
         sig_restart=.true.
         sig_usr1=.true.
         write(lunit,1)'USR1',timstp,'WRITING RESTART'
      else if(sig.eq.31) then
         sig_restart=.true.
         sig_usr2=.true.
         write(lunit,1)'USR2',timstp,'WRITING RESTART AND TERMINATING'
      else if(sig.eq.8) then
	 sig_restart=.true.
	 sig_usr2=.true.
	 write(lunit,1)'FPE', timstp,'Writing restart and terminating'
      endif
c
      call flush(lunit)
c
      return
      end
#else
      subroutine sig_hander(sig)
      integer sig
c
c no sig_handler function
c
      end
#endif
