#if defined(sun)
      integer function myhandler(sig,code,context)
c
c   $Id: myhandler.F,v 1.3 1998/07/24 23:15:34 dwd Exp $
c
      integer sig,code,context(5)
      real cputime,fctr
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemtim.inc"
#include "communicate.inc"
c
      integer taskmaster
c
      taskmaster = -10
c
c     write(*,*)'exception code ',code,' at ',context(4)
      myhandler=0
      call fld_write(timstp,1)
c
      cputime=veltim+prestim+advtim+updtetim+floattim+writetim
      fctr=1.0e6/(float((timstp-nstart+1))*float(idim*kdim*jdim))
      if(my_pe.eq.taskmaster) then
        write(lunit,63)veltim,veltim*100./cputime,
     *     prestim,prestim*100./cputime,advtim,advtim*100./cputime,
     *     updtetim,updtetim*100./cputime,floattim,
     *     floattim*100./cputime,writetim,writetim*100./cputime,
     *     cputime,cputime*fctr
 63   format(/'         Velocity calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '         Pressure calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '        Advection calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '           Update calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '            Float calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '               Writing Output ',f10.2,'  CPU seconds',
     *     f10.2,' percent'//
     *        '     Model execution required ',f10.2,'  CPU seconds'//
     *        '      per time step per point ',f11.3,
     *     ' CPU microseconds')
      endif
      call abort()
#if defined(USEMPI)
      call MPI_ABORT(MPI_COMM_WORLD,code)
#endif
      end
#else
      integer function myhandler(sig, code, context)
      integer sig,code,context(5)
c
c no myhandler function
c
      myhandler=0
      end
#endif
