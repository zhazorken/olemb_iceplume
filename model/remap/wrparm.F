      subroutine wrparm
c
c   $Id: wrparm.F,v 1.5 1998/11/04 00:25:35 dwd Exp $
c
c write the parameters for 2-D model
c 
      character code*50,str*50
      character slbx*12,slby*12,slbz*12
      integer sxlen, sylen, szlen
c
#include "param.inc"
#include "olemp.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "refdate.inc"
c
      integer i,j,len,lnblk
c 
c TITLE
c
      write (lunit,10)
 10   format (/'Titles:')
      do i=1,ntitle
         len=lnblk(title,80)
         write (lunit,20) title(i)(:len)
 20      format (1x,a)
      enddo
c
c TIME
c
      write(lunit,25) delt,ttime/3600.,ntime
 25   format(/'Time:'/
     *     5x,'1.        time step = ',  f12.1,' (sec)'/
     *     5x,'2.         end time = ',f12.3,' (hours)'//
     *     5x,'   total time steps = ',i12)
c
c REFDATE
c
      if(carddate) then
         call frmtdate(cardref, str)
         write(lunit,30)str, cardref(1), cardref(2)
 30      format(/'Reference Date from Input: ', a20, 2x, 2i10)
      else
         write(lunit,31)
 31      format(/'No Reference Date from Input')
      endif
c
c DATEFLAGS
c
      write(lunit,40)
 40   format(/'Date Flags:')
      if(usedate) then
         write(lunit,62) 1,'Use Date System'
         if(useicdate) then
            write(lunit,62) 2,'Use Initial Condition/Restart Date'
         else
            write(lunit,62) 2,
     *           'Do not use Initial Condition/Restart Date'
         endif
      else
         write(lunit,62) 1,'Do not use Date System'
      endif
c
c REFFRAME
c
      write(lunit,64)
 64   format(/'Reference Velocity:')
      if(usetopref)then
         write(lunit,62) 1, 'Use reference velocity from top'
      else
	 write(lunit,62) 0, 'Do not use reference velocity'
      endif
c
c MODELFLAG
c
      write (lunit,61)
 61   format(/'Model Flags:')
      if(advsal) then
         write(lunit,62) 1,'Advect Salinity'
 62      format(4x,i2,'. ',a)
      else
         write(lunit,62) 1,'Do not advect Salinity'
      endif
      write(lunit,63)filtord*2
 63   format(4x,' 2. Order of 2dx filter = ',i3)
      write(lunit,67)filtskp
 67   format(4x,' 3. Apply filter every ',i2,' time steps')
      if(filtvert)then
         write(lunit,62) 4,'Filter in the vertical'
      else
         write(lunit,62) 4,'Do not filter in the vertical'
      endif
      if(doanal)then
         write(lunit,62) 5,'Write analysis fields'
      else
         write(lunit,62) 5,'Do not write analysis fields'
      endif
      if(rotate)then
         write(lunit,62) 6,'Rotating Model Domain'
      else
         write(lunit,62) 6,'Non-rotating Model Domain'
      endif
      if(useomegy)then
         write(lunit,62) 7,'Use omega-y rotation component'
      else
         write(lunit,62) 7,'Do not use omega-y rotation component'
      endif
c
c MODELVALUE
c
      write (lunit,65)
 65   format(/'Model Values:')
      write (lunit,66)xlat,filtfac,tomegaz,tomegay
 66   format(4x,' 1.      latitude = ',  f12.3,' (degrees)'/
     *     4x,' 2. Filter factor = ',f12.5//
     *     4x,'      2 x Omega_z = ',1pe12.5,' (rad/sec)'/
     *     4x,'      2 x Omega_y = ',  e12.5,' (rad/sec)')
c
c FIELD
c
      write(lunit,68)
 68   format(/'Field Output:')
      do i=1,notimes
         if (nodt(i).eq.0) then
            write (lunit,70)i
 70         format (4x,i2,'. No Output Fields')
         else
            write (lunit,80)i,nostart(i)*delt/3600.,
     *           nostop(i)*delt/3600.,nodt(i)*delt/3600.
 80         format (4x,i2,'. Output from',f8.2,'(h) until',f8.2,
     *           '(h) every',f8.2,'(h)')
            do j=1,nofld(i)
               call varstr (ofield(i,j),code)
               write (lunit,90) code
 90            format (8x,a)
            enddo
         end if
      enddo
c
c SLAB
c
      write(lunit,239)
 239  format(/'Slab Output:')
      do i=1,slbtimes
         if (nslab(i).eq.0) then
            write (lunit,240)i
 240        format (4x,i2,'. No Slab Output')
         else
            if(slbtype(i).eq.SLAB_ALL)then
               str='x, y, and z'
               len=11
            else
               str='z'
               len=1
            endif
c
            if(siloc(i) .eq. -1) then
               write(slbx, 251) sxloc(i)
 251           format(f8.1,' (m)')
            else
               write(slbx, 252) siloc(i)
 252           format(i6)
            endif
            call squish(slbx,1,12)
            sxlen = lnblk(slbx,12)
c
            if(sjloc(i) .eq. -1) then
               write(slby, 251) syloc(i)
            else
               write(slby, 252) sjloc(i)
            endif
            call squish(slby,1,12)
            sylen = lnblk(slby,12)
c
            if(skloc(i) .eq. -1) then
               write(slbz, 251) szloc(i)
            else
               write(slbz, 252) skloc(i)
            endif
            call squish(slbz,1,12)
            szlen = lnblk(slbz,12)
c
            write (lunit,250)i,str(:len),slbx(1:sxlen),slby(1:sylen),
     *           slbz(1:szlen), nslab(i)*delt/60.
 250        format(4x,i2,'. Write ',a,' slabs at (',
     *           a,',',a,',',a,') every ',f8.2,' (min)')
            do j=1,noslb(i)
               call varstr(oslab(i,j),code)
               write(lunit,90)code
            enddo
         endif
      enddo
c
c DIAGNOSTIC
c
      write(lunit,109)
 109  format(/'Diagnostic Output:')
      if (ndiag.eq.0) then
          write (lunit,110)
 110      format (5x,'1. No Diagnostic Output')
      else
          write (lunit,120) diloc,djloc,dkloc,ndiag*delt/60.
 120      format (5x,'1. Write diagnostics at (',i3,',',i3,',',i3,
     *         ') every ',f8.2,' (min)')
      end if
c
c RESTART
c
      write(lunit,139)
 139  format(/'Restart Output:')
      if (nrstrt.eq.0)then
          write(lunit,140)
 140      format(5x,'1. No writing of restart file')
      else
          write(lunit,150)nrstrt,nrstrt*delt/3600.
 150      format(5x,'1. Write restart file every ',i5,' time steps (',
     1     f8.3,' hours)')
      endif
      if(enddumpon) then
         write(lunit,141)
 141     format(5x,'2. Force writing of restart file at end')
      else
         write(lunit,142)
 142     format(5x,'2. Don''t force writing of restart file at end')
      endif
c
c BCFLAG, BCVALUE, BCFORCE echoed in bcinit
c
c ICFLAG, ICVALUE echoed in icinit
c
c SGSFLAG, SGSVALUE echoed in sgsinit
c
c TRACER
c
      write (lunit,60) ntrace
 60   format (/'Tracer:'/
     *     5x,'1. number of tracers = ',i3)
c 
c IES80
c
      write(lunit,160)thtamin,thtamax,salmin,salmax
 160  format(/'IES-80 Table:'/
     *     5x,'1.    Theta min = ',f7.3/
     *     5x,'2.    Theta max = ',f7.3/
     *     5x,'3. Salinity min = ',f7.3/
     *     5x,'4. Salinity max = ',f7.3)
c
c APRROXEOS
c
      write(lunit,169)
 169  format(/'Approximate EOS:')
      if(appden)then
         write(lunit,170)centhta, censal
 170     format(5x,'1. Use approximate Equation of State'/
     *        5x,'2.    Center theta = ',f7.3,/
     *        5x,'3. Center salinity = ',f7.3)
         if(linden)then
            write(lunit,180)
 180        format(5x,'4. Use linear approximation')
         else
            write(lunit,190)
 190        format(5x,'4. Use quadradic approximation')
         endif
         if(potden)then
            write(lunit,200)
 200        format(5x,'5. Do not include pressure effect')
         else
            write(lunit,210)
 210        format(5x,'5. Include pressure effect')
         endif
      else
         write(lunit,220)
 220     format(5x,'1. Use full Equation of State')
      endif
c
c ZI
c
      write(lunit,300)
 300  format(/'Zi computation method:')
      if(zimethod.eq.0) then
         write(lunit,310)zicriteria
 310     format(5x,'abs(delta theta) < ',f8.4,' (degC)')
      else if(zimethod.eq.1) then
         write(lunit,320)zicriteria
 320     format(5x,'delta sigma-t < ',f8.4,' (kg m-3)')
      else
         write(lunit,330)
 330     format(5x,'Improper mixed layer depth method chosen'/
     *        5x,'Using default method')
         zimethod=0
         zicriteria=0.01
         write(lunit,310)zicriteria
      endif
c
      return
      end
