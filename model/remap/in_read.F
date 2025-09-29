      subroutine in_read(ier)
c
c   $Id: in_read.F,v 1.6 1998/11/04 00:25:33 dwd Exp $
c
c read the parameter file for 2-D model
c input cards independent of order,   dwd Apr 95
c
#include "param.inc"
#include "olemp.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "refdate.inc"
#include "communicate.inc"
c
      real dtemp,dsal
      integer tempsize,salsize
c
      common/iessize/dtemp,dsal,tempsize,salsize
c
      character str*80,value*80,name*15,icard*30
      logical flgs(3),mkstst
      integer unts(3),ier,i,nt,ns
      real xtemp,outstart(5),outstop(5),outdt(5),diagt,rstrt,slabt(5)
c 
c initialize strings and variables
c
      icvstr=' '
      icfstr=' '
      notimes=0
      slbtimes=0
      ntitle=0
      nrstrt=0
      ndiag=0
      zimethod=1
      zicriteria=0.01
      ttime=0.0
      diagt=0.0
      rstrt=0.0
c
      carddate  = .false.
      usedate   = .true.
      useicdate = .true.
      usetopref = .false.
      useray = .false.
c
      do i=1,5
         title(i)=' '
         nofld(i)=0
         noslb(i)=0
      enddo
      ier=0
 20   read (punit,30,end=60) icard
 30   format (a30)
      call upper(icard,30)
c 
      if (icard(:5).eq.'TITLE') then
c 
c title
c 
         ntitle=ntitle+1
         if (ntitle.gt.5) ntitle=5
         read (punit,40,err=50) title(ntitle)
 40      format (a80)
c
      else if (icard(:4).eq.'TIME') then
c     
c Time step and run time
c 
         read (punit,40) str
c
c delt (s), time (s)
c
         call parse(str,value)
         call mks(value,delt,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
         call parse(str,value)
         call mks(value,ttime,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
      else if (icard(:9).eq.'MODELFLAG') then
c 
c Model flags
c 
         read (punit,40) str
c
c advsal, filtord, filtskp, filtvert, doanal, rotate, useomegy
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(1)=int(xtemp)
         advsal=mdflgs(1).eq.1
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(2)=int(xtemp)
         filtord=mdflgs(2)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(3)=int(xtemp)
         filtskp=mdflgs(3)
c
c with Adams-Bashforth filtskp=50 is a reasonable number
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(4)=int(xtemp)
         filtvert=mdflgs(4).eq.1
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(5)=int(xtemp)
         doanal=mdflgs(5).eq.1
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(6)=int(xtemp)
         rotate=mdflgs(6).eq.1
c     
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         mdflgs(7)=int(xtemp)
         useomegy=mdflgs(7).eq.1
c
      else if (icard(:10).eq.'MODELVALUE') then
c 
c Model values
c 
         read (punit,40) str
c
c latitude, filtfac
c
         call parse(str,value)
         call mks(value,xlat,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c
c compute fz = 2 omega sin(latitude)
c         fy = 2 omega cos(latitude)
c
         if(xlat.gt.90. .or. xlat.lt.-90.)then
            tomegaz=0.0
            tomegay=0.0
         else
            tomegaz = 2.0*omega*sin(xlat*0.017453293)
            tomegay = 2.0*omega*cos(xlat*0.017453293)
         endif
c     
         call parse(str,value)
         call mks(value,filtfac,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000 
         mdvals(1)=filtfac
c     
      else if (icard(:5).eq.'FIELD') then
c 
c Output control
c 
         notimes=notimes+1
         if(notimes.gt.5)notimes=5
c
         read (punit,40) str
c
c start (s), stop (s), dt (s), nofld
c
         call parse(str,value)
         call mks(value,outstart(notimes),flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c     
         call parse(str,value)
         call mks(value,outstop(notimes),flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c     
         call parse(str,value)
         call mks(value,outdt(notimes),flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c     
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         nofld(notimes)=int(xtemp)
c     
         if (nofld(notimes).ne.0) then
            read(punit,40)str
c     
            do i=1,nofld(notimes)
               call parse(str,name)
               call varcde(name,ofield(notimes,i))
            enddo
         end if
c     
      else if (icard(:4).eq.'SLAB') then
c
c Slab output control
c
         slbtimes=slbtimes+1
         if(slbtimes.gt.5)slbtimes=5
c
         read (punit,40) str
c
c slabt (s), siloc, sjloc, skloc, type, noslb
c or
c slabt (s), sxloc (m), syloc (m), szloc (m), type, noslb
c
         call parse(str,value)
         call mks(value,slabt(slbtimes),flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
c get index from coordinate assuming cell-centered values
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0)) then
           if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0)) then
              goto 1000
           else
              sxloc(slbtimes) = xtemp
              siloc(slbtimes)= -1
           endif
         else
            siloc(slbtimes)=int(xtemp)
         endif
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0)) then
           if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0)) then
              goto 1000
           else
              syloc(slbtimes)= xtemp
              sjloc(slbtimes)= -1
           endif
         else
            sjloc(slbtimes)=int(xtemp)
         endif
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0)) then
           if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0)) then
              goto 1000
           else
              szloc(slbtimes) = xtemp
              skloc(slbtimes)= -1
           endif
         else
            skloc(slbtimes)=int(xtemp)
         endif
c
         call parse(str,value)
         call upper(value,80)
         if(value(:3).eq.'ALL')then
            slbtype(slbtimes)=SLAB_ALL
         else if(value(:4).eq.'VERT') then
            slbtype(slbtimes)=SLAB_VERT
         else
            slbtype(slbtimes)=SLAB_ALL
         endif
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         noslb(slbtimes)=int(xtemp)
c
         if (noslb(slbtimes).ne.0) then
            read(punit,40)str
c     
            do i=1,noslb(slbtimes)
               call parse(str,name)
               call varcde(name,oslab(slbtimes,i))
            enddo
         end if
c         
      else if (icard(:10).eq.'DIAGNOSTIC') then
c 
c Diagnositic control
c 
         read (punit,40) str
c
c diagt (s), diloc, djloc, dkloc
c
         call parse(str,value)
         call mks(value,diagt,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         diloc=int(xtemp)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         djloc=int(xtemp)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         dkloc=int(xtemp)      
c
      else if (icard(:7).eq.'RESTART') then
c
c Restart control
c
         read(punit,40)str
c
c rstrt (s), dump_end
c
         call parse(str,value)
         call mks(value,rstrt,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
         if(str.eq.' ') then
            xtemp = 0
         else
            call parse(str,value)
            call mks(value,xtemp,flgs,unts)
            if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         endif
         enddumpon=int(xtemp).eq.1
c
      else if (icard(:6).eq.'BCFLAG') then
c 
c B.C. flags
c 
         read (punit,40) bcfstr
c
      else if (icard(:7).eq.'BCVALUE') then 
c     
c B.C. values
c 
         read (punit,40) bcvstr
c
      else if (icard(:7).eq.'BCFORCE') then
c
c B.C. force
c
         read (punit,40) bcfrcstr
c
      else if (icard(:6).eq.'ICFLAG') then
c
c I.C. flags
c
         read (punit,40) icfstr
c     
      else if (icard(:7).eq.'ICVALUE') then
c
c I.C. values
c
         read (punit,40) icvstr
c
      else if (icard(:7).eq.'SGSFLAG') then
c
c SGS flags
c
         read (punit,40)sgsfstr
c
      else if (icard(:8).eq.'SGSVALUE') then
c
c SGS values
c
         read (punit,40)sgsvstr
c
c #if !defined(RADIATE)
      else if (icard(:8).eq.'RAYLEIGH') then
c
c do Rayleigh top?
c
         read (punit,40)str
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         useray=int(xtemp).eq.1
c
c taumax (s)
c
         call parse(str,value)
         call mks(value,taumax,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c #endif

      else if (icard(:6).eq.'TRACER') then 
c 
c tracers
c
         read (punit,*)ntrace
         if(ntrace.gt.trcdim)ntrace=trcdim
c
      else if (icard(:5).eq.'IES80') then
c
c iestbl (thtamin, thtamax, salmin, salmax)
c
         read (punit,40)str
c     
         call parse(str,value)
         call mks(value,thtamin,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c     
         call parse(str,value)
         call mks(value,thtamax,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,salmin,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,salmax,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c
c fix tempmax and salmax for actual table size
c     
         nt=min(tempsize,int((thtamax-thtamin)/dtemp)+1)
         ns=min(salsize,int((salmax-salmin)/dsal)+1)
         thtamax=thtamin+dtemp*(nt-1)
         salmax=salmin+dsal*(ns-1)
c
      else if (icard(:9).eq.'APPROXEOS') then
c
c approx equation of state (appden, linden, potden, centhta, censal)
c
         read (punit,40)str
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         appden=int(xtemp).eq.1
c     
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         linden=int(xtemp).eq.1
c     
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         potden=int(xtemp).eq.1
c     
         call parse(str,value)
         call mks(value,centhta,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c     
         call parse(str,value)
         call mks(value,censal,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
c
      else if (icard(:7).eq.'REFDATE') then
c
c reference date (YYYY-MM-DD hh:mm:ss)
c
         read (punit,40) str
         call parsedate(str,cardref)
         carddate = .true.
c
      else if (icard(:8).eq.'DATEFLAG') then
c
c Date flags
c
         read (punit,40) str
c
c usedate, useicdate
c
         call parse(str, value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         usedate = int(xtemp).eq.1
c
         call parse(str, value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         useicdate = int(xtemp).eq.1
c
      else if (icard(:8).eq.'REFFRAME') then
c
c Model reference frame, if 1 then subtract off top current velocities
c otherwise do nothing
c
         read (punit,40) str
c
c usetopref
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         usetopref=int(xtemp).eq.1
c
      else if (icard(:2).eq.'ZI') then
c
c Set the method to compute mixed layer depth (zi)
c
c
c presently two methods zimethod=0, abs(dt) < zicriteria and
c                                1, sigma-t < zicriteria
c
c default is zimethod=0, zicriteria=0.01
c
         read (punit,40) str
c
c zimethod, zicriteria (don't specify units!!!)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         zimethod=int(xtemp)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         zicriteria=xtemp
c
      else if (icard(:3).eq.'END') then
c 
c end of input, compute variables
c 
         ntime=int(ttime/delt+0.5)
c
         do i=1,notimes
            nostart(i)=int(outstart(i)/delt+0.5)
            nostop(i)=int(outstop(i)/delt+0.5)
            if(nofld(i).ne.0)then
               nodt(i)=int(outdt(i)/delt+0.5)
            else
               nodt(i)=0
            endif
         enddo
c
         if(diagt.ne.0.0)then
            ndiag=int(diagt/delt+0.5)
         else
            ndiag=0
         end if
c
         nrstrt=int(rstrt/delt+0.5)
c
         do i=1,slbtimes
            if(noslb(i).ne.0.0)then
               nslab(i)=int(slabt(i)/delt+0.5)
            else
               nslab(i)=0
            end if
         enddo
c         
         return
      else
         write(*,998)icard
 998     format('Unrecognized card type: ',a30)
         ier=1
      end if
      go to 20
 1000 write(*,999)icard,value(:30)
 999  format('Units error on card ',a30,' for entry: ',a)
 50   ier=1
 60   return
      end
