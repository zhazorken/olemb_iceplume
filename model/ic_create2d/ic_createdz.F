c---------------------------------------------------------------------------
c program ic_create
c 
c last modified: 7-27-92 by JR Bridge
c
c suroutines called: mkfield, ic_cdf
c functions called: timer, lnblk in ~golem/src/util/lnblk.F
c NetCDF calls are also made.
c---------------------------------------------------------------------------
      program ic_create
c
c     This program writes OLEM initial condition files in CDF format.
c     Currently,it assumes that the user only wants to specify initial
c     fields for temperature and salinity (this program also computes the
c     density based on these specifications). The user must also specify the 
c     spatial discretization (dx,dy,dz), and the grid size (imax,jmax,kmax).

c     The initial fields are set by linearly interpolating from a single 
c     vertical profile of temperature and salinity.  The user provides the 
c     name of this ascii file for which each record has depth, temperature and
c     salinity.  Therefore, this program assume that the user desires an 
c     initial field with no horizontal variations in temperature or salinity.

c     This program also assumes that the user wants the initial time to be 
c     set to zero.

c     In the following comments, (cdf) will signify a comment/routine that 
c     is directly used in a call to the netcdf function, IC will refer to 
c     'initial condition'.  Variable comments which involve 'dummy' refer to
c     those variables which have no significant impact on program performance
c     other than for the few surrounding lines.

c     The variable listing order is according the order in which the variables
c     first appear, used or declared, in this program.

c     i,j,k = int indices for x-, y-, and z-directions
c     icnt = int data record counter
c     ttle = char title (dummy)
c     holdi,holdj,holdk = char representations of imax,jmax,kmax
c     ver = char SCCS version number (cdf)
c     ofile = output file name (cdf)
c     ifile = input file name (ascii)
c     timer() = char function for getting system time and date
c     d,t,s = real depth, temperature and salinity data arrays
c     start = real model time start (0.0 for IC) (cdf)
c     ident = char date and time (cdf)
c     mver = char SCCS version number
c     com1,com2,com3 = char IC comments (cdf)                  (olemi.inc)
c     dx,dy,dz = real spatial discretization (cdf)
c     imax,jmax,kmax = int grid size (cdf)
c     nsize = int parameter for dimension of d,t,s 
c     th = real IC temperature array (cdf)
c     sal = real IC salinity array (cdf)

c     trcdim     max # of tracer files                         (param.inc)
c     nsize      max # of lines in profile data file           (param.inc)


#include "param.inc"
#include "ic.inc"
#include "olemp.inc"
#include "olemi.inc"
#include "olemtke.inc"
#include "netcdf.inc"

      integer psize
      parameter (psize=300)

c-----Declare variables
      real          d(psize),t(psize),s(psize),dtrc(psize),trace(psize)
      real          du(psize),uu(psize),vv(psize),hfinfact,hfoutfact
      real          lfoutfact,lfinfact,stoutfact,stinfact
      real          tpamp,tpdepth,z,ssal,sigthta
      logical       utest,sheartst,titletst,geotst,tstst,tprtb
      integer       i,j,k,icnt,itrc,n,icntu,stx,sty,strad
      integer	    hfx,hfy,hfrad,lfx,lfy,lfrad,m
      integer       tpseed
      character     ttle*80,timer*30,ifile*80,ofile*80,ufile*80,str*80
      character	    pattern*80,name*80,curline*80,tsfile*80
      character	    type*80
 
c
c setup error checking
c
#if defined(sun)
      integer ieeer,myhandler,ieee_handler
      external myhandler
#endif
#if defined(sun) && !defined(__APOGEE__)
c
c set up ieee handler
c
      ieeer = ieee_handler('set','common',myhandler)
#endif
c-----Initialize variables

      sheartst=.false.
      titletst=.false.
      geotst=.false.
      tstst=.false.
      hftst=.false.
      lftst=.false.
      strtst=.false.
      tprtb=.false.
      ufile='   '
      ntrace = 0
      start=0.0
      ident=timer()
      mver='$Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $'
      com2='initial temperature and salinity fields'
c
      write(*,999)idim,jdim,kdim
 999  format(' Maximum domain size = (',i3,'x',i3,'x',i3,')'/)
c-----Get input from user: input profile filename, output netCDF
c-------filename, dx,dy,dz,imax,jmax,kmax, number of tracers, and
c-------tracer profile filenames.
c **** read ic input file ****
      write(*,998)
 998  format(' Enter input configuration file name: ',$)
      read (*,'(a)') ifile
      m=index(ifile,'.')-1
      ofile=ifile(1:m)//'_ic.nc'
      write(*,*)' Making ',ofile

      open(1,file=ifile,status='old')
      
 200   read(1,300,end=100)str
 300   format(a)
      if(str(1:1).eq.'#')goto 200
      curline=str
      
      call newparse(str, pattern)
      
      if(pattern.eq.'title') then
	print *,'found title'
         titletst=.true.
         ttle = str
         
      else if(pattern.eq.'ts_profile') then
	print *,'found ts_profile'
         tstst=.true.
         call newparse(str, tsfile)
         
      else if(pattern.eq.'geometry') then
	print *,'found geometry'
         geotst=.true.
         read(str,*)dx,dy,dz,imax,jmax,kmax
         if(imax.lt.1 .or. imax.gt.idim .or.
     *        jmax.lt.1 .or. jmax.gt.jdim .or.
     *        kmax.lt.1 .or. kmax.gt.kdim) then
            print *,'imax,jmax,kmax must be <= ',idim,jdim,kdim
            call exit(1)
         endif
         imax1=imax+1
         jmax1=jmax+1
         
      else if(pattern.eq.'shear') then
	print *,'found shear'
         sheartst=.true.
         call newparse(str, ufile)
         
      else if(pattern.eq.'tracer') then
	print *,'found tracer'
         ntrace = ntrace+1
         call newparse(str, trfile(i))

      else if(pattern.eq.'t_perturb') then
         print *,'found temp perturbation'
          tprtb=.true.
          read(str,*)tpamp,tpdepth,tpseed

      else if(pattern.eq.'mask') then
	print *,'found mask'
         call newparse(str, name)
         if(name.eq.'hf')then
	print *,'found hf'
            hftst=.true.
            call newparse(str, type)
            if(type.eq.'rect') then
               read(str, *)hfx,hfy,hfinfact,hfoutfact
               
            else if(type.eq.'circle') then
               read(str, *)hfrad,hfinfact,hfoutfact
               
            else if(type.eq.'array') then
c     NOT IMPLEMENTED
            else
c     ERROR
            endif
         else if(name.eq.'lf') then
	print *,'found lf'
            lftst=.true.
            call newparse(str, type)
            if(type.eq.'rect') then
               read(str, *)lfx,lfy,lfinfact,lfoutfact
               
            else if(type.eq.'circle') then
               read(str, *)lfrad,lfinfact,lfoutfact
               
            else if(type.eq.'array') then
c     NOT IMPLEMENTED
            else
c     ERROR
            endif
         else if(name.eq.'str') then
	print *,'found str'
            strtst=.true.
            call newparse(str, type)
            if(type.eq.'rect') then
               read(str, *)stx,sty,stinfact,stoutfact
               
            else if(type.eq.'circle') then
               read(str, *)strad,stinfact,stoutfact
               
            else if(type.eq.'array') then
c     NOT IMPLEMENTED
            else
c     ERROR
            endif
         endif
      else
c     ERROR
      endif
      goto 200
      
 100  if((.not.tstst).or.(.not.geotst).or.(.not.titletst)) then
         print *,'TS profile, title, and geometry are required 
     &		in the ic file'
         goto 1000
      endif


c-----Make netcdf output title line, and comments 1 and 3
c-------Title line
      do 20 i=1,80
        if(ttle(i:i).eq.'_') ttle(i:i)=' '
20    continue
      title(1)=ttle
c-------Comments 1 and 3
      com1=tsfile
      write(com3,991) ntrace
991   format(i5,' tracer(s) present')

c-----Read comment 2 and temperature & salinity profile data
      open(1,file=tsfile,status='old')
      read(1,401)com2
401   format(a120)
      do 10 i=1,psize
        read(1,*,end=11)d(i),t(i),s(i)
10    continue
      i=psize+1
11    icnt=i-1
      close(1)
      utest = .FALSE.
      do n=1,80
         if((sheartst).and.(ufile(n:n).ne.' ')) then
            utest = .TRUE.
         endif
      enddo
      if(utest) then
         open(1,file=ufile,status='old')
         read(1,401)com3
         do i=1,psize
            read(1,*,end=12)du(i),uu(i),vv(i)
         enddo
         i=psize+1
12       icntu = i-1
         close(1)
      endif
c
c define dz 
c
      zz = 0.0
      do i=1,kdim
        if(i.gt.kdim/5) then
c          if(i.le.8) then
c            dzv(i) = 10.*exp(1.61*(i-1)/7.)
c          else
            dzv(i) = dz
c          endif
        else
          dzv(i) = dz+4.*dz*sin(1.571*(kdim/5-i+1)/float(kdim/5))
        endif
c        dzv(i) = dz
        zz = zz+dzv(i)
        write(*,*)"dzv, i",dzv(i),i,zz
      enddo
      dzv(0) = dzv(1)
      dzv(kdim1) = dzv(kdim)

 
c-----Make temperature and salinity initial condition fields
      call mkfield(th,idim,jdim,kdim,imax,jmax,kmax,d,t,dzv,icnt)
      if(tprtb)call addprtrb(th,tpamp,tpdepth,tpseed)

#if defined(SNGLSALT)
      call mkfield(sal,idim,jdim,kdim,imax,jmax,kmax,d,s,dzv,icnt)
#else
      call mkfielddp(sal,idim,jdim,kdim,imax,jmax,kmax,d,s,dz,icnt)
#endif
C-------Build and write hyperslabs of values into rhobm1.

c #if defined(RADTOP)
c
c for radtop, assume water is shallow, otherwise will need to
c input a depth offset
c
        z = dzv(kmax)*0.5
c       do 50 k=0,kmax
c          z=(kdim-k+0.5)*dz
c #else
      do 50 k=kmax,1,-1
c         z=(kmax-k+0.5)*dz
c #endif
         ssal=sal(1,k,1)
         rhobm1(k)=sigthta(ssal,th(1,k,1),z*0.1)
         thbar(k) = th(1,k,1)
         salbar(k) = sal(1,k,1)
         write(*,*)"t,s,rhobm1,z ", th(1,k,1),ssal,rhobm1(k),z
         z = z + (dzv(k)+dzv(k-1))*0.5
 50   continue
c      rhobm1(kmax+1) = rhobm1(kmax)
c      rhobm1(0) = rhobm1(1)


      if(utest) then
c         call mkfield(u,idim1,jdim1,kdim1,imax,jmax,kmax,du,uu,dz,icntu)
c         call mkfield(v,idim1,jdim1,kdim1,imax,jmax,kmax,du,vv,dz,icntu)
      endif
      write(*,*) ' Done making temperature and salinity fields.'
#if defined(TRACERS)
c-----For each tracer file get comment and data, and make tracer field
      do 30 itrc=1,ntrace 
c---------Get the tracer comment and data
        open (1,file=trfile(itrc),status='old')
        read(1,401) trccom(itrc)
        do 35 j=1,psize
          read(1,*,end=36) dtrc(j),trace(j) 
35      continue
        j=psize+1 
36      icnt=j-1
        close(1)
c---------Make the tracer field
        call mkfield(trc(0,0,0,itrc),idim1,jdim1,kdim1,imax,jmax,kmax,
     &  dtrc,trace,dz,icnt)
        write(*,37)itrc
37      format(' Done making tracer field #',i1,'.')
30    continue
#endif

c-----Subroutine to create masks if needed
      if(hftst)then
	call mkmask(hf_tmask,dx,dy,idim,jdim,hfrad,hfx,hfy,
     &		  hfinfact,hfoutfact)
      endif
      if(lftst)then
	call mkmask(lf_tmask,dx,dy,idim,jdim,lfrad,lfx,lfy,
     &		  lfinfact,lfoutfact)
      endif
      if(strtst)then
	call mkstmask(str_tmask,str_tbc,dx,dy,idim,jdim,strad,stx,sty,
     &		  stinfact,stoutfact)
      endif

c-----Subroutine to create the netcdf initial conditions output file.
      call ic_cdf(ofile,utest)
1000      end

c------------------------------------------------------------------------------
c subroutine ic_cdf
c
c called by: program ic_create in ~golem/src/bin/ic_create.F
c subroutines called: setvar in ~golem/src/util/varcde.F
c functions called: lnblk in ~golem/src/util/lnblk.F
c NetCDF calls are made.
c---------------------------------------------------------------------------
      subroutine ic_cdf(ifile,uvtest)

c     This routine writes a 3d OLEM IC cdf file used for input to OLEM.
c     The contents of the file are based upon the temperature (th), 
c     salinity (sal), and possibly tracer (trc) fields which were created,
c     and the grid information (dx, dy, dz, imax, jmax, and kmax) entered 
c     by the user in the calling program. The only two arguments sent to 
c     this routine is the IC cdf filename, ifile, and the tracer filenames.
c
c     input variable uvtest is a logical variable that tests true if
c     fields of u and v are input
c
c     cdfout = int IC output file identifier (cdf)
c     rcode = int error code (not equal to zero if error occurs) (cdf)
c     i,j,k = int x, y and z direction counter indices
c     imax,jmax,kmax = int grid size in x,y, and z directions
c
c     xu() = real x-direction spatial location for u computations (cdf)
c     xudim,xuid = int dimension and variable identifers for xu (cdf)
c
c     yu() = real y-direction spatial location for u computations (cdf)
c     yudim,yuid = int dimension and variable identifiers for yu (cdf)
c
c     zu() = real z-direction spatial location for u computations (cdf)
c     zudim,zuid = int dimension and variable identifiers for zu (cdf)
c
c     xw() = real x-direction spatial location for w computations (cdf)
c     xwdim,xwid = int dimension and variable identifiers for xw (cdf)
c
c     zw() = real z-direction spatial location for w computations (cdf)
c     zwdim,zwid = int dimension and variable identifiers for zw (cdf)
c
c     yv() = real y-direction spatial location for v computations (cdf)     
c     yvdim,yvid = int dimension and variable identifiers for yv (cdf)
c
c     xdim,xid = 
c     ydim,yid =
c     zdim,zid = 
c
c     varid(1:4) = int xwdim,yudim,zudim,tdim (cell-centered) identifier (cdf)
c     thid = int temperature field identifier (cdf)
c     salid = int salinity field identifier (cdf)
c     rhoid = int density field identifier (cdf)
c     rhobm1(0:kdim1) = from olemf.inc
c     strtf(1:4),countf(1:4) = int counters for writing to cdfout (cdf)
c     
c     ilen, iii = char counters
c
c
c   
c 

#include "param.inc"
#include "olemp.inc"
#include "olemi.inc"
#include "olemtke.inc"
#include "netcdf.inc"
#include "ic.inc"
      integer psize
      parameter (psize=300)

      character ifile*80,name*15,hold*2
      integer cdfout
      integer lnblk,k,iii,j,ilen,i,j1,k1,itrc
      integer rcode,varid(4),varidr(4),varidu(4),varidv(4),varidm(4)
      integer xudim(1),xwdim(1),ydim(1),zudim(1),zwdim(1)
      integer tdim(1),xdim(1),zdim(1)
      integer xuid,xwid,yid,zuid,zwid,xid,zid,salid,thid,tid,rhoid
      integer salbarid,thbarid,dzvid
      integer hfid,lfid,strmid,strbid
      integer uid,vid,iret
      integer yudim(1),yuid,yvdim(1),yvid,trcid(trcdim)
      real xu(idim),xw(idim),zu(kdim)
      real zw(kdim),yv(jdim),yu(jdim),dzvv(kdim)
#if !defined(SNGLSALT)
      real ssal(0:idim1)
#endif
      real vals(1)
      integer corner(4),count(4)
      integer strtf(4),countf(4),countr(4),countm(4)
      logical uvtest

#include "moddef.inc"

C-----Create cdf file (CDF)
c      cdfout=nccre(ifile,ncclob,rcode)
      iret = nf_create(ifile,IOR(NF_CLOBBER,NF_64BIT_OFFSET),
     *        cdfout)

C-----Define dimensions
c
      xudim(1)=ncddef(cdfout,'xu',imax,rcode)
      xuid =ncvdef(cdfout,'xu',NCFLOAT,1,xudim,rcode)
      call ncapt(cdfout, xuid, 'epic_code', NCLONG, 1, 501, rcode)
      call ncaptc(cdfout, xuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yudim(1)=ncddef(cdfout,'yu',jmax,rcode)
      yuid=ncvdef(cdfout,'yu',NCFLOAT,1,yudim,rcode)
      call ncapt(cdfout, yuid, 'epic_code', NCLONG, 1, 500, rcode)
      call ncaptc(cdfout, yuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c      
      zudim(1)=ncddef(cdfout,'zu',kmax,rcode)
      zuid  = ncvdef(cdfout, 'zu', NCFLOAT, 1, zudim, rcode)
      call ncapt(cdfout, zuid, 'epic_code', NCLONG, 1, 1, rcode)
      call ncaptc(cdfout, zuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xwdim(1)=ncddef(cdfout,'xw',imax,rcode)
      xwid  = ncvdef(cdfout, 'xw', NCFLOAT, 1, xwdim, rcode)
      call ncapt(cdfout, xwid, 'epic_code', NCLONG, 1, 501, rcode)
      call ncaptc(cdfout, xwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yvdim(1)=ncddef(cdfout,'yv',jmax,rcode)
      yvid  = ncvdef(cdfout, 'yv', NCFLOAT, 1, yvdim, rcode)
      call ncapt(cdfout, yvid, 'epic_code', NCLONG, 1, 500, rcode)
      call ncaptc(cdfout, yvid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yvid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yvid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zwdim(1)=ncddef(cdfout,'zw',kmax,rcode)
      zwid  = ncvdef(cdfout, 'zw', NCFLOAT, 1, zwdim, rcode)
      call ncapt(cdfout, zwid, 'epic_code', NCLONG, 1, 1, rcode)
      call ncaptc(cdfout, zwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xdim(1)=ncddef(cdfout,'x',1,rcode)
      xid   = ncvdef(cdfout, 'x', NCFLOAT, 1, xdim, rcode)
      call ncapt(cdfout, xid, 'epic_code', NCLONG, 1, 501, rcode)
      call ncaptc(cdfout, xid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      ydim(1)=ncddef(cdfout,'y',1,rcode)
      yid   = ncvdef(cdfout, 'y', NCFLOAT, 1, ydim, rcode)
      call ncapt(cdfout, yid, 'epic_code', NCLONG, 1, 500, rcode)
      call ncaptc(cdfout, yid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zdim(1)=ncddef(cdfout,'z',1,rcode)
      zid   = ncvdef(cdfout, 'z', NCFLOAT, 1, zdim, rcode)
      call ncapt(cdfout, zid, 'epic_code', NCLONG, 1, 1, rcode)
      call ncaptc(cdfout, zid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      tdim(1)=ncddef(cdfout,'time',NCUNLIM,rcode)
      tid=ncvdef(cdfout, 'time', NCLONG, 1, tdim, rcode)
      call ncapt(cdfout, tid, 'epic_code', NCLONG, 1, 625, rcode)
      call ncaptc(cdfout, tid, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, tid, 'units', NCCHAR, 7, 'seconds',rcode)
      call ncaptc(cdfout, tid, 'type' , NCCHAR, 1, ' ', rcode)
c

C-----Define cell-centered variables (CDF).
c-------These are the dimensions of the MDTEMP and MDSAL variables.
      varid(1)=xwdim(1)
      varid(2)=yudim(1)
      varid(3)=zudim(1)
      varid(4)=tdim(1)
c-------These are the dimensions of the mask variables.
      varidm(1)=xwdim(1)
      varidm(2)=yudim(1)
      varidm(3)=zdim(1)
      varidm(4)=tdim(1)
c-------These are the dimensions of the MDRHOBM1 variable.
      varidr(1)=xdim(1)
      varidr(2)=ydim(1)
      varidr(3)=zudim(1)
      varidr(4)=tdim(1)
c
      varidu(1)=xudim(1)
      varidu(2)=yudim(1)
      varidu(3)=zudim(1)
      varidu(4)=tdim(1)
c
      varidv(1)=xwdim(1)
      varidv(2)=yvdim(1)
      varidv(3)=zudim(1)
      varidv(4)=tdim(1)

C-----Set attributes.

c-------Obtain the var name corresponding to the const MDTEMP. (varnme)
c-------Obtain the cdf id number for this var & define it. (ncvdef)
c-------Set the attributes of this variable:  units, long name, 
c-------fortran format, and missing value. (setvar)
      call varnme(MDTEMP,name)
      ilen=lnblk(name,15)
      thid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
      call setvar(MDTEMP,cdfout,thid)

c-------Repeat the same steps for MDSAL.
      call varnme(MDSAL,name)
      ilen=lnblk(name,15)
      salid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
      call setvar(MDSAL,cdfout,salid)
c
      call varnme(MDDZV,name)
      ilen=lnblk(name,15)
      dzvid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidr,rcode)
      call setvar(MDDZV,cdfout,dzvid)
c
      if(uvtest) then
         call varnme(MDUVEL,name)
         ilen=lnblk(name,15)
         uid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidu,rcode)
         call setvar(MDUVEL,cdfout,uid)
c
         call varnme(MDVVEL,name)
         ilen=lnblk(name,15)
         vid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidv,rcode)
         call setvar(MDVVEL,cdfout,vid)
      endif

c-------Same steps for each tracer.
      do 75 i=1,ntrace
        call varnme(MDTRCR+i,name)
        ilen=lnblk(name,15)
        trcid(i)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
        call setvar(MDTRCR+i,cdfout,trcid(i))
75    continue

c-------Same steps for MDRHOBM1.
      call varnme(MDRHOBM1,name)
      ilen=lnblk(name,15)
      rhoid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidr,rcode)
      call setvar(MDRHOBM1,cdfout,rhoid)
c
c-------Same steps for MDSALBAR
      call varnme(MDSALBAR,name)
      ilen=lnblk(name,15)
      salbarid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidr,rcode)
      call setvar(MDSALBAR,cdfout,salbarid)
c
c-------Same steps for MDTHBAR
      call varnme(MDTHBAR,name)
      ilen=lnblk(name,15)
      thbarid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidr,rcode)
      call setvar(MDTHBAR,cdfout,thbarid)


      
c-------Same steps for MDHFTMASK.
      if(hftst) then
         call varnme(MDHFTMASK,name)
         ilen=lnblk(name,15)
         hfid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDHFTMASK,cdfout,hfid)
      endif
      
c-------Same steps for MDLFTMASK.
      if(lftst) then
         call varnme(MDLFTMASK,name)
         ilen=lnblk(name,15)
         lfid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDLFTMASK,cdfout,lfid)
      endif
      
c-------Same steps for MDSTRTMASK.
      if(strtst) then
         call varnme(MDSTRTMASK,name)
         ilen=lnblk(name,15)
         strmid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDSTRTMASK,cdfout,strmid)
         
         call varnme(MDSTRTBC,name)
         ilen=lnblk(name,15)
         strbid=ncvdef(cdfout,name(:ilen),NCLONG,4,varidm,rcode)
         call setvar(MDSTRTBC,cdfout,strbid)
      endif
      
C-------Setglobal attributes.
      iii=lnblk(ident,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'creation_date',NCCHAR,
     *     iii,ident,rcode)
      iii=lnblk(mver,80)
      call ncaptc(cdfout,NCGLOBAL,'version',NCCHAR,iii,mver,rcode)
      iii=lnblk(title(1),80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'title',NCCHAR,
     *     iii,title(1),rcode)
      call ncaptc(cdfout,NCGLOBAL,'DATA_TYPE',NCCHAR,5,'MODEL',rcode)
      call ncaptc(cdfout,NCGLOBAL,'COORD_SYSTEM',NCCHAR,5,'LOCAL',
     *     rcode)
      call ncaptc(cdfout,NCGLOBAL,'FILE_TYPE',NCCHAR,5,'START',rcode)
      call ncapt(cdfout,NCGLOBAL,'dx',NCFLOAT,1,dx,rcode)
      call ncapt(cdfout,NCGLOBAL,'dy',NCFLOAT,1,dy,rcode)
      call ncapt(cdfout,NCGLOBAL,'dz',NCFLOAT,1,dz,rcode)
      call ncapt(cdfout,NCGLOBAL,'imax',NCLONG,1,imax,rcode)
      call ncapt(cdfout,NCGLOBAL,'jmax',NCLONG,1,jmax,rcode)
      call ncapt(cdfout,NCGLOBAL,'kmax',NCLONG,1,kmax,rcode)
      call ncapt(cdfout,NCGLOBAL,'field_time',NCFLOAT,1,start,rcode)
      iii=lnblk(com1,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_1',NCCHAR,
     *     iii,com1,rcode)
      iii=lnblk(com2,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_2',NCCHAR,
     *     iii,com2,rcode)
      iii=lnblk(com3,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_3',NCCHAR,
     *     iii,com3,rcode)
      call ncapt(cdfout,NCGLOBAL,'ntrace',NCLONG,1,ntrace,rcode)
      do 10 itrc=1,ntrace
         write(hold,'(i1)')itrc
         iii=lnblk(trccom(itrc),120)
         call ncaptc(cdfout,NCGLOBAL,'trccom_'//hold,NCCHAR,
     *        iii,trccom(itrc),rcode)
         iii=lnblk(trfile(itrc),80)
         call ncaptc(cdfout,NCGLOBAL,'trcfile_'//hold,NCCHAR,
     *        iii,trfile(itrc),rcode)
 10   continue
      call ncendf(cdfout,rcode)  
      
C-------Buildthe axes.
      do 300 i=1,imax
         xu(i)=dx*float(i-1)
         xw(i)=dx*float(i-1)+dx*0.5
 300  continue
      do 310 i=1,kmax
         zu(i)=dz*float(kmax-i) + dz*0.5
         zw(i)=dz*float(kmax-i+1)
 310  continue
      do 320 i=1,jmax
         yv(i)=dy*float(i-1)
         yu(i)=yv(i)+dy*0.5
 320  continue
      
C-------Write a hyperslab of values into the axes variables using ncvpt.
      corner(1)=1
      count(1)=imax
      call ncvpt(cdfout, xuid, corner, count, xu, rcode)
      call ncvpt(cdfout, xwid, corner, count, xw, rcode)
      count(1)=jmax
      call ncvpt(cdfout, yuid, corner, count, yu, rcode)
      call ncvpt(cdfout, yvid, corner, count, yv, rcode)
      count(1)=kmax
      call ncvpt(cdfout, zuid, corner, count, zu, rcode)
      call ncvpt(cdfout, zwid, corner, count, zw, rcode)
      count(1)=1
      vals(1)=0.
      call ncvpt(cdfout, xid, corner, count, vals, rcode)
      call ncvpt(cdfout, yid, corner, count, vals, rcode)
      call ncvpt(cdfout, zid, corner, count, vals, rcode)
      call ncvpt(cdfout, tid, corner, count, vals, rcode)
      
C-------Tells ncvpt (for th, sal, rhobm1) where to start writing values.
      strtf(1)=1
      strtf(2)=1
      strtf(3)=1
      strtf(4)=1
C-------Tells ncvpt (for th, sal) how many values to write.
      countf(1)=imax
      countf(2)=1
      countf(3)=1
      countf(4)=1
C-------Tells ncvpt (for hf, lf, str) how many values to write
      countm(1)=imax
      countm(2)=1
      countm(3)=1
      countm(4)=1
C-------Tells ncvpt (for rhobm1) how many values to write.
      countr(1)=1
      countr(2)=1
      countr(3)=kmax
      countr(4)=1
      
c-------write out rhobm1,salbar,thbar
      call ncvpt(cdfout,rhoid,strtf,countr,rhobm1(1),rcode)
      call ncvpt(cdfout,salbarid,strtf,countr,salbar(1),rcode)
      call ncvpt(cdfout,thbarid,strtf,countr,thbar(1),rcode)
c---
      do k=1,kdim
        dzvv(k) = dzv(k)
      enddo
      call ncvpt(cdfout,dzvid,strtf,countr,dzvv,rcode)

      
C-------Write hyperslabs of values into the temp and salin variables.
      do k=1,kmax
         strtf(3)=k
         k1=k
         do j=1,jmax
            strtf(2)=j
            j1=j
            call ncvpt(cdfout,thid,strtf,countf,th(1,k1,j1),rcode)
#if defined(SNGLSALT)
            call ncvpt(cdfout,salid,strtf,countf,sal(1,k1,j1),rcode)
#else
            do i=1,imax
               ssal(i)=sal(i,k1,j1)
            enddo
            call ncvpt(cdfout,salid,strtf,countf,ssal,rcode)
#endif
#if defined(VELOCITIES)
            if(uvtest) then
               call ncvpt(cdfout,uid,strtf,countf,u(1,k1,j1),rcode)
               call ncvpt(cdfout,vid,strtf,countf,v(1,k1,j1),rcode)
            endif
#endif
         enddo
      enddo

C-------For each mask, write hyperslabs of values into its variable.
      strtf(3)=1
      strtf(4)=1
      if(hftst) then
         do j=1,jmax
            strtf(2)=j
            call ncvpt(cdfout,hfid,strtf,countm,hf_tmask(1,j),rcode)
         enddo
      endif
      
      if(lftst) then
         do j=1,jmax
            strtf(2)=j
            call ncvpt(cdfout,lfid,strtf,countm,lf_tmask(1,j),rcode)
         enddo
      endif
      
      if(strtst) then
         do j=1,jmax
            strtf(2)=j
            call ncvpt(cdfout,strmid,strtf,countm,str_tmask(1,j),rcode)
            call ncvpt(cdfout,strbid,strtf,countm,str_tbc(1,j),rcode)
         enddo
      endif
      
#if defined(TRACERS)
C-------For each tracer, write hyperslabs of values into its variable.
      do itrc=1,ntrace
         do 65 j=1,jmax
            strtf(2)=j
            j1=j
            do 60 k=1,kmax
               strtf(3)=k
               k1=k
               call ncvpt(cdfout,trcid(itrc),strtf,countf,
     &              trc(1,k1,j1,itrc),rcode)
 60         continue
 65      continue
      enddo
#endif
      
      call ncsnc(cdfout,rcode)
      write(*,*) ' Done making cdf file.'
      return
      end

c---------------------------------------------------------------------------
c subroutine mkfielddp  double precision version
c
c called by: program ic_create in ~golem/src/bin/ic_create.F
c---------------------------------------------------------------------------
      subroutine mkfielddp(fld,idm,jdm,kdm,im,jm,km,dd,f,dz,icnt)

#include "param.inc"
      integer psize
      parameter (psize=300)

      integer idm,jdm,im,jm,icnt,i,k,j,kk,km,kdm
      real*8 fld(0:idm,0:kdm,0:jdm)
      real dd(icnt),f(icnt),dz
      logical flip
      real*8 tmp,ddz,slope,z,d(psize)


      do 10 i=1,icnt
         d(i)=dd(i)
 10   continue

C-----Bubble sort data
20    flip=.false.
      do 30 i=2,icnt
         if(d(i).lt.d(i-1))then
            flip=.true.
            tmp=d(i)
            d(i-1)=d(i)
            d(i)=tmp
            tmp=f(i)
            f(i-1)=f(i)
            f(i)=tmp
         endif
30    continue
      if(flip)goto 20

C-----Interpolate data to grid points
      if(dz*0.5.lt.d(1) .or. (km-0.5)*dz.gt.d(icnt))then
         write(*,*)' Field not made: Profile file does not cover ',
     &             'depths.'
         stop
      endif

#if defined(RADTOP)
      do 40 k=1,km
         z=(k-0.5)*dz
#else
      do 40 k=km,1,-1
         z=(km-k+0.5)*dz
#endif
         do 45 kk=1,icnt-1
            if(z.ge.d(kk) .and. z.lt.d(kk+1))goto 46
45       continue

C--------Interpolate temperature and remove background
46       ddz=z-d(kk)
         slope=(f(kk)-f(kk+1))/(d(kk)-d(kk+1))
         tmp=f(kk)+slope*ddz
         do 50 i=1,im
           do 60 j=1,jm
             fld(i,k,j)=tmp
60         continue
50       continue
40    continue
      return
      end

c-----------------------------------------------------------------------
      subroutine mkmask(fld,dx,dy,im,jm,irad,mx,my,infact,outfact)

#include "param.inc"
#include "cyclic.inc"
      integer psize
      parameter (psize=300)

      integer im,jm,i,j,irad,mx,my
      real fld(0:idim1,0:jdim1),infact,outfact,tmpx, tmpy, cir, r
      real dx, dy
      integer stx, sty, endx, endy, cx, cy

      
      if(irad.eq.0)then
         tmpx=mx/dx
         tmpy=my/dy
         if(tmpx.gt.im)then
            stx = 0
            endx = im
         else
            stx = im - tmpx
            stx = stx/2
            endx = im-stx
         endif 
         if(tmpy.gt.jm)then
            sty = 0
            endy = jm
         else
            sty = jm - tmpy
            sty = sty/2
            endy = jm-sty
         endif 
         do i=0,im
            do j=0,jm
               if((i.lt.stx).or.(j.lt.sty).or.(i.gt.endx)
     &              .or.(j.gt.endy))then
                  fld(i,j)=outfact
               else
                  fld(i,j)=infact
               endif
            enddo
         enddo
         return
      else
         cx = im/2
         cy = jm/2
	 r = irad/dx
	 cir=(r*r)-(cx*cx)-(cy*cy)
         do i=0,im
            do j=0,jm
	       if((r.ge.cx).and.(r.ge.cy)) then
		  fld(i,j)=infact
	       else if(((i*i)+(j*j)-(2*i*cx)-(2*j*cy)).lt.cir) then
                  fld(i,j)=infact
               else
		  fld(i,j)=outfact
               endif
            enddo
         enddo              
         return
      endif
      end
      
c---------------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine mkstmask(fld,stfld,dx,dy,im,jm,irad,mx,my,
     &			infact,outfact)
      
#include "param.inc"
#include "cyclic.inc"
            integer psize
      parameter (psize=300)

      integer im,jm,i,j,stfld(0:idim1,0:jdim1),irad,mx,my
      real fld(0:idim1,0:jdim1),infact,outfact,cir,r,tmpx,tmpy
      real dx,dy
      integer stx, sty, endx, endy, cx, cy
      
      
      if(irad.eq.0)then
         tmpx=mx/dx
         tmpy=my/dy
         if(tmpx.gt.im)then
            stx = 0
            endx = im
         else
            stx = im - tmpx
            stx = stx/2
            endx = im-stx
         endif 
         if(tmpy.gt.jm)then
            sty = 0
            endy = jm
         else
            sty = jm - tmpy
            sty = sty/2
            endy = jm-sty
         endif 
         do i=0,im
            do j=0,jm
               if((i.lt.stx).or.(j.lt.sty).or.(i.gt.endx)
     &              .or.(j.gt.endy))then
                  fld(i,j)=outfact
		  stfld(i,j) = 1
               else
                  fld(i,j)=infact
		  stfld(i,j) = 1
               endif
            enddo
         enddo
      else
         cx = im/2
         cy = jm/2
	 r = irad/dx
	 cir=(r*r)-(cx*cx)-(cy*cy)
         do i=0,im
            do j=0,jm
	       if((r.ge.cx).and.(r.ge.cy)) then
		  fld(i,j)=infact
		  stfld(i,j) = 1
	       else if(((i*i)+(j*j)-(2*i*cx)-(2*j*cy)).lt.cir) then
                  fld(i,j)=infact
		  stfld(i,j) = 1
               else
		  fld(i,j)=outfact
		  stfld(i,j) = 1
               endif
            enddo
         enddo              
      endif
      
      if((infact.eq.0).or.(outfact.eq.0)) then
         do i=0,im
            do j=0,jm
               if(fld(i,j).eq.0) then
                  stfld(i,j) = 0
               endif
            enddo
         enddo
      endif
      
      return
      end
      
c---------------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine mkfield(fld,idm,jdm,kdm,im,jm,km,dd,f,dzv,icnt)

#include "param.inc"
      integer psize
      parameter (psize=300)

      integer idm,jdm,im,jm,icnt,i,k,j,kk,km,kdm
      real fld(idm,kdm,jdm),dd(icnt),f(icnt),dzv(kdm)
      logical flip
      real dz,tmp,ddz,slope,z,d(psize)


      do 10 i=1,icnt
         d(i)=dd(i)
 10   continue

C-----Bubble sort data
20    flip=.false.
      do 30 i=2,icnt
         if(d(i).lt.d(i-1))then
            flip=.true.
            tmp=d(i)
            d(i-1)=d(i)
            d(i)=tmp
            tmp=f(i)
            f(i-1)=f(i)
            f(i)=tmp
         endif
30    continue
      if(flip)goto 20

C-----Interpolate data to grid points
      if(dz*0.5.lt.d(1) .or. (km-0.5)*dz.gt.d(icnt))then
         write(*,*)' Field not made: Profile file does not cover ',
     &             'depths.'
         stop
      endif
c
c if using RADTOP then the profiles should be referenced from
c the bottom, not the top:
c
c    z   t
c    0   5
c   20   6
c
c  where z is now distance from the bottom and not depth.
c
c
#if defined(RADTOP)
      do 40 k=1,km
         z=(k-0.5)*dz
#else
      z = 0.5*dzv(km)
      do 40 k=km,1,-1
c         z=(km-k+0.5)*dz
#endif
         do 45 kk=1,icnt-1
            if(z.ge.d(kk) .and. z.lt.d(kk+1))goto 46
45       continue

C--------Interpolate temperature and remove background
46       ddz=z-d(kk)
         slope=(f(kk)-f(kk+1))/(d(kk)-d(kk+1))
         tmp=f(kk)+slope*ddz
         do 50 i=1,im
           do 60 j=1,jm
             fld(i,k,j)=tmp
60         continue
50       continue
         z = z+0.5*(dzv(k-1)+dzv(k))
40    continue
      return
      end

c---------------------------------------------------------------------------
c function timer    Assigns the date character string to the variable timer.
c---------------------------------------------------------------------------
      character*30 function timer()
#if defined(Win32)
      use dfport
#endif
      character dte*24
      call fdate(dte)
      timer=dte
      return
      end


      subroutine newparse(str,value)
      character str*(*),value*(*)
c
c   @(#)newparse.F     Version 2.1     7/13/92
c
c get next value from string  value is terminated by a comma,:, or space
c
      character tmp*120
      integer icma,icolin,ispace,ilen,i
c
c remove leading spaces
c
      ilen=len(str)
      do 10 i=1,ilen
	  if(str(i:i).eq.'#')goto 12
          if(str(i:i).ne.' ')goto 11
10    continue
12    value=' '
      return
11    tmp=str(i:)
      str=tmp
c
      icma=index(str,',')
      icolin=index(str,':')
      ispace=index(str,' ')
      if(icma.eq.0) then
	if(icolin.ne.0) then
	  if(ispace.ne.0) then
		if(icolin.gt.ispace) then
		  icma = ispace
		else
		  icma = icolin
		endif	
	  else
	    icma = icolin
	  endif
	else
	  if(ispace.ne.0) then
		icma = ispace
	  endif
	endif
      else
      if((icma.gt.icolin).and.(icolin.ne.0))icma=icolin
      if((icma.gt.ispace).and.(ispace.ne.0))icma=ispace
	endif
      if(icma.eq.0)then
          value=str
          str=' '
      else
          value=str(:icma-1)
          tmp=str(icma+1:)
          str=tmp
      endif
      return
      end

      subroutine addprtrb(fld,amp,depth,seed)
c
#include "param.inc"
#include "olemp.inc"
#include "ic.inc"
c
      integer psize
      parameter (psize=300)

      integer seed
      real fld(0:idim1,0:kdim1,0:jdim1),amp,depth
      real ranfld(0:idim1,0:kdim1,0:jdim1),ran2
      real filtfld(0:idim1,0:kdim1,0:jdim1)
      real*8 filtave(0:kdim1)
      integer i,k,j,kstart
c
c create random field
c
      do 10 j=1,jmax
         do 20 k=1,kmax
            do 30 i=1,imax
               ranfld(i,k,j)=(1.-2.*ran2(seed))*amp
 30         continue
 20      continue
 10   continue
c
c filter field
c
      filtord=2
      filtfac=0.5
      filtvert=.true.
      call filter(ranfld,filtfld)
c
c compute filtered average perturbation
c
      do 70 k=1,kmax
         filtave(k)=0.0
         do 80 j=1,jmax
            do 90 i=1,imax
               filtave(k)=filtave(k)+filtfld(i,k,j)
 90         continue
 80      continue
         filtave(k)=filtave(k)/(jmax*imax)
 70   continue
c
c combine fields
c
      kstart=kmax-int(depth/dz)+1
      kstart=max0(1,kstart)
      kstart=min0(kmax,kstart)
      do 40 j=1,jmax
         do 50 k=kstart,kmax
            do 60 i=1,imax
               fld(i,k,j)=fld(i,k,j)+filtfld(i,k,j)-filtave(k)
 60         continue
 50      continue
 40   continue
      return
      end

      subroutine filter(zz,fzz)
c
c --- zz = field to be filtered
c --- nord = order of the filter
c
c --- this subroutine filters a field of data. it works in combination
c --- with subroutine shpiro. filter picks the column or row of data,
c --- calls shpiro, and then adds a factor to the data (as calculated
c --- by shpiro.
c
c returns the filter coefficients of zz in fzz, should be added to zz!
c
#include "param.inc"
#include "olemp.inc"
#include "ic.inc"
c
      integer psize
      parameter (psize=300)

      real zz(0:idim1,0:kdim1,0:jdim1),fzz(0:idim1,0:kdim1,0:jdim1)
      real horzx(0:idim1),horzy(0:jdim1),vert(0:kdim1)
      real fac,fsign,facv
      integer nord,iodev,i,j,k,kord
c
      nord=filtord
c
      iodev=(nord+1)/2 -nord/2
      fac=-1.+2.*float(iodev)
      fac=fac/2.**(2*nord)
c
      fac=filtfac*fac
c
c filter horizontally in x direction
c
      do 90 j=1,jmax
         do 100 k=1,kmax
            do 200 i=1,imax
               horzx(i)=zz(i,k,j)
 200        continue
c
            call shpiro_x(horzx,nord)
c
            do 300 i=1,imax
               fzz(i,k,j)=fac*horzx(i)
 300        continue
 100     continue
 90   continue
c
c filter horizontally in y direction
c
c first make sure its not really 2-D
c
      if(jmax.eq.1) goto 401
c
c
      do 400 i=1,imax
         do 500 k=1,kmax
            do 600 j=1,jmax
               horzy(j)=zz(i,k,j) + fzz(i,k,j)
 600        continue
c
            call shpiro_y(horzy,nord)
c     
            do 700 j=1,jmax
               fzz(i,k,j)=fzz(i,k,j) + fac*horzy(j)
 700        continue    
 500     continue  
 400  continue
c
 401  if(.not.filtvert) return
c
c filter vertically
c
      do 1400 i=1,imax
         do 1500 j=1,jmax
            do 1600 k=1,kmax
               vert(k)=zz(i,k,j) + fzz(i,k,j)
 1600       continue
c
            call shpiro_v(vert,nord)
c
c do not use the del2 (nord=1) filterd points
c ie, skip vert(2) and vert(kmax-1).  vert(1) and vert(kmax)
c are zero (unfiltered).
c
            fsign=-1.
            facv=fsign*filtfac/16.0
c
c the filter increases in order by 1 for each point into the domain
c
            do 1700 kord=2,nord-1
               fzz(i,kord+1,j)=fzz(i,kord+1,j)+facv*vert(kord+1)
               fzz(i,kmax-kord,j)=fzz(i,kmax-kord,j)+
     *              facv*vert(kmax-kord)
c
               fsign=-fsign
               facv=fsign*filtfac/2.**(2*(kord+1))
 1700       continue
c
c apply nord filter to the rest of the points
c
            do 1800 k=nord+1,kmax-nord
               fzz(i,k,j)=fzz(i,k,j)+fac*vert(k)
 1800       continue
 1500    continue
 1400 continue
      return
      end
c
      subroutine shpiro_x(g,nord)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the x direction
c
#include "param.inc"
#include "olemp.inc"
#include "ic.inc"
c
      integer psize
      parameter (psize=300)

      real g(0:idim1)
      integer nord,i,kord
c
      do 100 kord=1,nord
         g(0)=g(imax)
         g(imax1)=g(1)
c
         do 200 i=0,imax
            g(i)=g(i+1)-g(i)
 200     continue
c
         do 300 i=imax,1,-1
            g(i)=g(i)-g(i-1)
 300     continue
 100  continue
c
      return
      end
      subroutine shpiro_y(g,nord)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter in the y direction
c
#include "param.inc"
#include "olemp.inc"
#include "ic.inc"
c
      integer psize
      parameter (psize=300)

      real g(0:jdim1)
      integer nord,i,kord
c
      do 100 kord=1,nord
         g(0)=g(jmax)
         g(jmax1)=g(1)
c
         do 200 i=0,jmax
            g(i)=g(i+1)-g(i)
 200     continue
c
         do 300 i=jmax,1,-1
            g(i)=g(i)-g(i-1)
 300     continue
 100  continue
c
      return
      end
c
      subroutine shpiro_v(g,nord)
c
c --- g = portion of a larger array to be filtered
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
c shpiro_v does not use boundary conditions in the vertical
c it increases the order of the filter by 1 each point into
c the domain.  Thus, g(1),g(kmax)= 0., g(2),g(kmax-1) is order 1
c g(3),g(kmax-2) is order 2, etc...
c
#include "param.inc"
#include "olemp.inc"
#include "ic.inc"
c
      integer psize
      parameter (psize=300)

      real g(0:kdim1),gsave
      integer nord,i,kord
      integer mm,mm1,ms,ms1
c
      mm=kmax
      mm1=kmax-1
      ms=1
      ms1=2
c
      do 100 kord=1,nord
c
c save first point, first pass overwrites
c last point is OK
c
         gsave=g(ms)
         do 10 i=ms,mm1
            g(i)=g(i+1)-g(i)
 10      continue
         do 20 i=mm1,ms1,-1
            g(i)=g(i)-g(i-1)
 20      continue
c
c restore g(ms)
c
         g(ms)=gsave
c
c redefine start and stop points for next application
c of the del2 filter
c
         mm = mm-1
         mm1= mm-1
         ms = ms+1
         ms1= ms+1
 100  continue
      return
      end
#if defined(sun)
      integer function myhandler(sig,code,context)
c
c   @(#)myhandler.F     Version 3.1     11/30/93
c
      integer sig,code,context(5)
      real cputime,fctr
c
#include "param.inc"
#include "olemp.inc"
#include "olemtim.inc"
c
      integer psize
      parameter (psize=300)

c     write(*,*)'exception code ',code,' at ',context(4)
      myhandler=0
c     call wrtflds(timstp,1)
c
c     fctr=1.e6/float(timstp*imax*kmax*jmax)
c     cputime=veltim+prestim+advtim+updtetim
c     write(*,*)'Model execution ',cputime,' seconds'
c     write(*,*)'Velocity calc   ',veltim,' seconds'
c     write(*,*)'Pressure calc   ',prestim,' seconds'
c     write(*,*)'Advection calc  ',advtim,' seconds'
c     write(*,*)'Update calc     ',updtetim,' seconds'
c     write(*,*)cputime*fctr,' microseconds per step-point'
      write(*,*)'ic_create - error'
c
      call abort()
      end
#endif

      integer function lnblk (str,ilen)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
      character str*(*)
      integer ilen,i
c 
      do 10 i=ilen,1,-1
          if (str(i:i).ne.' ') go to 20
 10   continue
      i=0
 20   lnblk=i
      return
      end

c-----------------------------------------------------------------------
      subroutine varcde (str,icode)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
      character*15 str
      integer icode,i,k
 
c     ivar  =  variable name (input)
c     icode =  variable code number (output)

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
#include "varinit.inc"
       integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(str.eq.mdvnme(i))goto 11
         if(str(1:4).eq.'trcr')goto 12
 10   continue
      icode=-1
      return
 11   icode=mdvcde(i)
      return
 12   read(str(5:5),'(I1)')k
      icode=MDTRCR+k
      return
      end

c-----------------------------------------------------------------------
      subroutine varnme(icode,str)
      character*15 str
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if(icode.ge.MDTRCR)then
         write(str,999)icode-MDTRCR
 999     format('trcr',i1)
      else
         str=' '
      endif
      return
 11   str=mdvnme(i)
      return
      end

c-----------------------------------------------------------------------
c  Assigns 4 attributes to the specified variable:  units, long_name,
c    FORTRAN_format, and missing_value.
c  Called by winit3, a subroutine of ~golem/src/bin/winit3.F,
c    and other programs...
c  Subroutines called: varunt, varstr, varfmt, which are in this file.
c  Functions called: lnblk from ~golem/src/include/lnblk.F
c  NetCDF subroutines are also called.
c-----------------------------------------------------------------------
      subroutine setvar(fld,cdfid,varid)

      integer fld,cdfid,varid
      integer ilen,lnblk,rcode
      character str*50,unt*20,fmt*10

#include "netcdf.inc"

      call varunt(fld,unt)
      ilen=max0(lnblk(unt,20),1)
      call ncaptc(cdfid,varid,'units',NCCHAR,ilen,unt,rcode)

      call varstr(fld,str)
      ilen=max0(lnblk(str,50),1)
      call ncaptc(cdfid,varid,'long_name',NCCHAR,ilen,str,rcode)

      call varfmt(fld,fmt)
      ilen=max0(lnblk(fmt,10),1)
      call ncaptc(cdfid,varid,'FORTRAN_format',NCCHAR,ilen,fmt,rcode)

      call ncapt(cdfid,varid,'missing_value',NCFLOAT,1,1.e35,rcode)

      return
      end

c-----------------------------------------------------------------------
      subroutine varstr(icode,str)
      character*50 str
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if(icode.ge.MDTRCR)then
         write(str,999)icode-MDTRCR
 999     format('Tracer ',i1)
      else
         str=' '
      endif
      return
 11   str=mdvstr(i)
      return
      end

c-----------------------------------------------------------------------
      subroutine varunt(icode,unt)
      character*20 unt
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      unt=' '
      return
 11   unt=mdvunt(i)
      return
      end

c-----------------------------------------------------------------------
c  modified for tracers by jrb on 7-20-92
c-----------------------------------------------------------------------
      subroutine varfmt(icode,fmt)
      character*10 fmt
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if (icode.ge.MDTRCR) then
         fmt='t6.2'
      else
         fmt='e10.2'
      endif
      return
 11   fmt=mdvfmt(i)
      return
      end

      real function ran2 (idum)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
c        returns a uniform random deviate between 0.0 and 1.0.  Set
c        IDUM to any negative value to initialize or reinitialize
c        the sequence.
c
      save
      integer m,ia,ic
      real rm 
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      integer ir(97),idum,iff,j,mod,iy
      data iff/0/
      if (idum.lt.0.or.iff.eq.0) then
          iff=1
          idum=mod(ic-idum,m)
          do 10 j=1,97
              idum=mod(ia*idum+ic,m)
              ir(j)=idum
 10       continue
          idum=mod(ia*idum+ic,m)
          iy=idum
      end if
      j=1+(97*iy)/m
      if (j.gt.97.or.j.lt.1)then
         write(*,*)'Problem with random number'
         return
      endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
