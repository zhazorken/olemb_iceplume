      subroutine bc_read(ifile)
c
c   $Id: bc_read.F,v 1.5 1998/08/04 20:58:10 dwd Exp $
c
#include "moddef.inc"
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "olemi.inc" 
#include "olemtke.inc"
#include "netcdf.inc"
#include "cyclic.inc"
#include "olemcdf.inc"
#include "refdate.inc"
#include "communicate.inc"
c 
      character timnam*30
      integer tstrt,tcount
      integer fstrt(4),fcount(4)
      character ifile*80,name*15
      character tunits*80,str*20
      integer sindx,i
      integer cdf,timout
      integer lnblk,ilen,n
      integer rcode
      integer tdim
      data fstrt/1,1,1,1/
      data fcount/1,1,1,1/

c
C-----open cdf file (CDF)
c
      cdf=ncopn(ifile,ncnowrit,rcode)
      if(iope) then
        if(rcode.ne.0.0) then
           write(lunit,*)' Boundary condition file non-existant'
     *       ,' use default conditions'
           return
        else
           write(lunit,*)' Boundary conditions cdf file opened'
        endif
      endif

c-----get global attributes
 
      call ncagtc(cdf,NCGLOBAL,'creation_date',bctime,30,rcode)
      call ncagtc(cdf,NCGLOBAL,'version',bcver,80,rcode)
      call ncagtc(cdf,NCGLOBAL,'bccom_1',bccom1,120,rcode)
      call ncagtc(cdf,NCGLOBAL,'bccom_2',bccom2,120,rcode)
      call ncagtc(cdf,NCGLOBAL,'bccom_3',bccom3,120,rcode)

      if(iope) then
        write(lunit,802)bccom2
802     format(' B.C. File title: ',a80)
      endif

c
C-----define dimensions

      tdim=ncdid(cdf,'time',rcode)
      call ncdinq(cdf,tdim,timnam,timsiz,rcode)
      bctid = ncvid(cdf,'time',rcode)
      tcount = timsiz
      tstrt = 1
      if(iope) then
        write(lunit,*)timsiz,' values in b.c. file'
      endif
c
c read units and try to find reference time
c
      bcdate = .false.
c
      call ncagtc(cdf,bctid,'units', tunits, 80, rcode)
      sindx = index(tunits,'since')
      if(sindx .gt. 0) then
         do i=sindx+5,len(tunits)
            if(tunits(i:i).ne.' ') then
               call parsedate(tunits(i:), bcref)
               bcdate = .true.
               goto 100
            endif
         enddo
      endif
 100  if(bcdate.and.iope) then
         call frmtdate(bcref, str)
         write(lunit,*) 'B.C. reference date = ', str
      endif
c
c read in the time axis values into tt
c
      call ncvgt(cdf,bctid,tstrt,tcount,tt,rcode)
c
c------Turn off netCDF error checking. Otherwise the routine
c      will bomb when a missing variable is encountered
c
      call ncpopt(0)



      call varnme(MDHFT,name)
      ilen=lnblk(name,15)
      hftid=ncvid(cdf,name(:ilen),rcode)
      do 10 n=1,timsiz
         heat(n) = -200.0
10    continue
      if(rcode.eq.0)then
	heaton = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in heat flux '
        endif
        call ncvgt(cdf,hftid,fstrt,fcount,heat,rcode)
      else
	heaton = .false.
      endif

      call varnme(MDUSTT,name)
      ilen=lnblk(name,15)
      ustrtid=ncvid(cdf,name(:ilen),rcode)
      do 11 n=1,timsiz
         ust(n) =  0.0
11    continue
      if(rcode.eq.0)then
	uston = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in u component stress'
        endif
        call ncvgt(cdf,ustrtid,fstrt,fcount,ust,rcode)
      else
	uston = .false.
      endif
        
      call varnme(MDVSTT,name)
      ilen=lnblk(name,15)
      vstrtid=ncvid(cdf,name(:ilen),rcode)
      do 12 n=1,timsiz
         vst(n) =  0.0
12    continue
      if(rcode.eq.0)then 
	vston = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in v component stress'
        endif
        call ncvgt(cdf,vstrtid,fstrt,fcount,vst,rcode) 
      else
	vston = .false.
      endif

      call varnme(MDWAVEL,name)
      ilen=lnblk(name,15)
      wavelid=ncvid(cdf,name(:ilen),rcode)
      do 13 n=1,timsiz
          wavel(n) =  0.0
13    continue
      if(rcode.eq.0)then 
        wavelon = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in surface wavelength'
        endif
        call ncvgt(cdf,wavelid,fstrt,fcount,wavel,rcode) 
      else
	wavelon = .false.
      endif
 
      call varnme(MDWAVEH,name)
      ilen=lnblk(name,15)
      wavehid=ncvid(cdf,name(:ilen),rcode)
      do 14 n=1,timsiz
          waveh(n) =  0.0
14    continue
      if(rcode.eq.0)then
        wavehon = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in surface wave height'
        endif
        call ncvgt(cdf,wavehid,fstrt,fcount,waveh,rcode)
      else
	wavehon = .false.
      endif

      call varnme(MDANGLE,name)
      ilen=lnblk(name,15)
      angleid=ncvid(cdf,name(:ilen),rcode)
      do 15 n=1,timsiz  
          angle(n) =  0.0
15    continue
      if(rcode.eq.0)then
        angleon = .true.
        fcount(4) = timsiz
        if(iope)then
          write(lunit,*)'Read in surface wave directional angle'
        endif
        call ncvgt(cdf,angleid,fstrt,fcount,angle,rcode)
      else
        angleon = .false.
      endif
c
c read in the surface short wave heat flux
c
      call varnme(MDSWHFT,name)
      ilen=lnblk(name,15)
      swhftid=ncvid(cdf,name(:ilen),rcode)
      do 16 n=1,timsiz
         swheat(n) = 0.0
16    continue
      if(rcode.eq.0)then
        swheaton = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in short wave heat flux '
        endif
        call ncvgt(cdf,swhftid,fstrt,fcount,swheat,rcode)
        if(iope) then
          write(lunit,*)'Initialize depth dependent irradiance'
        endif
        call initrad

      else
        swheaton = .false.
      endif
c
c read in the surface rain flux
c
      call varnme(MDRAIN,name)
      ilen=lnblk(name,15)
      rainid=ncvid(cdf,name(:ilen),rcode)
      do 17 n=1,timsiz
         rain(n) = 0.0
17    continue
      if(rcode.eq.0)then
        rainon = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in rain flux '
        endif
        call ncvgt(cdf,rainid,fstrt,fcount,rain,rcode)
      else
        rainon = .false.
      endif
c
c read in the surface latent heat flux, used to calculate salinity
c changes
c
      call varnme(MDLHFT,name)
      ilen=lnblk(name,15)
      lhftid=ncvid(cdf,name(:ilen),rcode)
      do 18 n=1,timsiz
         lheat(n) = 0.0
18    continue
      if(rcode.eq.0)then
        lheaton = .true.
        fcount(4) = timsiz
        if(iope) then
          write(lunit,*)'Read in latent heat flux '
        endif
        call ncvgt(cdf,lhftid,fstrt,fcount,lheat,rcode)
      else
        lheaton = .false.
      endif

c------Turn netCDF error checking back on.
      call ncpopt(NCVERBOS+NCFATAL)

c      timout = min(timsiz,10)
      timout = timsiz
      if(iope) then
        write(lunit,*)'Top boundary forcing'
        write(lunit,800)
        write(lunit,801)(tt(n)/3600.,heat(n),lheat(n),swheat(n),
     *      ust(n),vst(n),
     *      wavel(n),waveh(n),angle(n),
     *      rain(n),n=1,timout)
800     format(3x,'time',4x,'heat',2x,'latent',1x,'sw heat',1x,
     *      'ustress',1x,
     *      'vstress',2x,'wave l',
     *      2x,'wave h',2x,'wave a',4x,'rain')
801     format(1x,f6.1,3f8.1,2f8.3,2f8.2,2f8.1)
      endif
c
c all done 
c
      call ncclos(cdf,rcode)
      if(nprocs.gt.1)
     * write(*,*)' done reading ', ifile(1:60),' processor ',my_pe

      return
      end
