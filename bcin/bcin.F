      program bcin 
c
#include "param.inc"
#include "bc.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
c
      integer i,j,icnt,ival(1)
      character timer*30
      character*80 ofile
      character str*160, name*15
      character refdate*80
      logical datemode

      real t(nsize)
      real dat_array(20,nsize)
      integer time,strtt(4),countt(4), vars(20)
      integer rcode, ofield(20), nfld
      integer srt(1),cnt(1)
      logical first
c
      data countt/1,1,1,1/
      data strtt/1,1,1,1/
      data first/.true./

      bver="$Id: bcin.F,v 4.5 1998/11/09 18:25:22 eds Exp $"
      bident=timer()

      write(*,995)
995   format(' Enter input b.c. file name: ')
      read(*,998)ifile
      write(*,996)
996   format(' Enter output b.c. file name: ')
      read(*,998)ofile

      open(1,file=ifile,status='old')
      read(1,998)com2
      read(1,998)refdate
      read(1,'(a)')str
998   format(a80)
      com3=str
      datemode = refdate(1:4) .ne. '    '
      do 13 i=1,20
	call parse(str, name)
	call varcde(name, ofield(i))
	nfld=i-1
	if((ofield(i).eq.-1).or.(ofield(i).eq.999))goto 100
 13   continue
 100  do 10 i=1,nsize
          read(1,*,end=11)t(i),(dat_array(j,i),j=1,nfld)
10    continue
      i=nsize+1
11    icnt=i-1
      close(1)

      call initbc(ofile,ofield,nfld,vars,datemode, refdate)

      do 12 i=1,icnt
         time = int(t(i)*3600)
         ival(1) = time
         if(first)then
            first=.false.
            tindext = 1
            strtt(4) = 1
            srt(1) = 1
            cnt(1) = 1
            call ncvpt(bcid, bctid, srt, cnt, ival,rcode)
         else
            tindext = tindext + 1
            strtt(4)=tindext
            srt(1) = tindext
            cnt(1) = 1
            call ncvpt(bcid, bctid, srt, cnt, ival, rcode)
         endif
c
c  
c
	do 14 j=1,nfld
         call ncvpt(bcid,vars(j),strtt,countt,dat_array(j,i),rcode)
14	continue
c
c sync the netCDF file
c
         call ncsnc(bcid, rcode)
12    continue

      end

      subroutine initbc(froot,ofield,nfld,vars,datemode,refdate)
c
c
#include "param.inc"
#include "bc.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
c 
      character froot*80,name*15
      integer cdfout,kmax
      integer i,lnblk,iii,j,ilen
      integer rcode,varidp(4),vars(20)
      integer ydim(1),tdim(1),xdim(1),zdim(1)
      integer yid,xid,zid,ofield(20),nfld
      integer srt(1),cnt(1)
      real xval(1)
      logical datemode
      character refdate*80,tunits*80
      integer tlen
c
#include "moddef.inc"
c
      kmax=kdim
c
C-----create cdf file (CDF)
c
      j=lnblk(froot,80)
      bcid=nccre(froot(1:j)//' ',ncclob,rcode)
c
C-----define dimensions
c
      cdfout=bcid
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
      bctid   = ncvdef(cdfout, 'time', NCLONG, 1, tdim, rcode)
      if(.not.datemode)then
         call ncapt(cdfout, bctid, 'epic_code', NCLONG, 
     *        1, 625, rcode)
      endif
      call ncaptc(cdfout, bctid, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      if(datemode) then
         tunits = 'seconds since '// refdate
         tlen = lnblk(tunits,80)
      else
         tunits = 'seconds'
         tlen = 7
      endif
      call ncaptc(cdfout, bctid, 'units', NCCHAR, 
     *     tlen, tunits,rcode)
      call ncaptc(cdfout, bctid, 'type' , NCCHAR, 1, ' ', rcode)
c
c
C-----define non-coordinate variables
c
      varidp(1)=xdim(1)
      varidp(2)=ydim(1)
      varidp(3)=zdim(1)
      varidp(4)=tdim(1)
c
C---- set attributes

      do 15 i=1,nfld
	call varnme(ofield(i),name)
	ilen=lnblk(name,15)
        write(*,'("Processing: ",a15)')name(:ilen)
	vars(i)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
	call setvar(ofield(i),cdfout,vars(i))
15    continue
c
c global attributes for the boundary condition file
c
      iii=lnblk(bident,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'creation_date',NCCHAR,
     *     iii,bident,rcode)
      iii=lnblk(bver,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'version',NCCHAR,
     *     iii,bver,rcode)
      iii=lnblk(ifile,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_1',NCCHAR,
     *            iii,ifile,rcode)
      iii=lnblk(com2,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_2',NCCHAR,
     *            iii,com2,rcode)
      iii=lnblk(com3,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_3',NCCHAR,
     *            iii,com3,rcode)
c
c mark the end of file
c
      call ncendf(cdfout,rcode)
c 
c load the axes variable values
c
      xval(1) = 0.0
      srt(1)=1
      cnt(1)=1
      call ncvpt(cdfout,  xid, srt, cnt, xval, rcode)
      call ncvpt(cdfout,  yid, srt, cnt, xval, rcode)
      call ncvpt(cdfout,  zid, srt, cnt, xval, rcode)

c 
c sync the netcdf file
c
      call ncsnc(cdfout, rcode)

      return
      end
c
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

