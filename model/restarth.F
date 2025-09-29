      subroutine restart(n,froot,ifile,bcfile)
c
c   $Id: restart.F,v 1.11 1998/11/04 00:25:34 dwd Exp $
c 
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "moddef.inc"
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemi.inc"
#include "olemtke.inc"
#include "netcdf.inc"
#include "olemcdf.inc"
#include "cyclic.inc"
#include "refdate.inc"
#include "communicate.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif

      character ifile*80,froot*80,hold*2
      character rfile*80,bcfile*80,tunits*80,str*30
      integer cdfout,omode,iret

      integer lnblk,k,iii,ilen,i,n,tlen
      integer rcode,varid(4),varidu(4),varidw(4),varidv(4)
      integer varidr(4),countr(4),varidm(4),countm(4)
      integer countall(4)
      integer xudim,xwdim,ydim,zudim,zwdim,tdim,xdim,zdim
      integer xuid,xwid,yid,zuid,zwid,xid,zid,tid
      integer yudim,yuid,yvdim,yvid
      real xu(0:idim1),xw(0:idim1),zu(0:kdim1)
      real zw(0:kdim1),yv(0:jtotal1),yu(0:jtotal1)
      real tim
c
c temporary space for MPI i/o
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real phi2(0:idim1,0:jtotal1)
      integer startle(4),county(4),itrc,trcid(trcdim)
      integer prid,salid,thid,tkeid,rhoid,khid,kmid,rhopid
      integer umid,ummid,wmid,wmmid,vmid,vmmid
      integer hfid,lfid,strid,tbcid
#if defined(ICE)
      integer hicefid,ticefid,sicefid,vbrfid,cfrfid,icehfid,icesfid
#endif
      integer wallhfid,wallsfid
#if defined(HBOT)
      integer hbotfid
#endif
#if !defined(SNGLSALT)
c      real ssal(0:idim1)
c      integer k1,j1,j
#endif
      integer eplat,eplon,eppres,eptime
      data eplat,eplon,eppres,eptime/500,501,1,625/
c
      if(iope) then
c
c------Create cdf file (CDF).
      k=lnblk(froot,80)
      write(rfile,999)froot(1:k),int(n*delt+0.5)
 999  format(a,'_',i10,'_rst.nc')
      call squish(rfile,1,80)

c      cdfout=nccre(rfile,ncclob,rcode)
      iret = nf_create(rfile,IOR(NF_CLOBBER,NF_64BIT_OFFSET),
     *        cdfout)


      omode=ncsfil(cdfout,NCNOFILL,rcode)
c
      if(datemode) then
         call frmtdate(refdate, str)
         tunits = 'seconds' // ' since ' // str
         tlen = lnblk(tunits, 80)
      else
         tunits = 'seconds'
         tlen = 7
      endif
c
c------Define dimensions.
      xudim=ncddef(cdfout,'xu',idim+2,rcode)
      xuid =ncvdef(cdfout,'xu',NCFLOAT,1,xudim,rcode)
      call ncapt(cdfout, xuid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yudim=ncddef(cdfout,'yu',jtotal+2,rcode)
      yuid=ncvdef(cdfout,'yu',NCFLOAT,1,yudim,rcode)
      call ncapt(cdfout, yuid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c      
      zudim=ncddef(cdfout,'zu',kdim+2,rcode)
      zuid  = ncvdef(cdfout, 'zu', NCFLOAT, 1, zudim, rcode)
      call ncapt(cdfout, zuid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xwdim=ncddef(cdfout,'xw',idim+2,rcode)
      xwid  = ncvdef(cdfout, 'xw', NCFLOAT, 1, xwdim, rcode)
      call ncapt(cdfout, xwid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yvdim=ncddef(cdfout,'yv',jtotal+2,rcode)
      yvid  = ncvdef(cdfout, 'yv', NCFLOAT, 1, yvdim, rcode)
      call ncapt(cdfout, yvid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yvid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yvid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yvid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zwdim=ncddef(cdfout,'zw',kdim+2,rcode)
      zwid  = ncvdef(cdfout, 'zw', NCFLOAT, 1, zwdim, rcode)
      call ncapt(cdfout, zwid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xdim=ncddef(cdfout,'x',1,rcode)
      xid   = ncvdef(cdfout, 'x', NCFLOAT, 1, xdim, rcode)
      call ncapt(cdfout, xid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      ydim=ncddef(cdfout,'y',1,rcode)
      yid   = ncvdef(cdfout, 'y', NCFLOAT, 1, ydim, rcode)
      call ncapt(cdfout, yid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zdim=ncddef(cdfout,'z',1,rcode)
      zid   = ncvdef(cdfout, 'z', NCFLOAT, 1, zdim, rcode)
      call ncapt(cdfout, zid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      tdim=ncddef(cdfout,'time',NCUNLIM,rcode)
      tid=ncvdef(cdfout, 'time', NCLONG, 1, tdim, rcode)
      if(.not.datemode)then
         call ncapt(cdfout, tid, 'epic_code', NCLONG, 1, 
     *        EPTIME, rcode)
      endif
      call ncaptc(cdfout, tid, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, tid, 'units', NCCHAR, tlen, tunits,rcode)
      call ncaptc(cdfout, tid, 'type' , NCCHAR, 1, ' ', rcode)

      tim=n*delt


c------Define cell-centered variables (CDF).
      varid(1)=xwdim
      varid(2)=yudim
      varid(3)=zudim
      varid(4)=tdim

      varidu(1)=xudim
      varidu(2)=yudim
      varidu(3)=zudim
      varidu(4)=tdim

      varidw(1)=xwdim
      varidw(2)=yudim
      varidw(3)=zwdim
      varidw(4)=tdim

      varidv(1)=xwdim
      varidv(2)=yvdim
      varidv(3)=zudim
      varidv(4)=tdim

      varidr(1)=xdim
      varidr(2)=ydim
      varidr(3)=zudim
      varidr(4)=tdim

      varidm(1)=xwdim
      varidm(2)=yudim
      varidm(3)=zdim
      varidm(4)=tdim

c------Define variables and set attributes.
      call setfld3(MDUM,umid,varidu,cdfout,4)
c      call setfld3(MDUMM,ummid,varidu,cdfout,4)

      call setfld3(MDWM,wmid,varidw,cdfout,4)
c      call setfld3(MDWMM,wmmid,varidw,cdfout,4)

      call setfld3(MDVM,vmid,varidv,cdfout,4)
c      call setfld3(MDVMM,vmmid,varidv,cdfout,4)

#if defined(LEAPF)
      call setfld3(MDPRES,prid,varid,cdfout,4)
#endif
      call setfld3(MDTEMP,thid,varid,cdfout,4)
      call setfld3(MDSAL,salid,varid,cdfout,4)
c      call setfld3(MDTKE,tkeid,varid,cdfout,4)
c      call setfld3(MDKH,khid,varid,cdfout,4)
c      call setfld3(MDKM,kmid,varid,cdfout,4)
      call setfld3(MDRHO,rhopid,varid,cdfout,4)
      call setfld3(MDRHOBM1,rhoid,varidr,cdfout,4)

      if(hfton)call setfld3(MDHFTMASK,hfid,varidm,cdfout,4)
      if(lfton)call setfld3(MDLFTMASK,lfid,varidm,cdfout,4)
      if(strton) then
         call setfld3(MDSTRTMASK,strid,varidm,cdfout,4)
         call setfld4(MDSTRTBC,tbcid,varidm,cdfout,4)
      endif

      do 20 i=1,ntrace
         call setfld3(MDTRCR+i,trcid(i),varid,cdfout,4)
 20   continue
#if defined(ICEWALL) || defined(MPOND)
       call setfld3(MDWHF,wallhfid,varid,cdfout,4)
       call setfld3(MDWSF,wallsfid,varid,cdfout,4)
#endif
#if defined(ICE)
       call setfld3(MDHICE,hicefid,varidm,cdfout,4)
       call setfld3(MDTICE,ticefid,varidm,cdfout,4)
       call setfld3(MDSICE,sicefid,varidm,cdfout,4)
c       call setfld3(MDVBR,vbrfid,varidm,cdfout,4)
#if !defined(WARMICE)
       call setfld3(MDCFR,cfrfid,varid,cdfout,4)
#endif
c       call setfld3(MDFICE,icehfid,varidm,cdfout,4)
c       call setfld3(MDSFICE,icesfid,varid,cdfout,4)
#endif
#if defined(HBOT)
       call setfld3(MDHBOT,hbotfid,varidm,cdfout,4)
#endif

      iii=lnblk(ident,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'creation_date',NCCHAR,
     *     iii,ident,rcode)
      iii=lnblk(mver,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'version',NCCHAR,iii,
     *     mver,rcode)

      do i=1,4
         iii=lnblk(modflags(i),80)
         if(iii.gt.0) then
            write(hold,'(i1)')i
            call ncaptc(cdfout,NCGLOBAL,'model_flags'//hold,
     *           NCCHAR, iii, modflags(i), rcode)
         endif
      enddo

      iii=lnblk(title,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'title',NCCHAR,iii,
     *     title,rcode)
      call ncaptc(cdfout,NCGLOBAL,'DATA_TYPE',NCCHAR,5,'MODEL',rcode)
      call ncaptc(cdfout,NCGLOBAL,'COORD_SYSTEM',NCCHAR,5,'LOCAL',rcode)
      call ncaptc(cdfout,NCGLOBAL,'FILE_TYPE',NCCHAR,7,'RESTART',rcode)
      ilen=lnblk(ifile,80)
      if(ilen.gt.0)call ncaptc(cdfout,NCGLOBAL,'IC_FILE_NAME',NCCHAR,
     *     ilen,ifile,rcode)
      ilen=lnblk(bcfile,80)
      if(ilen.gt.0)call ncaptc(cdfout,NCGLOBAL,'BC_FILE_NAME',NCCHAR,
     *     ilen,bcfile,rcode)
      call ncapt(cdfout,NCGLOBAL,'dx',NCFLOAT,1,dx,rcode)
      call ncapt(cdfout,NCGLOBAL,'dy',NCFLOAT,1,dy,rcode)
      call ncapt(cdfout,NCGLOBAL,'dz',NCFLOAT,1,dz,rcode)
      call ncapt(cdfout,NCGLOBAL,'imax',NCLONG,1,idim,rcode)
      call ncapt(cdfout,NCGLOBAL,'jmax',NCLONG,1,jtotal,rcode)
      call ncapt(cdfout,NCGLOBAL,'kmax',NCLONG,1,kdim,rcode)
      call ncapt(cdfout,NCGLOBAL,'field_time',NCFLOAT,1,tim,rcode)
      call ncapt(cdfout,NCGLOBAL,'ncool',NCLONG,1,ncool,rcode)
      call ncapt(cdfout,NCGLOBAL,'heatf',NCFLOAT,1,heatf,rcode)
      iii=lnblk(com1,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_1',NCCHAR,
     *                         iii,com1,rcode)
      iii=lnblk(com2,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_2',NCCHAR,
     *                         iii,com2,rcode)
      iii=lnblk(com3,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_3',NCCHAR,
     *                         iii,com3,rcode)
      call ncapt(cdfout,NCGLOBAL,'ntrace',NCLONG,1,ntrace,rcode)
      do 10 itrc=1,ntrace
        write(hold,'(i1)')itrc
        iii=lnblk(trccom(itrc),120)
        if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'trccom_'//hold,
     *       NCCHAR,iii,trccom(itrc),rcode)
        iii=lnblk(trfile(itrc),80)
        if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'trcfile_'//hold,
     *       NCCHAR,iii,trfile(itrc),rcode)
10    continue
      iii=lnblk(sgscom1,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'sgscom_1',NCCHAR,
     *     iii,sgscom1,rcode)
      iii=lnblk(sgscom2,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'sgscom_2',NCCHAR,
     *     iii,sgscom2,rcode)
      iii=lnblk(sgscom3,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'sgscom_3',NCCHAR,
     *     iii,sgscom3,rcode)
      iii=lnblk(bccom1,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'bccom_1',NCCHAR,
     *     iii,bccom1,rcode)
      iii=lnblk(bccom2,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'bccom_2',NCCHAR,
     *     iii,bccom2,rcode)
      iii=lnblk(bccom3,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'bccom_3',NCCHAR,
     *     iii,bccom3,rcode)
      call ncapt(cdfout,NCGLOBAL,'mdflgs',NCLONG,20,mdflgs,rcode)
      call ncapt(cdfout,NCGLOBAL,'mdvals',NCFLOAT,20,mdvals,rcode)
      call ncapt(cdfout,NCGLOBAL,'notimes',NCLONG,1,notimes,rcode)
      call ncapt(cdfout,NCGLOBAL,'nofld',NCLONG,1,nofld,rcode)
      call ncapt(cdfout,NCGLOBAL,'nostart',NCLONG,1,nostart,rcode)
      call ncapt(cdfout,NCGLOBAL,'nostop',NCLONG,1,nostop,rcode)
      call ncapt(cdfout,NCGLOBAL,'nodt',NCLONG,1,nodt,rcode)
      call ncapt(cdfout,NCGLOBAL,'ndiag',NCLONG,1,ndiag,rcode)
      call ncapt(cdfout,NCGLOBAL,'diloc',NCLONG,1,diloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'djloc',NCLONG,1,djloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'dkloc',NCLONG,1,dkloc,rcode)
      call ncendf(cdfout,rcode)

c------Build the axes.
      do 300 i=0,idim1
        xu(i)=dx*float(i-1)
        xw(i)=dx*float(i-1)+dx*0.5
 300  continue
#if defined(RADTOP)
      do 310 i=0,kdim1
c         zu(i)=dz*float(kdim-i) + dz*0.5
         zu(i)=dz*float(i-1) + dz*0.5
c         zw(i)=dz*float(kdim-i+1)
         zw(i)=dz*float(i-1)
 310  continue
#else
      do 310 i=0,kdim1
         zu(i)=dz*float(kdim-i) + dz*0.5
         zw(i)=dz*float(kdim-i+1)
#endif
 310  continue
      do 320 i=0,jtotal+1
         yv(i)=dy*float(i-1)
         yu(i)=yv(i)+dy*0.5
320   continue

      yv(jtotal+1)=dy*jtotal

c------Load the axes variable values.

      call ncvpt(cdfout, xuid, 1, idim+2, xu, rcode)
      call ncvpt(cdfout, xwid, 1, idim+2, xw, rcode)
      call ncvpt(cdfout, yuid, 1, jtotal+2, yu, rcode)
      call ncvpt(cdfout, yvid, 1, jtotal+2, yv, rcode)
      call ncvpt(cdfout, zuid, 1, kdim+2, zu, rcode)
      call ncvpt(cdfout, zwid, 1, kdim+2, zw, rcode)
      call ncvpt(cdfout, xid, 1, 1, 0., rcode)
      call ncvpt(cdfout, yid, 1, 1, 0., rcode)
      call ncvpt(cdfout, zid, 1, 1, 0., rcode)
      call ncvpt(cdfout, tid, 1, 1,int(tim), rcode)

      county(1)=idim+2
      county(2)=1 
      county(3)=1 
      county(4)=1

      countr(1)=1
      countr(2)=1
      countr(3)=kdim+2
      countr(4)=1


      countm(1)=idim+2
      countm(2)=jtotal+2
      countm(3)=1
      countm(4)=1

      startle(1)=1
      startle(2)=1
      startle(3)=1
      startle(4)=1
c
c write fields directly
c
      countall(1)=idim+2
      countall(2)=jtotal+2
      countall(3)=kdim+2
      countall(4)=1
c
      endif
c
      call write3d(cdfout,umid,startle,countall,um,phi)
c      call write3d(cdfout,ummid,startle,countall,umm,phi)
      call write3d(cdfout,vmid,startle,countall,vm,phi)
c      call write3d(cdfout,vmmid,startle,countall,vmm,phi)
      call write3d(cdfout,wmid,startle,countall,wm,phi)
c      call write3d(cdfout,wmmid,startle,countall,wmm,phi)
      call write3d(cdfout,thid,startle,countall,th,phi)
#if defined(LEAPF)
      call write3d(cdfout,prid,startle,countall,p,phi)
#endif
c      call write3d(cdfout,tkeid,startle,countall,tke,phi)
c      call write3d(cdfout,khid,startle,countall,kheat,phi)
c      call write3d(cdfout,kmid,startle,countall,kmom,phi)
      call write3d(cdfout,rhopid,startle,countall,rhop,phi)
c
c output bottom terrain
c
#if defined(HBOT)
      countall(3) = 1
      call write2d(cdfout,hbotfid,startle,countall,hbot,phi2)
      countall(3)=kdim+2
#endif
#if defined(ICEWALL) || defined(MPOND)
      call write3d(cdfout,wallhfid,startle,countall,wallhf,phi)
      call write3d(cdfout,wallsfid,startle,countall,wallsf,phi)
#endif

#if defined(ICE)
c
#if !defined(WARMICE)
      call write3d(cdfout,cfrfid,startle,countall,cfr,phi)
#endif
c
c change vertical count to 1
c
c      countall(3) = 1
c
c output ice thickness
c
      call write2d(cdfout,hicefid,startle,countm,hice,phi2)
c
c output ice temperature
c
      call write2d(cdfout,ticefid,startle,countm,tice,phi2)
c
c output ice salinity
c
      call write2d(cdfout,sicefid,startle,countm,sice,phi2)
c
c output brine volume
c
c      call write2d(cdfout,vbrfid,startle,countm,vbr,phi2)
c
c output heat flux from ice into upper ocean
c this is not needed for restart, but is for reference
c
c       call write2d(cdfout,icehfid,startle,countm,iceh,phi2)
c
c
c
c output salt flux from ice into ocean
c this is not needed for restart, but is for reference
c
c      call imag2d(icesal)
c      call write2d(cdfout,icesfid,startle,countm,icesal,phi)
c
c change vertical count back kdim+2
c
c      countall(3) = kdim+2
c
#endif
c
#if defined(TRACER)
      do i=1,ntrace
         call write3d(cdfout,trcid(i),startle,countall,trc(0,0,0,i),phi)
      enddo
#endif
#if defined(SNGLSALT)
      call write3d(cdfout,salid,startle,countall,sal,phi)
#else 
      call write3ddp(cdfout,salid,startle,countall,sal,phi)
c
c do not have double precision salt with MPI
c
c      do k=1,kdim+2
c         startle(3)=k
c         k1=k-1
c         do j=1,jdim+2
c            startle(2)=j
c            j1=j-1
c            do i=0,idim1
c              ssal(i)=sal(i,j1,k1)
c            enddo
c            call ncvpt(cdfout, salid,startle,county,  ssal,rcode)
c            
c         enddo
c      enddo
#endif

C-------For each mask, write hyperslabs of values into its variable.
      startle(3)=1
      startle(4)=1

      if(hfton) then
        call write2d(cdfout,hfid,startle,countm,hf_tmask,phi2)
      endif
      if(lfton) then
        call write2d(cdfout,lfid,startle,countm,lf_tmask,phi2)
      endif
      if(strton) then
        call write2d(cdfout,strid,startle,countm,str_tmask,phi2)
        call write2d(cdfout,tbcid,startle,countm,str_tbc,phi2)
      endif
      if(iope) then     
        startle(1)=1
        startle(2)=1
        startle(3)=1
        startle(4)=1
        call ncvpt(cdfout,rhoid,startle,countr,rhobm1(0),rcode)
  
        call ncclos(cdfout,rcode)
      endif

      return
      end


c***********************************************************************
      subroutine setfld3(ivarib,id,varid,cdfout,ndim)

#include "netcdf.inc"

      integer ivarib,id,varid(4),ndim,cdfout,lnblk,rcode,ilen
      character name*15

      call varnme(ivarib,name)
      ilen=lnblk(name,15)
      id=ncvdef(cdfout,name(:ilen),NCFLOAT,ndim,varid,rcode)
      call setvar(ivarib,cdfout,id)
      return
      end

c***********************************************************************
      subroutine setfld4(ivarib,id,varid,cdfout,ndim)

#include "netcdf.inc"

      integer ivarib,id,varid(4),ndim,cdfout,lnblk,rcode,ilen
      character name*15

      call varnme(ivarib,name)
      ilen=lnblk(name,15)
      id=ncvdef(cdfout,name(:ilen),NCLONG,ndim,varid,rcode)
      call setvar(ivarib,cdfout,id)
      return
      end

#if !defined(SNGLSALT)
      subroutine write3ddp(cdfout,cdfid,start,count,fld,phi)
c
c this routine outputs a three dimensional field including the
c image points (for use in restart). If using MPI, it does a
c collection from all processors.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      real*8 fld(0:idim1,0:jdim1,0:kdim1)
      real phi(0:idim1,0:jtotal1,0:kdim1)
      integer cdfout,cdfid,start(4),count(4),rcode
      integer i,j,k
#if defined(USEMPI)
      real*8 phi8(0:idim1,0:jtotal1,0:kdim1)
      integer nerr
#endif
c
#if defined(USEMPI)
      call MPI_BARRIER(comm,nerr)
      call collect_3ddp(fld,phi8)
      do k=0,kdim1
         do j=0,jtotal1
            do i=0,idim1
               phi(i,j,k) = phi8(i,j,k)
            enddo
         enddo
      enddo
      if(iope) then
        call ncvpt(cdfout,  cdfid, start, count,   phi,rcode)
      endif
#else
      do k=0,kdim1
         do j=0,jdim1
            do i=0,idim1
               phi(i,j,k) = fld(i,j,k)
            enddo
         enddo
      enddo
      call ncvpt(cdfout,  cdfid, start, count, phi,rcode)
#endif
      return
      end
c
#endif
c
      subroutine write2d(cdfout,cdfid,start,count,fld,phi)
c
c this routine outputs a two dimensional field including the
c image points (for use in restart). If using MPI, it does a
c collection from all processors.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
c
      real fld(0:idim1,0:jdim1)
      real phi(0:idim1,0:jtotal1)
      integer cdfout,cdfid,start(4),count(4),rcode
#if defined(USEMPI)
      integer nerr
#endif
c
#if defined(USEMPI)
      call MPI_BARRIER(comm,nerr)
      call collect_2d(fld,phi)
      if(iope) then
        call ncvpt(cdfout,  cdfid, start, count,   phi,rcode)
      endif
#else
      call ncvpt(cdfout,  cdfid, start, count,   fld,rcode)
#endif
      return
      end
c
c
      subroutine write3d(cdfout,cdfid,rstart,rcount,fld,phi)
c
c Routine to write out netcdf file in chunks to avoid large array
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
c

      real fld(0:idim1,0:jdim1,0:kdim1)
      real phi(0:idim1,0:jdim1,0:kdim1)
      integer cdfout,cdfid,start(4),count(4)
      integer rstart(4),rcount(4),rcode
      integer np,ier,recv_req(nprocs),send_req,k,nprec,nps
      integer recv_status(MPI_STATUS_SIZE,nprocs),
     &             status(MPI_STATUS_SIZE)


c
      count(1) = idim+2
      count(2) = jdim+2
      count(3) = kdim+2
      count(4) = 1
      if(iope) then
c        write(*,*)"Processing data chunks on iope"
        do np=1,nprocs
          start(1) = 1
          if(np.eq.1) then
            start(2) = 1
          else
            start(2) = jdim*(np-1)+1
          endif
          start(3) = 1
          start(4) = 1
          do nps=1,nprocs
            if(blockcord(nps).eq.np-1) then
               nprec = nps-1
            endif
          enddo
c          write(*,*)"Post Recv ",nprec
          if(nprec.ne.0)then
            call MPI_RECV(phi(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $         nprec,mpitag_io,comm,status,ier)
            call ncvpt(cdfout, cdfid, start, count, phi, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite", nprec
          else
c            write(*,*)"writing iope chunk"
            call ncvpt(cdfout, cdfid, start, count, fld, rcode)
            if(rcode.ne.0) write(*,*)"Error in ncwrite on iope"
          endif            
c          write(*,*)"Output phi ",phi(1,1,1),nprec,start(2)
        enddo
      else
c        write(*,*)"Post Sends ",my_pe
        call MPI_SEND(fld(0,0,0),idim2*jdim2*kdim2,MPI_REAL,
     $               mastertask, mpitag_io,comm,ier)
      endif
      call MPI_BARRIER(comm,ier)
      return
      end
c
c
      subroutine read3d(cdfin,cdfid,rstart,rcount,fld,phi)
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      integer idimh,jdimh,kdimh,idimh1,jdimh1,kdimh1
      integer idimh2,jdimh2,kdimh2
      parameter(idimh=idim/2,jdimh=jdim/2,kdimh=kdim/2,
     *          idimh1=idimh+1,jdimh1=jdimh+1,kdimh1=kdimh+1,
     *          idimh2=idimh+2,jdimh2=jdimh+2,kdimh2=kdimh+2)
      real fld(0:idimh1,0:jdimh1,0:kdimh1)
      real phi(0:idimh1,0:jdimh1,0:kdimh1)
      integer cdfin,cdfid,start(4),count(4),rcode
      integer rstart(4),rcount(4)
      integer np,ier,recv_req, send_req,i,k,nprec,nps
      integer status(MPI_STATUS_SIZE)


c
      count(1) = idimh+2
      count(2) = jdimh+2
      count(3) = kdimh+2
      count(4) = 1
      if(iope) then
        do np=1,nprocs
          start(1) = 1
          if(np.eq.1) then
            start(2) = 1
          else
            start(2) = (jdimh)*(np-1)+1
          endif
          start(3) = 1
          start(4) = 1
          do nps=1,nprocs
            if(blockcord(nps).eq.np-1) then
               nprec = nps-1
            endif
          enddo
          if(nprec.ne.0) then
c            write(*,*)"getting phi ",nprec,count(1),count(2),count(3),
c     *          start(2)
            call ncvgt(cdfin,cdfid,start,count,phi,rcode)
            call MPI_SEND(phi(0,0,0),idimh2*jdimh2*kdimh2,MPI_REAL,
     $         nprec,
     $         mpitag_io,comm,ier)
          else
c            write(*,*)"getting fld ",nprec,count(1),count(2),count(3),
c     *          start(2)
            call ncvgt(cdfin,cdfid,start,count,fld,rcode)
          endif
        enddo
      else
        call MPI_RECV(fld(0,0,0),idimh2*jdimh2*kdimh2,MPI_REAL,
     $                mastertask,mpitag_io,comm,status,ier)
      endif
      call MPI_BARRIER(comm,ier)
c      write(*,*)"my pe ",my_pe,phi(1,1,1)
      return
      end



