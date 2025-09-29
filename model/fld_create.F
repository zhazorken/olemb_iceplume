      subroutine fld_create(ifile, bcfile, froot, ko)
c
c   $Id: fld_create.F,v 1.7 1998/10/15 16:47:18 dwd Exp $
c
#include "moddef.inc"
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "olemi.inc"
#include "cyclic.inc"
#include "olemtke.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
#include "refdate.inc"
c 
      character ifile*80,froot*80,temp*80,name*15,str*30
      character hold*2,bcfile*80,tunits*80
      integer ko,cdfout
      integer lnblk,k,iii,j,ilen,i,tlen,ival
      integer rcode,varid(4),varidu(4),varidw(4),varidp(4)
      integer varidm(4),countm(4),strtm(4)
      integer varidv(4)
      integer xudim,xwdim,ydim,zudim,zwdim,tdim,xdim,zdim
      integer xuid,xwid,yid,zuid,zwid,xid,zid
      integer yudim,yvdim,yuid,yvid
      integer hffid,lffid,strfid,strtbcfid,t1dim,t1id
      real xu(idim1),xw(idim),zu(kdim),zw(kdim1),yu(jtotal),yv(jtotal1)
#if defined(ALLOC_DISK)
      integer tstart,tstop,tpoints
#else
      integer omode
#endif
      integer eplat,eplon,eppres,eptime
      data eplat,eplon,eppres,eptime/500,501,1,625/
c
C-----create cdf file (CDF)
c
      ilen=lnblk(froot,80)
      if(ko.eq.1)then
c         fldid(ko)=nccre(froot(1:ilen)//'_fld.nc',ncclob,rcode)
        rcode = nf_create(froot(1:ilen)//'_fld.nc',
     *        IOR(NF_CLOBBER,NF_64BIT_OFFSET),fldid(ko))
      else
         write(temp,999)froot(1:ilen),ko
 999     format(a,'_',i1)
         ilen=lnblk(temp,80)
c         fldid(ko)=nccre(temp(1:ilen)//'_fld.nc',ncclob,rcode)
        rcode = nf_create(froot(1:ilen)//'_fld.nc',
     *        IOR(NF_CLOBBER,NF_64BIT_OFFSET),fldid(ko))

      endif
#if !defined(ALLOC_DISK)
c
c don't use fills
c
      omode=ncsfil(fldid(ko),NCNOFILL,rcode)
#endif
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
C-----define dimensions
c
      cdfout=fldid(ko)
      xudim=ncddef(cdfout,'xu',idim1,rcode)
      xuid =ncvdef(cdfout,'xu',NCFLOAT,1,xudim,rcode)
      call ncapt(cdfout, xuid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yudim=ncddef(cdfout,'yu',jtotal,rcode)
      yuid  = ncvdef(cdfout, 'yu', NCFLOAT, 1, yudim, rcode)
      call ncapt(cdfout, yuid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zudim=ncddef(cdfout,'zu',kdim,rcode)
      zuid  = ncvdef(cdfout, 'zu', NCFLOAT, 1, zudim, rcode)
      call ncapt(cdfout, zuid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zuid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zuid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zuid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xwdim=ncddef(cdfout,'xw',idim,rcode)
      xwid  = ncvdef(cdfout, 'xw', NCFLOAT, 1, xwdim, rcode)
      call ncapt(cdfout, xwid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yvdim=ncddef(cdfout,'yv',jtotal1,rcode)
      yvid  = ncvdef(cdfout, 'yv', NCFLOAT, 1, yvdim, rcode)
      call ncapt(cdfout, yvid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yvid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yvid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yvid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zwdim=ncddef(cdfout,'zw',kdim1,rcode)
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
      t1dim=ncddef(cdfout,'t1',1,rcode)
      t1id = ncvdef(cdfout, 't1', NCLONG, 1, t1dim, rcode)
      if(.not.datemode) then
         call ncapt(cdfout, t1id, 'epic_code', NCLONG, 1, 
     *        EPTIME, rcode)
      endif
      call ncaptc(cdfout, t1id, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, t1id, 'units', NCCHAR, tlen, tunits, rcode)
      call ncaptc(cdfout, t1id, 'type' , NCCHAR, 1, ' ', rcode)
c
#if defined(ALLOC_DISK)
      tstart=max0(nostart(ko),nstart-1)
      tstop=min0(nostop(ko),ntime)
      tpoints=(tstop-tstart)/nodt(ko)+1
      write(*,*)'Number of Records ',tpoints
      tdim=ncddef(cdfout,'time',tpoints,rcode)
#else
      tdim=ncddef(cdfout,'time',NCUNLIM,rcode)
#endif
      tfid(ko)   = ncvdef(cdfout, 'time', NCLONG, 1, tdim, rcode)
      if(.not.datemode) then
         call ncapt(cdfout, tfid(ko), 'epic_code', NCLONG, 1, 
     *        EPTIME, rcode)
      endif
      call ncaptc(cdfout, tfid(ko), 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, tfid(ko), 'units', NCCHAR,tlen,tunits,rcode)
      call ncaptc(cdfout, tfid(ko), 'type' , NCCHAR, 1, ' ', rcode)
c
C-----define cell-centered variables (CDF)
      varid(1)=xwdim
      varid(2)=yudim
      varid(3)=zudim
      varid(4)=tdim

C-----define u variable (shifted to left of center) (CDF)
      varidu(1)=xudim
      varidu(2)=yudim
      varidu(3)=zudim
      varidu(4)=tdim

C-----define v variable
      varidv(1)=xwdim
      varidv(2)=yvdim
      varidv(3)=zudim
      varidv(4)=tdim

C-----define w variable (shifted to bottom of center) (CDF)
      varidw(1)=xwdim
      varidw(2)=yudim
      varidw(3)=zwdim
      varidw(4)=tdim

C-----define non-grid variables
      varidp(1)=xdim
      varidp(2)=ydim
      varidp(3)=zdim
      varidp(4)=tdim

C-----define mask variables
      varidm(1)=xwdim
      varidm(2)=yudim
      varidm(3)=zdim
      varidm(4)=t1dim

C---- set attributes
      do 980 j=1,nofld(ko)
         call varnme(ofield(ko,j),name)
         ilen=lnblk(name,15)
         if(ofield(ko,j).eq.MDUVEL) then
c
C----     set attributes for u
c
            ufid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidu,rcode)
            call setvar(ofield(ko,j),cdfout,ufid(ko))
         elseif(ofield(ko,j).eq.MDVVEL) then
c
C----     set attributes for v
c
            vfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidv,rcode)
            call setvar(ofield(ko,j),cdfout,vfid(ko))
c
         elseif(ofield(ko,j).eq.MDWVEL) then
c
C----    set attributes for w
c
            wfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidw,rcode)
            call setvar(ofield(ko,j),cdfout,wfid(ko))
         elseif(ofield(ko,j).eq.MDPRES) then
c
C----    set attributes for p_p
c
            pfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,pfid(ko))
         elseif(ofield(ko,j).eq.MDTEMP) then
c
C----    set attributes for t
c
            thfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,thfid(ko))
         elseif(ofield(ko,j).eq.MDSAL) then
c
C----    set attributes for s
c
            salfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,salfid(ko))
         elseif(ofield(ko,j).eq.MDTKE) then
c
C----    set attributes for tke
c
            tkefid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,tkefid(ko))
         elseif(ofield(ko,j).eq.MDRHO) then
c
C----    set attributes for rho_p
c
            rho_pfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,
     *           rcode)
            call setvar(ofield(ko,j),cdfout,rho_pfid(ko))
         else if(ofield(ko,j).eq.MDKH)then
c
c set attributes for kh
c
            khfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,khfid(ko))
         else if(ofield(ko,j).eq.MDKM)then
c
c set attributes for km
c
            kmfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,kmfid(ko))
         else if(ofield(ko,j).eq.MDRHOT)then
c
c set attributes for rho
c
            rhofid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,rhofid(ko))
         else if(ofield(ko,j).eq.MDHFSGS)then
c
c set attrigutes for hf_sgs  (SGS heat flux)
c
            hf_sgsid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidw,
     *           rcode)
            call setvar(ofield(ko,j),cdfout,hf_sgsid(ko))
         else if(ofield(ko,j).gt.MDTRCR)then
            k=ofield(ko,j)-MDTRCR
            write(str,994)k
 994        format('trcr',i1)
            trcrfid(ko,k)=ncvdef(cdfout,str,NCFLOAT,4,varid,rcode)
            call setvar(ofield(ko,j),cdfout,trcrfid(ko,k))
         endif
 980  continue
c
C----    set attributes for hf_tmask
c
      if(hfton)then
         call varnme(MDHFTMASK,name)
         ilen=lnblk(name,15)
         hffid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDHFTMASK,cdfout,hffid)
      endif
c
C----    set attributes for lf_tmask
c
      if(lfton)then
         call varnme(MDLFTMASK,name)
         ilen=lnblk(name,15)
         lffid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDLFTMASK,cdfout,lffid)
      endif
c
C----    set attributes for str_tmask
c
      if(strton)then
         call varnme(MDSTRTMASK,name)
         ilen=lnblk(name,15)
         strfid=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidm,rcode)
         call setvar(MDSTRTMASK,cdfout,strfid)
c
C----    set attributes for str_tbc
c
         call varnme(MDSTRTBC,name)
         ilen=lnblk(name,15)
         strtbcfid=ncvdef(cdfout,name(:ilen),NCLONG,4,varidm,rcode)
         call setvar(MDSTRTBC,cdfout,strtbcfid)
      endif
c
c set attributes for non-grid variable (allways created)
c
c set attributes for mixed layer depth
c
      call varnme(MDZI,name)
      ilen=lnblk(name,15)
      zifid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
      call setvar(MDZI,cdfout,zifid(ko))
c
c set attributes for w_star
c
      call varnme(MDWSTAR,name)
      ilen=lnblk(name,15)
      w_starfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
      call setvar(MDWSTAR,cdfout,w_starfid(ko))
c
c set attributes for u_star
c
      call varnme(MDUSTAR,name)
      ilen=lnblk(name,15)
      u_starfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
      call setvar(MDUSTAR,cdfout,u_starfid(ko))
c
c set attributes for t_star
c
      call varnme(MDTSTAR,name)
      ilen=lnblk(name,15)
      t_starfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
      call setvar(MDTSTAR,cdfout,t_starfid(ko))
c
c set attributes for q
c
      call varnme(MDQ,name)
      ilen=lnblk(name,15)
      qfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,4,varidp,rcode)
      call setvar(MDQ,cdfout,qfid(ko))
c
c set attributes for boundary condition variables
c
c      if(bcon) then
c
c set attributes for hf_top
c
         if(heaton) then
            call varnme(MDHFT,name)
            ilen=lnblk(name,15)
            hf_topfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDHFT,cdfout,hf_topfid(ko))
         endif
c
c set attributes for ustr_t
c
         if(uston) then
            call varnme(MDUSTT,name)
            ilen=lnblk(name,15)
            ustr_tfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDUSTT,cdfout,ustr_tfid(ko))
         endif
c
c set attributes for vstr_t
c
         if(vston) then
            call varnme(MDVSTT,name)
            ilen=lnblk(name,15)
            vstr_tfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDVSTT,cdfout,vstr_tfid(ko))
         endif
c
c set attributes for swf_top
c
         if(swheaton) then
            call varnme(MDSWHFT,name)
            ilen=lnblk(name,15)
            swf_topfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDSWHFT,cdfout,swf_topfid(ko))
         endif
c
c set attributes for rain
c
         if(rainon) then
            call varnme(MDRAIN,name)
            ilen=lnblk(name,15)
            rainfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDRAIN,cdfout,rainfid(ko))
         endif
c
c set attributes for lhf_top
c
         if(lheaton) then
            call varnme(MDLHFT,name)
            ilen=lnblk(name,15)
            lhf_topfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDLHFT,cdfout,lhf_topfid(ko))
         endif
#if defined(STOKES)
c
c set attributes for wave_l
c
         if(wavelon) then
            call varnme(MDWAVEL,name)
            ilen=lnblk(name,15)
            wave_lfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDWAVEL,cdfout,wave_lfid(ko))
         endif
c
c set attributes for wave_h
c
         if(wavehon) then
            call varnme(MDWAVEH,name)
            ilen=lnblk(name,15)
            wave_hfid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDWAVEH,cdfout,wave_hfid(ko))
         endif
c
c set attributes for w_angle
c
         if(angleon) then
            call varnme(MDANGLE,name)
            ilen=lnblk(name,15)
            w_anglefid(ko)=ncvdef(cdfout,name(:ilen),NCFLOAT,
     *           4,varidp,rcode)
            call setvar(MDANGLE,cdfout,w_anglefid(ko))
         endif
#endif
c      endif
c
      iii=lnblk(ident,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'identifier',NCCHAR,
     *            iii,ident,rcode)

      iii=lnblk(mver,80)
      if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'model_version',NCCHAR,
     *            iii,mver,rcode)

      do i=1,4
         iii=lnblk(modflags(i),80)
         if(iii.gt.0) then
            write(hold,'(i1)')i
            call ncaptc(cdfout,NCGLOBAL,'model_flags'//hold,
     *           NCCHAR, iii, modflags(i), rcode)
         endif
      enddo

      iii=lnblk(title(1),80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'title',NCCHAR,
     *            iii,title(1),rcode)
      do i=2,ntitle
         iii=lnblk(title(i),80)
         write(hold,'(i1)')i
         if(iii.gt.0)call ncaptc(cdfout,NCGLOBAL,'title_'//hold,
     *        NCCHAR,iii,title(i),rcode)
      enddo
      call ncaptc(cdfout,NCGLOBAL,'DATA_TYPE',NCCHAR,
     *     5,'MODEL',rcode)
      call ncaptc(cdfout,NCGLOBAL,'COORD_SYSTEM',NCCHAR,
     *     5,'LOCAL',rcode)
      call ncapt(cdfout,NCGLOBAL,'dx',NCFLOAT,1,dxx,rcode)
      call ncapt(cdfout,NCGLOBAL,'dy',NCFLOAT,1,dy,rcode)
      call ncapt(cdfout,NCGLOBAL,'dz',NCFLOAT,1,dz,rcode)
      ival = idim
      call ncapt(cdfout,NCGLOBAL,'imax',NCLONG,1,ival,rcode)
      ival = jtotal
      call ncapt(cdfout,NCGLOBAL,'jmax',NCLONG,1,ival,rcode)
      ival = kdim
      call ncapt(cdfout,NCGLOBAL,'kmax',NCLONG,1,ival,rcode)
      call ncapt(cdfout,NCGLOBAL,'latitude',NCFLOAT,1,xlat,rcode)
c 8-21-92 jrb  "ntrace" added to cdf output file
      call ncapt(cdfout,NCGLOBAL,'ntrace',NCLONG,1,ntrace,rcode)
      call ncapt(cdfout,NCGLOBAL,'time_start',NCFLOAT,1,start,rcode)
      call ncapt(cdfout,NCGLOBAL,'time_stop',NCFLOAT,1,stop,rcode)
      call ncapt(cdfout,NCGLOBAL,'time_step',NCFLOAT,1,delt,rcode)
      call ncapt(cdfout,NCGLOBAL,'restart_int',NCLONG,1,nrstrt,rcode)
      call ncapt(cdfout,NCGLOBAL,'mdflgs',NCLONG,20,mdflgs,rcode)
      call ncapt(cdfout,NCGLOBAL,'mdvals',NCFLOAT,20,mdvals,rcode)
      iii=lnblk(sgscom1,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'sgscom1',NCCHAR,
     *            iii,sgscom1,rcode)
      iii=lnblk(sgscom2,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'sgscom2',NCCHAR,
     *            iii,sgscom2,rcode)
      iii=lnblk(sgscom3,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'sgscom3',NCCHAR,
     *            iii,sgscom3,rcode)
      call ncapt(cdfout,NCGLOBAL,'sgsflgs',NCLONG,20,sgsflgs,rcode)
      call ncapt(cdfout,NCGLOBAL,'sgsvals',NCFLOAT,20,sgsvals,rcode)
      call ncapt(cdfout,NCGLOBAL,'notimes',NCLONG,1,notimes,rcode)
      call ncapt(cdfout,NCGLOBAL,'field',NCLONG,1,ko,rcode)
      call ncapt(cdfout,NCGLOBAL,'nofld',NCLONG,5,nofld,rcode)
      call ncapt(cdfout,NCGLOBAL,'ofield',NCLONG,5*varoutdim,ofield,
     *     rcode)
      call ncapt(cdfout,NCGLOBAL,'nostart',NCLONG,5,nostart,rcode)
      call ncapt(cdfout,NCGLOBAL,'nostop',NCLONG,5,nostop,rcode)
      call ncapt(cdfout,NCGLOBAL,'nodt',NCLONG,5,nodt,rcode)
      call ncapt(cdfout,NCGLOBAL,'ndiag',NCLONG,1,ndiag,rcode)
      call ncapt(cdfout,NCGLOBAL,'diloc',NCLONG,1,diloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'djloc',NCLONG,1,djloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'dkloc',NCLONG,1,dkloc,rcode)
c
      call ncapt(cdfout,NCGLOBAL,'slbtimes',NCLONG,1,slbtimes,rcode)
      call ncapt(cdfout,NCGLOBAL,'nslab',NCLONG,5,nslab,rcode)
      call ncapt(cdfout,NCGLOBAL,'siloc',NCLONG,5,siloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'sjloc',NCLONG,5,sjloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'skloc',NCLONG,5,skloc,rcode)
      call ncapt(cdfout,NCGLOBAL,'slbtype',NCLONG,5,slbtype,rcode)
      call ncapt(cdfout,NCGLOBAL,'noslb',NCLONG,5,noslb,rcode)
      call ncapt(cdfout,NCGLOBAL,'oslab',NCLONG,5*varoutdim,oslab,
     *     rcode)
c
      call ncapt(cdfout,NCGLOBAL,'thtamin',NCFLOAT,1,thtamin,rcode)
      call ncapt(cdfout,NCGLOBAL,'thtamax',NCFLOAT,1,thtamax,rcode)
      call ncapt(cdfout,NCGLOBAL,'salmin',NCFLOAT,1,salmin,rcode)
      call ncapt(cdfout,NCGLOBAL,'salmax',NCFLOAT,1,salmax,rcode)
c
      call ncapt(cdfout,NCGLOBAL,'appden',NCLONG,1,appden,rcode)
      call ncapt(cdfout,NCGLOBAL,'linden',NCLONG,1,linden,rcode)
      call ncapt(cdfout,NCGLOBAL,'potden',NCLONG,1,potden,rcode)
      call ncapt(cdfout,NCGLOBAL,'centhta',NCFLOAT,1,centhta,rcode)
      call ncapt(cdfout,NCGLOBAL,'censal',NCFLOAT,1,censal,rcode)
c
      iii=lnblk(ifile,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'ifile',NCCHAR,
     *            iii,ifile,rcode)
      iii=lnblk(ititle,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'ititle',NCCHAR,
     *            iii,ititle,rcode)
      iii=lnblk(itime,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'itime',NCCHAR,
     *            iii,itime,rcode)
      iii=lnblk(iver,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iver',NCCHAR,
     *            iii,iver,rcode)
      iii=lnblk(com1,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_1',NCCHAR,
     *            iii,com1,rcode)
      iii=lnblk(com2,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_2',NCCHAR,
     *            iii,com2,rcode)
      iii=lnblk(com3,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'iccom_3',NCCHAR,
     *            iii,com3,rcode)
c
      call ncapt(cdfout,NCGLOBAL,'icflgs',NCLONG,20,icflgs,rcode)
      call ncapt(cdfout,NCGLOBAL,'icvals',NCFLOAT,20,icvals,rcode)
c 8-21-92 jrb added tracer loop
      if (ntrace.gt.0) then
         do 295 i=1,ntrace
            write(hold,'(i1)')i
            iii=lnblk(trfile(i),80)
            if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'trcfile_'//hold,
     *            NCCHAR,iii,trfile(i),rcode)
            iii=lnblk(trccom(i),120)
            if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'trccom_'//hold,
     *            NCCHAR,iii,trccom(i),rcode)
 295     continue
      endif

      iii=lnblk(bcfile,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bcfile',NCCHAR,
     *            iii,bcfile,rcode)
      iii=lnblk(bctime,30)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bctime',NCCHAR,
     *            iii,bctime,rcode)
      iii=lnblk(bcver,80)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bcver',NCCHAR,
     *            iii,bcver,rcode)
      iii=lnblk(bccom1,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_1',NCCHAR,
     *            iii,bccom1,rcode)
      iii=lnblk(bccom2,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_2',NCCHAR,
     *            iii,bccom2,rcode)
      iii=lnblk(bccom3,120)
      if(iii.gt.0) call ncaptc(cdfout,NCGLOBAL,'bccom_3',NCCHAR,
     *            iii,bccom3,rcode)
      call ncapt(cdfout,NCGLOBAL,'bcflgs',NCLONG,20,bcflgs,rcode)
      call ncapt(cdfout,NCGLOBAL,'bcvals',NCFLOAT,20,bcvals,rcode)
      call ncapt(cdfout,NCGLOBAL,'bcfrce',NCFLOAT,20,bcfrce,rcode)
    
      call ncendf(cdfout,rcode)
c 
c build the axes
c
      do 300 i=1,idim
        xu(i)=dxx*float(i-1)
        xw(i)=dxx*float(i-1)+dxx*0.5
 300  continue
      xu(idim1)=idim*dxx
      do 310 i=1,kdim
         zu(i)=dz*float(kdim-i) + dz*0.5
         zw(i)=dz*float(kdim-i+1)
 310  continue
      zw(kdim1)=0.0
      do 320 i=1,jtotal
         yv(i)=dy*float(i-1)
         yu(i)=dy*float(i-1) + dy*0.5
 320  continue
      yv(jtotal1)=dy*jtotal
c
c load the axes variable values
c
      call ncvpt(cdfout, xuid, 1, idim1, xu, rcode)
      call ncvpt(cdfout, xwid, 1, idim, xw, rcode)
      call ncvpt(cdfout, yuid, 1, jtotal, yu, rcode)
      call ncvpt(cdfout, yvid, 1, jtotal1, yv, rcode)
      call ncvpt(cdfout, zuid, 1, kdim, zu, rcode)
      call ncvpt(cdfout, zwid, 1, kdim1, zw, rcode)
      call ncvpt(cdfout, xid, 1, 1, 0., rcode)
      call ncvpt(cdfout, yid, 1, 1, 0., rcode)
      call ncvpt(cdfout, zid, 1, 1, 0., rcode)
      call ncvpt(cdfout, t1id, 1, 1, 0, rcode)
c 
C-------Tells ncvpt (for th, sal, rhobm1) where to start writing values.
      strtm(1)=1
      strtm(2)=1
      strtm(3)=1
      strtm(4)=1

C-------Tells ncvpt (for hf, lf, str) how many values to write
      countm(1)=idim
      countm(2)=1
      countm(3)=1
      countm(4)=1

C-------For each mask, write hyperslabs of values into its variable.
      do j=1,jtotal
         strtm(2)=j
         if(hfton)call ncvpt(cdfout,hffid,
     *        strtm,countm,hf_tmask(1,j),rcode)
         if(lfton)call ncvpt(cdfout,lffid,
     *        strtm,countm,lf_tmask(1,j),rcode)
         if(strton)then
            call ncvpt(cdfout,strfid,
     *           strtm,countm,str_tmask(1,j),rcode)
            call ncvpt(cdfout,strtbcfid,
     *           strtm,countm,str_tbc(1,j),rcode)
         endif
      enddo
         
      call ncsnc(cdfout, rcode)
      return
      end
