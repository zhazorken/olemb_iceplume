      subroutine slb_create(ifile, bcfile, froot, ko)
c
c   $Id: slb_create.F,v 1.8 1999/08/16 20:22:32 eds Exp $
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
      character ifile*80,froot*80,temp*80,name*15,namx*15,namy*15
      character namz*15,hold*2,bcfile*80,str*30,tunits*80
      integer ko,cdfout
      integer lnblk,iii,j,ilen,i,k,tlen,index
      integer rcode,varidx(4),varidux(4),varidvx(4)
      integer varidy(4),variduy(4),varidvy(4)
      integer varidwx(4),varidwy(4)
      integer variduz(4),varidvz(4),varidz(4)
      integer varidwz(4)
      integer xudim,xwdim,zudim,zwdim,tdim
      integer xupdim,yupdim,xwpdim,yvpdim
      integer zpdim,zwpdim
      integer xupid,yupid,xwpid,yvpid,zupid,zwpid
      integer xuid,xwid,zuid,zwid
      integer yudim,yvdim,yuid,yvid
      integer t1dim,t1id
      real xu(itotal),xw(itotal),zu(kdim),zw(kdim1)
      real yu(jtotal),yv(jtotal)
      real xup,yup,xwp,yvp,zup,zwp,dxt(itotal)
      integer type
#if defined(ALLOC_DISK)
      integer tstart,tstop,tpoints
#else
      integer omode
#endif
c
      integer eplat,eplon,eppres,eptime
      data eplat,eplon,eppres,eptime/500,501,1,625/
c
c for now just do the first output list
c
C-----create slab cdf file (CDF)
c
      ilen=lnblk(froot,80)
      if(ko.eq.1)then
         slbid(ko)=nccre(froot(1:ilen)//'_slb.nc',ncclob,rcode)
      else
         write(temp,999)froot(1:ilen),ko
 999     format(a,'_',i1)
         ilen=lnblk(temp,80)
         slbid(ko)=nccre(temp(1:ilen)//'_slb.nc',ncclob,rcode)
      endif
#if !defined(ALLOC_DISK)
c
c don't use fills
c
      omode=ncsfil(slbid(ko),NCNOFILL,rcode)
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
      cdfout=slbid(ko)
      type=slbtype(ko)
c
      xudim=ncddef(cdfout,'xu',itotal,rcode)
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
      xwdim=ncddef(cdfout,'xw',itotal,rcode)
      xwid  = ncvdef(cdfout, 'xw', NCFLOAT, 1, xwdim, rcode)
      call ncapt(cdfout, xwid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xwid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xwid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xwid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yvdim=ncddef(cdfout,'yv',jtotal,rcode)
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
      t1dim=ncddef(cdfout,'t1',1,rcode)
      t1id = ncvdef(cdfout, 't1', NCLONG, 1, t1dim, rcode)
      if(.not.datemode)then
         call ncapt(cdfout, t1id, 'epic_code', NCLONG, 1, 
     *        EPTIME, rcode)
      endif
      call ncaptc(cdfout, t1id, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, t1id, 'units', NCCHAR, tlen, tunits,rcode)
      call ncaptc(cdfout, t1id, 'type' , NCCHAR, 1, ' ', rcode)
c
      xupdim=ncddef(cdfout,'xup',1,rcode)
      xupid =ncvdef(cdfout,'xup',NCFLOAT,1,xupdim,rcode)
      call ncapt(cdfout, xupid, 'epic_code', NCLONG, 1, EPLON, rcode)
      call ncaptc(cdfout, xupid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xupid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xupid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yupdim=ncddef(cdfout,'yup',1,rcode)
      yupid  = ncvdef(cdfout, 'yup', NCFLOAT, 1, yupdim, rcode)
      call ncapt(cdfout, yupid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yupid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yupid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yupid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      xwpdim=ncddef(cdfout,'xwp',1,rcode)
      xwpid  = ncvdef(cdfout, 'xwp', NCFLOAT, 1, xwpdim, rcode)
      call ncapt(cdfout, xwpid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, xwpid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, xwpid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, xwpid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      yvpdim=ncddef(cdfout,'yvp',1,rcode)
      yvpid  = ncvdef(cdfout, 'yvp', NCFLOAT, 1, yvpdim, rcode)
      call ncapt(cdfout, yvpid, 'epic_code', NCLONG, 1, EPLAT, rcode)
      call ncaptc(cdfout, yvpid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, yvpid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, yvpid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zpdim=ncddef(cdfout,'zup',1,rcode)
      zupid =ncvdef(cdfout,'zup',NCFLOAT,1,zpdim,rcode)
      call ncapt(cdfout, zupid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zupid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zupid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zupid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
      zwpdim=ncddef(cdfout,'zwp',1,rcode)
      zwpid =ncvdef(cdfout,'zwp',NCFLOAT,1,zwpdim,rcode)
      call ncapt(cdfout, zwpid, 'epic_code', NCLONG, 1, EPPRES, rcode)
      call ncaptc(cdfout, zwpid, 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', rcode)
      call ncaptc(cdfout, zwpid, 'units', NCCHAR, 1, 'm', rcode)
      call ncaptc(cdfout, zwpid, 'type' , NCCHAR, 4, 'EVEN', rcode)
c
#if defined(ALLOC_DISK)
      tpoints=(ntime-nstart+1)/nslab(ko)+1
      tdim=ncddef(cdfout,'time',tpoints,rcode)
#else
      tdim=ncddef(cdfout,'time',NCUNLIM,rcode)
#endif
      tsid(ko) = ncvdef(cdfout, 'time', NCLONG, 1, tdim, rcode)
      if(.not.datemode) then
         call ncapt(cdfout, tsid, 'epic_code', NCLONG, 1, 
     *        EPTIME, rcode)
      endif
      call ncaptc(cdfout, tsid, 'FORTRAN_format', NCCHAR,
     *     3,'i10', rcode)
      call ncaptc(cdfout, tsid, 'units', NCCHAR, tlen, tunits,rcode)
      call ncaptc(cdfout, tsid, 'type' , NCCHAR, 1, ' ', rcode)
c
C-----define cell-centered variables (CDF)
      varidx(1)=xwdim
      varidx(2)=yupdim
      varidx(3)=zudim
      varidx(4)=tdim
c
      varidy(1)=xwpdim
      varidy(2)=yudim
      varidy(3)=zudim
      varidy(4)=tdim
c
c
C-----define u variable (x-z plot) (shifted to left of center) (CDF)
      varidux(1)=xudim
      varidux(2)=yupdim
      varidux(3)=zudim
      varidux(4)=tdim
c
C-----define u variable (y-z plot) (shifted to left of center) (CDF)
c
      variduy(1)=xupdim
      variduy(2)=yudim
      variduy(3)=zudim
      variduy(4)=tdim

C-----define v variable (x-z plot)
      varidvx(1)=xwdim
      varidvx(2)=yvpdim
      varidvx(3)=zudim
      varidvx(4)=tdim
c
C-----define v variable (y-z plot)
      varidvy(1)=xwpdim
      varidvy(2)=yvdim
      varidvy(3)=zudim
      varidvy(4)=tdim

C-----define w variable (x-z plot) (shifted to bottom of center) (CDF)
      varidwx(1)=xwdim
      varidwx(2)=yupdim
      varidwx(3)=zwdim
      varidwx(4)=tdim
c
C-----define w variable (y-z plot) (shifted to bottom of center) (CDF)
      varidwy(1)=xwpdim
      varidwy(2)=yudim
      varidwy(3)=zwdim
      varidwy(4)=tdim
c
      variduz(1)=xudim
      variduz(2)=yudim
      variduz(3)=zpdim
      variduz(4)=tdim
c 
      varidvz(1)=xwdim 
      varidvz(2)=yvdim 
      varidvz(3)=zpdim
      varidvz(4)=tdim
c 
      varidwz(1)=xwdim 
      varidwz(2)=yudim 
      varidwz(3)=zwpdim
      varidwz(4)=tdim
c
c 
      varidz(1)=xwdim 
      varidz(2)=yudim 
      varidz(3)=zpdim
      varidz(4)=tdim

C---- set attributes
      do 980 j=1,noslb(ko)
         call varnme(oslab(ko,j),name)
         ilen=lnblk(name,15)
c
c define names to agree with olem90
c
         namx = name(:ilen)//'_y'
         namy = name(:ilen)//'_x'
         namz = name(:ilen)//'_z'
         ilen=ilen+2
         if(oslab(ko,j).eq.MDHICE) then
c
c for sea ice (or other surface array) only define the vertical xsect
c
           hifidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,hifidz(ko))
         else if(oslab(ko,j).eq.MDFICE) then
           hfifidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,hfifidz(ko))
         else if(oslab(ko,j).eq.MDSFICE) then
           sfifidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,sfifidz(ko))
         else if(oslab(ko,j).eq.MDUF) then
           uffidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,uffidz(ko))
         else if(oslab(ko,j).eq.MDVF) then
           vffidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,vffidz(ko))
         else if(oslab(ko,j).eq.MDWF) then
           wffidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,wffidz(ko))
         else if(oslab(ko,j).eq.MDRHOF) then
           rffidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,rffidz(ko))
c
c
#if defined(HBOT)
         else if(oslab(ko,j).eq.MDHBOT) then
c
c for bottom (or other surface array) only define the vertical xsect
c
           botfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
           call setvar(oslab(ko,j),cdfout,botfidz(ko))
c
c
#endif

         else if(oslab(ko,j).eq.MDUVEL) then
            if(type.eq.SLAB_ALL)then
c
               ufidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidux,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,ufidx(ko))
c           
               ufidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,variduy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,ufidy(ko))
c
            endif
c
            ufidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,variduz,rcode)
            call setvar(oslab(ko,j),cdfout,ufidz(ko))
c
         elseif(oslab(ko,j).eq.MDVVEL) then
            if(type.eq.SLAB_ALL)then
c
               vfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidvx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,vfidx(ko))
c
               vfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidvy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,vfidy(ko))
c
            endif
c
            vfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidvz,rcode)
            call setvar(oslab(ko,j),cdfout,vfidz(ko))
c
         elseif(oslab(ko,j).eq.MDWVEL) then
            if(type.eq.SLAB_ALL)then
c
               wfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidwx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,wfidx(ko))
c
               wfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidwy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,wfidy(ko))
c
            endif
c
            wfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidwz,rcode)
            call setvar(oslab(ko,j),cdfout,wfidz(ko))
c
         elseif(oslab(ko,j).eq.MDVADV) then
            if(type.eq.SLAB_ALL)then
c
               vafidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,vafidx(ko))
c
               vafidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,vafidy(ko))
c           
            endif
c           
            vafidz(ko)=
     *         ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,vafidz(ko))
c        
         elseif(oslab(ko,j).eq.MDDCOR) then
            if(type.eq.SLAB_ALL)then
c
               dcfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,dcfidx(ko))
c
               dcfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,dcfidy(ko))
c           
            endif
c           
            dcfidz(ko)=
     *         ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,dcfidz(ko))
c        
         elseif(oslab(ko,j).eq.MDDPRES) then
            if(type.eq.SLAB_ALL)then
c
               dpfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,dpfidx(ko))
c
               dpfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,dpfidy(ko))
c           
            endif
c           
            dpfidz(ko)=
     *        ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,dpfidz(ko))
c        
c        
         elseif(oslab(ko,j).eq.MDSGF) then
            if(type.eq.SLAB_ALL)then
c
               sgfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,sgfidx(ko))
c
               sgfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,sgfidy(ko))
c           
            endif
c           
            sgfidz(ko)=
     *        ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,sgfidz(ko))
c        


         elseif(oslab(ko,j).eq.MDPRES) then
            if(type.eq.SLAB_ALL)then
c
               pfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,pfidx(ko))
c
               pfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,pfidy(ko))
c
            endif
c
            pfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,pfidz(ko))
c
         elseif(oslab(ko,j).eq.MDTEMP) then
            if(type.eq.SLAB_ALL)then
c
               thfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,thfidx(ko))
c
               thfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,thfidy(ko))
c
            endif
c
            thfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,thfidz(ko))
c
         elseif(oslab(ko,j).eq.MDSAL) then
            if(type.eq.SLAB_ALL)then
c
               salfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,salfidx(ko))
c
               salfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,salfidy(ko))
c 
            endif
c
            salfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,
     *           rcode)
            call setvar(oslab(ko,j),cdfout,salfidz(ko))
         elseif(oslab(ko,j).eq.MDTKE) then
            if(type.eq.SLAB_ALL)then
c
               tkefidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,tkefidx(ko))
c
               tkefidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,tkefidy(ko))
c
            endif
c
            tkefidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,
     *           rcode)
            call setvar(oslab(ko,j),cdfout,tkefidz(ko))
c
         elseif(oslab(ko,j).eq.MDRHO) then
            if(type.eq.SLAB_ALL)then
c
               rho_pfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,rho_pfidx(ko))
c
               rho_pfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,rho_pfidy(ko))
c
            endif
c     
            rho_pfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,
     *           rcode)
            call setvar(oslab(ko,j),cdfout,rho_pfidz(ko))
         else if(oslab(ko,j).eq.MDKH)then
            if(type.eq.SLAB_ALL)then
c
               khfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,khfidx(ko))
c
               khfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,khfidy(ko))
c
            endif
c
            khfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,khfidz(ko))
         else if(oslab(ko,j).eq.MDKM)then
            if(type.eq.SLAB_ALL)then
c
               kmfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,kmfidx(ko))
c
               kmfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,kmfidy(ko))
c
            endif
c
            kmfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,rcode)
            call setvar(oslab(ko,j),cdfout,kmfidz(ko))
         else if(oslab(ko,j).eq.MDRHOT)then
            if(type.eq.SLAB_ALL)then
c
               rhofidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,varidx,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,rhofidx(ko))
c
               rhofidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,varidy,
     *              rcode)
               call setvar(oslab(ko,j),cdfout,rhofidy(ko))
c
            endif
c
            rhofidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,
     *           rcode)
            call setvar(oslab(ko,j),cdfout,rhofidz(ko))
         else if(oslab(ko,j).eq.MDHFSGS)then
            if(type.eq.SLAB_ALL)then
c
               hf_sgsfidx(ko)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,
     *              varidx,rcode)
               call setvar(oslab(ko,j),cdfout,hf_sgsfidx(ko))
c
               hf_sgsfidy(ko)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,
     *              varidy,rcode)
               call setvar(oslab(ko,j),cdfout,hf_sgsfidy(ko))
c
            endif
c
            hf_sgsfidz(ko)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,
     *           varidz,rcode)
            call setvar(oslab(ko,j),cdfout,hf_sgsfidz(ko))
         else if(oslab(ko,j).eq.MDTRCR)then
            k=oslab(ko,j)-MDTRCR
            k=1
            write(str,994)k
 994        format('trcr',i1)
            if(type.eq.SLAB_ALL)then
               trcrfidx(ko,k)=ncvdef(cdfout,namx(:ilen),NCFLOAT,4,
     *              varidx,rcode)
               call setvar(oslab(ko,j),cdfout,trcrfidx(ko,k))
c
               trcrfidy(ko,k)=ncvdef(cdfout,namy(:ilen),NCFLOAT,4,
     *              varidy,rcode)
               call setvar(oslab(ko,j),cdfout,trcrfidy(ko,k))
c
            endif
c
            trcrfidz(ko,k)=ncvdef(cdfout,namz(:ilen),NCFLOAT,4,varidz,
     *           rcode)
            call setvar(oslab(ko,j),cdfout,trcrfidz(ko,k))
         endif
 980  continue
c
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
            call ncaptc(cdfout, NCGLOBAL,'model_flags'//hold,
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
      call ncapt(cdfout,NCGLOBAL,'dx',NCFLOAT,1,dx,rcode)
      call ncapt(cdfout,NCGLOBAL,'dy',NCFLOAT,1,dy,rcode)
      call ncapt(cdfout,NCGLOBAL,'dz',NCFLOAT,1,dz,rcode)
      call ncapt(cdfout,NCGLOBAL,'imax',NCLONG,1,itotal,rcode)
      call ncapt(cdfout,NCGLOBAL,'jmax',NCLONG,1,jtotal,rcode)
      call ncapt(cdfout,NCGLOBAL,'kmax',NCLONG,1,kdim,rcode)
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
      call ncapt(cdfout,NCGLOBAL,'slab',NCLONG,1,ko,rcode)
      call ncapt(cdfout,NCGLOBAL,'nslab',NCLONG,1,nslab(ko),rcode)
      if(siloc(ko) .eq. -1) then
      call ncapt(cdfout,NCGLOBAL,'sxloc',NCFLOAT,1,sxloc(ko),rcode)
      else
      call ncapt(cdfout,NCGLOBAL,'siloc',NCLONG,1,siloc(ko),rcode)
      endif
      if(sjloc(ko) .eq. -1) then
      call ncapt(cdfout,NCGLOBAL,'syloc',NCFLOAT,1,syloc(ko),rcode)
      else
      call ncapt(cdfout,NCGLOBAL,'sjloc',NCLONG,1,sjloc(ko),rcode)
      endif
      if(skloc(ko) .eq. -1) then
      call ncapt(cdfout,NCGLOBAL,'szloc',NCFLOAT,1,szloc(ko),rcode)
      else
      call ncapt(cdfout,NCGLOBAL,'skloc',NCLONG,1,skloc(ko),rcode)
      endif
      call ncapt(cdfout,NCGLOBAL,'slbtype',NCLONG,1,slbtype(ko),rcode)
      call ncapt(cdfout,NCGLOBAL,'noslb',NCLONG,1,noslb(ko),rcode)
      call ncapt(cdfout,NCGLOBAL,'oslab',NCLONG,5*varoutdim,
     *     oslab,rcode)
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
c MJH TODO: BEGIN take this out - moved to main
c      if(siloc(ko) .eq. -1) then
c         index = int(sxloc(ko)/dx + 0.5)
c         if(index .lt. 1) index = 1
c         if(index .gt. idim) index = idim
c         siloc(ko) = index
c      endif
c MJH TODO: END take this out - moved to main
      xup=dxx*float(siloc(ko)-1)
      xwp=xup+0.5*dxx
c
c MJH TODO: BEGIN take this out - moved to main
c      if(sjloc(ko) .eq. -1) then
c         index = int(syloc(ko)/dy + 0.5)
c         if(index .lt. 1) index = 1
c         if(index .gt. jtotal) index = idim
c         sjloc(ko) = index
c      endif
c MJH TODO: END take this out - moved to main
      yvp=dy*float(sjloc(ko)-1)
      yup=yvp+0.5*dy
c
c MJH TODO: BEGIN take this out - moved to main
c      if(skloc(ko) .eq. -1) then
c         index = kdim - int(szloc(ko)/dz + 0.5)
c         if(index .lt. 1) index = 1
c         if(index .gt. kdim) index = kdim
c         skloc(ko) = index
c      endif
c MJH TODO: END take this out - moved to main
      zup=dz*float(kdim-skloc(ko)) + dz*0.5
      zwp=dz*float(kdim-skloc(ko)+1)
#if defined(RADTOP)
      do 300 i=1,itotal
        xu(i)=dx*float(i-1)
        xw(i)=dx*float(i-1)+dx*0.5
 300  continue
      do 310 i=1,kdim
         zu(i)=dz*float(i-1) + dz*0.5
         zw(i)=dz*float(i)
 310  continue
      zw(kdim1)=dz*kdim1
      do 320 i=1,jtotal
         yv(i)=dy*float(i-1)
         yu(i)=dy*float(i-1) + dy*0.5
 320  continue
#else

      xu(1) = 0.0
      xw(1) = 0.5*dxx
      do i=1,itotal
        if(i.gt.itotal/2) then
          dxt(i) = dxx+6*dxx*sin(3.14159*(i-itotal/2)/float(itotal))
        else
          dxt(i) = dxx
        endif
      enddo 
      do i=2,itotal
        xu(i)=xu(i-1)+dxt(i)
        xw(i)=xw(i-1)+0.5*(dxt(i)+dxt(i-1))
      enddo
 300  continue
      do 310 i=1,kdim
         zu(i)=dz*float(kdim-i) + dz*0.5
         zw(i)=dz*float(kdim-i+1)
 310  continue
      zw(kdim1)=0.0
      do 320 i=1,jtotal
         yv(i)=dy*float(i-1)
         yu(i)=dy*float(i-1) + dy*0.5
 320  continue
#endif
c
c load the axes variable values
c
      call ncvpt(cdfout, xuid, 1, itotal, xu, rcode)
      call ncvpt(cdfout, xwid, 1, itotal, xw, rcode)
      call ncvpt(cdfout, yuid, 1, jtotal, yu, rcode)
      call ncvpt(cdfout, yvid, 1, jtotal, yv, rcode)
      call ncvpt(cdfout, zuid, 1, kdim, zu, rcode)
      call ncvpt(cdfout, zwid, 1, kdim1, zw, rcode)
      call ncvpt(cdfout, xupid, 1, 1, xup, rcode)
      call ncvpt(cdfout, xwpid, 1, 1, xwp, rcode)
      call ncvpt(cdfout, yupid, 1, 1, yup,  rcode)
      call ncvpt(cdfout, yvpid, 1, 1, yvp, rcode)
      call ncvpt(cdfout, zupid, 1, 1, zup, rcode)
      call ncvpt(cdfout, zwpid, 1, 1, zwp, rcode)
      call ncvpt(cdfout, t1id, 1, 1, 0, rcode)
c 
      call ncsnc(cdfout, rcode)
      return
      end
