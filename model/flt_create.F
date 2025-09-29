      subroutine flt_create(froot)
c
c  $Id: flt_create.F,v 1.8 1999/08/05 15:51:26 dwd Exp $
c
#include "moddef.inc"
c
#include "param.inc"
#include "netcdf.inc"
#include "olemp.inc"
#include "olemi.inc"
#include "float.inc"
#include "refdate.inc"
c
      character ffile*80,froot*80,hold*2
c
      integer  iret,xdim, ydim, zdim, timedim
      integer  xid, yid, zid, i,len,lnblk,tlen,iii
      real xx(fltdim)
      integer  xddim(dscdim),xdid(dscdim),nn,rcode
      integer dims(4), corner(4), edges(4)
      character name*30,tunits*80,str*30
#if defined(ALLOC_DISK)
      integer tpoints
#endif
c data variables
      real zero
      integer eplat,eplon,eppres,eptime
      data zero/0/
      data eplat,eplon,eppres,eptime/500,501,1,625/
c
      len=lnblk(froot,80)
      ffile=froot(1:len)//'_flt.nc'
c
c enter define mode
      fltid = nccre (ffile, NCCLOB, iret)
c
      if(datemode) then
         call frmtdate(refdate, str)
         tunits = 'seconds' // ' since ' // str
         tlen = lnblk(tunits, 80)
      else
         tunits = 'seconds'
         tlen = 7
      endif
c define dimensions
      xdim = ncddef(fltid, 'x', 1, iret)
      ydim = ncddef(fltid, 'y', 1, iret)
      zdim = ncddef(fltid, 'z', 1, iret)
#if defined(ALLOC_DISK)
      tpoints=(ntime-nstart+1)/fltdt+1
      timedim = ncddef(fltid, 'time', tpoints, iret)
#else
      timedim = ncddef(fltid, 'time', NCUNLIM, iret)
#endif
      write(*,*)"Define axis"
c define axis variables
      dims(1) = xdim
      xid = ncvdef (fltid, 'x', NCFLOAT, 1, dims, iret)
      dims(1) = ydim
      yid = ncvdef (fltid, 'y', NCFLOAT, 1, dims, iret)
      dims(1) = zdim
      zid = ncvdef (fltid, 'z', NCFLOAT, 1, dims, iret)
      dims(1) = timedim
      flttime_id = ncvdef (fltid, 'time', NCLONG, 1, dims, iret)
      write(*,*)"Define axis attrib"
c assign axis attributes
      call ncapt(fltid, xid, 'epic_code', NCLONG, 1, EPLON, iret)
      call ncaptc(fltid, xid, 'FORTRAN_format', NCCHAR, 4, 'f6.2', iret)
      call ncaptc(fltid, xid, 'units', NCCHAR, 1, 'm', iret)
      call ncaptc(fltid, xid, 'type', NCCHAR, 4, 'EVEN', iret)
c
      call ncapt(fltid, yid, 'epic_code', NCLONG, 1, EPLAT, iret)
      call ncaptc(fltid, yid, 'FORTRAN_format', NCCHAR, 4, 'f6.2', iret)
      call ncaptc(fltid, yid, 'units', NCCHAR, 1, 'm', iret)
      call ncaptc(fltid, yid, 'type', NCCHAR, 4, 'EVEN', iret)
c
      call ncapt(fltid, zid, 'epic_code', NCLONG, 1, EPPRES, iret)
      call ncaptc(fltid, zid, 'FORTRAN_format', NCCHAR, 4, 'f6.2', iret)
      call ncaptc(fltid, zid, 'units', NCCHAR, 1, 'm', iret)
      call ncaptc(fltid, zid, 'type', NCCHAR, 4, 'EVEN', iret)
c
      write(*,*)"Define flttime"
      if(.not.datemode) then
         call ncapt(fltid, flttime_id, 'epic_code', NCLONG, 1, 
     *        EPTIME, iret)
      endif
      call ncaptc(fltid, flttime_id,'FORTRAN_format',NCCHAR, 3, 
     *     'i10', iret)
         call ncaptc(fltid, flttime_id, 'units', NCCHAR, tlen, 
     *        tunits, iret)
      call ncaptc(fltid, flttime_id, 'type', NCCHAR, 1, 'EVEN', iret)
c
c define special float number axis
c
      do 100 i = 1, ndesc
         len=lnblk(fltname(i),20)
         name=fltname(i)(1:len)//'x'
         iret = NF_DEF_DIM(fltid,name,fltdescsize(i),xddim(i))
         write(*,*)"ncddef return ",iret
         xdid(i) =ncvdef(fltid,name,NCFLOAT,1,xddim(i),iret)
         write(*,*)"ncvdef return",iret
         call ncapt(fltid, xdid(i), 'epic_code', NCLONG, 1, EPLON, iret)
         call ncaptc(fltid, xdid(i), 'FORTRAN_format', NCCHAR,
     *     4,'f6.2', iret)
         call ncaptc(fltid, xdid(i), 'units', NCCHAR, 1, 'm', iret)
         call ncaptc(fltid, xdid(i), 'type' , NCCHAR, 4, 'EVEN', iret)
         do 110 nn=1,fltdescsize(i)
            xx(nn)=float(nn)
110      continue
100   continue
c
      write(*,*)"Define global attrib"
c global attributes
      len=lnblk(ident,30)
      if(len.gt.0)call ncaptc(fltid, NCGLOBAL, 'identifier', 
     *     NCCHAR, len, ident, iret)
      len=lnblk(mver,80)
      if(len.gt.0)call ncaptc(fltid, NCGLOBAL, 'model_version', 
     *     NCCHAR, len, mver, iret)

      do i=1,4
         iii=lnblk(modflags(i),80)
         if(iii.gt.0) then
            write(hold,'(i1)')i
            call ncaptc(fltid, NCGLOBAL,'model_flags'//hold,
     *           NCCHAR, iii, modflags(i), iret)
         endif
      enddo

      len=lnblk(title,80)
      if(len.gt.0)call ncaptc(fltid, NCGLOBAL, 'title', NCCHAR, 
     *     len, title, iret)
      len=lnblk(fltcomment,80)
      if(len.gt.0)call ncaptc(fltid, NCGLOBAL, 'flt_comment', NCCHAR, 
     *     len, fltcomment, iret)
c
      call ncaptc(fltid, NCGLOBAL, 'DATA_TYPE', NCCHAR, 5, 
     *     'MODEL', iret)
      call ncaptc(fltid, NCGLOBAL, 'COORD_SYSTEM', NCCHAR, 5, 
     *     'LOCAL', iret)
c define variables
      write(*,*)"Define Variables"
      dims(4) = timedim
      dims(3) = zdim
      dims(2) = ydim
      do 10 i=1,ndesc
         dims(1) = xddim(i)
         len=lnblk(fltname(i),20)
         name=fltname(i)(1:len)//'_x'
         fltx_id(i) = ncvdef (fltid, name, NCFLOAT, 4, dims, iret)
         name=fltname(i)(1:len)//'_y'
         flty_id(i) = ncvdef (fltid, name, NCFLOAT, 4, dims, iret)
         name=fltname(i)(1:len)//'_z'
         fltz_id(i) = ncvdef (fltid, name, NCFLOAT, 4, dims, iret)
         name=fltname(i)(1:len)//'_t'
         fltt_id(i) = ncvdef (fltid, name, NCFLOAT, 4, dims, iret)
         name=fltname(i)(1:len)//'_s'
         flts_id(i) = ncvdef (fltid, name, NCFLOAT, 4, dims, iret)
c assign attributes
         call ncaptc(fltid, fltx_id(i), 'units', NCCHAR, 1, 'm', iret)
         call ncaptc(fltid, fltx_id(i), 'long_name', NCCHAR, 18, 
     *        'east-west position', iret)
         call ncaptc(fltid, fltx_id(i), 'FORTRAN_format', NCCHAR, 4, 
     *        'f8.2', iret)
         call ncapt(fltid, fltx_id(i), '_FillValue', NCFLOAT, 1, 
     *        1e35, iret)
c     
         call ncaptc(fltid, flty_id(i), 'units', NCCHAR, 1, 'm', iret)
         call ncaptc(fltid, flty_id(i), 'long_name', NCCHAR, 20, 
     *        'north-south position', iret)
         call ncaptc(fltid, flty_id(i), 'FORTRAN_format', NCCHAR, 4, 
     *        'f8.2', iret)
         call ncapt(fltid, flty_id(i), '_FillValue', NCFLOAT, 1, 
     *        1e35, iret)
c     
         call ncaptc(fltid, fltz_id(i), 'units', NCCHAR, 1, 'm', iret)
         call ncaptc(fltid, fltz_id(i), 'long_name', NCCHAR,5, 
     *        'depth', iret)
         call ncaptc(fltid, fltz_id(i), 'FORTRAN_format', NCCHAR, 4, 
     *        'f8.2', iret)
         call ncapt(fltid, fltz_id(i), '_FillValue', NCFLOAT, 1, 
     *        1e35, iret)
c     
         call ncaptc(fltid, fltt_id(i), 'units', NCCHAR, 4, 
     *        'degC', iret)
         call ncaptc(fltid, fltt_id(i), 'long_name', NCCHAR,11, 
     *        'Temperature', iret)
         call ncaptc(fltid, fltt_id(i), 'FORTRAN_format', NCCHAR, 4, 
     *        'f8.2', iret)
         call ncapt(fltid, fltt_id(i), '_FillValue', NCFLOAT, 1, 
     *        1e35, iret)
c     
         call ncaptc(fltid, flts_id(i), 'units', NCCHAR, 3, 'PSU', iret)
         call ncaptc(fltid, flts_id(i), 'long_name', NCCHAR,8, 
     *        'Salinity', iret)
         call ncaptc(fltid, flts_id(i), 'FORTRAN_format', NCCHAR, 4, 
     *        'f8.2', iret)
         call ncapt(fltid, flts_id(i), '_FillValue', NCFLOAT, 1, 
     *        1e35, iret)
 10   continue
      do 30 i=1,ndesc
c assign global variable attributes
         len=lnblk(fltroot(i),20)
         name=fltroot(i)(1:len)//'_points'
         call ncapt(fltid,NCGLOBAL,name, NCLONG, 1, 
     *        fltdescsize(i), iret)
 30   continue
c leave define mode
      write(*,*)"Leave define mode"
      call ncendf(fltid, iret)
c store x
      corner(1) = 1
      edges(1) = 1
      call ncvpt(fltid, xid, corner, edges, zero, iret)
c store y
      call ncvpt(fltid, yid, corner, edges, zero, iret)
c store z
      call ncvpt(fltid, zid, corner, edges, zero, iret)
c store float number
      write(*,*)"Store float number"
      do 120 i=1,ndesc
         call ncvpt(fltid, xdid(i), 1, fltdescsize(i), xx, iret)
120   continue
c
      flttimecnt = 0
      write(*,*)"Sync file"
      call ncsnc(fltid, rcode)    
      write(*,*)"Exit"
c
      return
      end
