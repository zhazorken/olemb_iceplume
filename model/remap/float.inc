c
c   $Id: float.inc,v 1.1.1.1 1998/07/20 22:03:48 eds Exp $
c 
      real fltx,flty,fltz
      integer fltstrt,fltstop,fltdt,flttype,nflt
      integer ndesc,fltdescsize,flttimecnt
      integer fltid,fltx_id,flty_id,fltz_id,flttime_id
      integer fltt_id,flts_id
      character fltname*20,fltcomment*80,fltroot*20
c
      common/fltcom/fltx(dscdim,fltdim),flty(dscdim,fltdim),
     *              fltz(dscdim,fltdim),
     *              fltstrt(dscdim),fltstop(dscdim),
     *              fltdt,flttype(dscdim),
     *              nflt,fltid,flttime_id,
     *              fltx_id(dscdim),flty_id(dscdim),fltz_id(dscdim),
     *              fltt_id(dscdim),flts_id(dscdim),
     *              ndesc,fltdescsize(dscdim),flttimecnt
c
      common/fltchr/fltname(dscdim),fltcomment,fltroot(dscdim)
