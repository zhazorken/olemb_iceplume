c     
c	$Id: ic.inc.in,v 4.3 1997/11/17 18:24:58 dwd Exp $
c     
      real th(idim,kdim,jdim)
#if defined(SNGLSALT)
      real sal(idim,kdim,jdim)
#else
      real*8 sal(idim,kdim,jdim)
#endif
      real rhobm1(kdim)
      real salbar(kdim)
      real thbar(kdim)
c      real u(idim,kdim,jdim),v(idim,kdim,jdim)
c      real trc(idim,kdim,jdim,trcdim)
      logical hftst,lftst,strtst
      real hf_tmask(idim,jdim),lf_tmask(idim,jdim)
      real str_tmask(idim,jdim)
      integer str_tbc(idim,jdim)
      integer imax,jmax,kmax,imax1,jmax1
c      
      common/comt/ th,sal,rhobm1,hf_tmask,lf_tmask,
     *     str_tmask,str_tbc,thbar,salbar
      common/comid/  hftst,lftst,strtst
      common/commax/imax,jmax,kmax,imax1,jmax1
