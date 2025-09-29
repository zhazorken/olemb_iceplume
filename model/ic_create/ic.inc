c     
c	$Id: ic.inc.in,v 4.3 1997/11/17 18:24:58 dwd Exp $
c     
      real th(0:idim1,0:kdim1,0:jdim1)
#if defined(SNGLSALT)
      real sal(0:idim1,0:kdim1,0:jdim1)
#else
      real*8 sal(0:idim1,0:kdim1,0:jdim1)
#endif
      real rhobm1(0:kdim1)
      real salbar(0:kdim1)
      real thbar(0:kdim1)
c      real u(0:idim1,0:kdim1,0:jdim1),v(0:idim1,0:kdim1,0:jdim1)
c      real trc(0:idim1,0:kdim1,0:jdim1,trcdim)
      logical hftst,lftst,strtst
      real hf_tmask(0:idim1,0:jdim1),lf_tmask(0:idim1,0:jdim1)
      real str_tmask(0:idim1,0:jdim1)
      integer str_tbc(0:idim1,0:jdim1)
      integer imax,jmax,kmax,imax1,jmax1
c      
      common/comt/ th,sal,rhobm1,hf_tmask,lf_tmask,
     *     str_tmask,str_tbc,thbar,salbar
      common/comid/  hftst,lftst,strtst
      common/commax/imax,jmax,kmax,imax1,jmax1
