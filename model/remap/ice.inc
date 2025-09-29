c
c   $Id: ice.inc,v 1.2 1999/08/16 20:22:31 eds Exp $
c
c ice variables: salinity,depth,brine volume, temperature, and
c frazil concentration
c
      real sice(0:idim1,0:jdim1),hice(0:idim1,0:jdim1)
      real vbr(0:idim1,0:jdim1),iceh(0:idim1,0:jdim1)
      real tice(0:idim1,0:jdim1),cfr(0:idim1,0:jdim1,0:kdim1)
#if defined(ICEWALL) || defined(MPOND)
      real wallsf(0:idim1,0:jdim1,0:kdim1)
      real wallhf(0:idim1,0:jdim1,0:kdim1)
      real*8 wsf,whf,wallarea,bottoma,wvst,bvst
#endif
      real icesal(0:idim1,0:jdim1),uice
      real vice
      real*8 cheat(kdim),hicebc,sicebc,sflux(kdim)
c
      common/ice/sice,hice,vbr,tice,cfr,iceh,icesal,uice,
     *           vice
      common/ice8/cheat,hicebc,sicebc,sflux
#if defined(ICEWALL) || defined(MPOND)
      common/icewall/wallsf,wallhf,wsf,whf,wallarea,bottoma,
     *               wvst,bvst
#endif

