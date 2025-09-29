c
c   $Id: olemtke.inc,v 1.1.1.1 1998/07/20 22:03:49 eds Exp $
c
      real kmom(0:idim1,0:jdim1,0:kdim1),kheat(0:idim1,0:jdim1,0:kdim1)
      real kmomp(0:idim1,0:jdim1,0:kdim1)
      real kmb,khb,kms,nu0
      logical usesgs,usedeardorff,advtke
c
      common/tkediff/kmom,kmomp,kheat,kmb,khb,kms,usesgs,
     *     usedeardorff,nu0,advtke
