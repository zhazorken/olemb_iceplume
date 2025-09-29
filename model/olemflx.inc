c
c   $Id: olemflx.inc,v 1.1.1.1 1998/07/20 22:03:49 eds Exp $
c
      real uaa(0:idim1,0:jdim1,0:kdim1),vaa(0:idim1,0:jdim1,0:kdim1)
      real waa(0:idim1,0:jdim1,0:kdim1)
      real taa(0:idim1,0:jdim1,0:kdim1)
c
      real saa(0:idim1,0:jdim1,0:kdim1)
c

c
      real uwa(0:idim1,0:jdim1,0:kdim1),uva(0:idim1,0:jdim1,0:kdim1)
      real vwa(0:idim1,0:jdim1,0:kdim1)
c
      real uta(0:idim1,0:jdim1,0:kdim1),wta(0:idim1,0:jdim1,0:kdim1)
      real vta(0:idim1,0:jdim1,0:kdim1)
c
      real usa(0:idim1,0:jdim1,0:kdim1),wsa(0:idim1,0:jdim1,0:kdim1)
      real vsa(0:idim1,0:jdim1,0:kdim1)
c
      common/aveflux/uaa,vaa,waa,uwa,uva,vwa,uta,vta,wta,usa,vsa,wsa,
     *       taa,saa
