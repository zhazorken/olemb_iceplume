c
c   $Id: finitv.inc,v 1.2 1999/08/16 20:22:31 eds Exp $
c
c
c volume and cell boundary areas
c
      real ax(0:idim1,0:jdim1,0:kdim1)
      real ay(0:idim1,0:jdim1,0:kdim1)
      real bay(0:idim1,0:jdim1,0:kdim1)
      real vay(0:idim1,0:jdim1,0:kdim1)
      real uax(0:idim1,0:jdim1,0:kdim1)
      real bax(0:idim1,0:jdim1,0:kdim1)
      real way(0:idim1,0:jdim1,0:kdim1)
      real dhdy(0:idim1,0:jdim1)
      real dhdx(0:idim1,0:jdim1)
      real calpha(0:idim1,0:jdim1)
      real salpha(0:idim1,0:jdim1)
      integer kdist(0:idim1,0:jdim1)
      real frac(0:idim1,0:jdim1)
      real harea(0:idim1,0:jdim1)

      common/finitv/ax,ay,way,bay,bax,uax,vay,dhdx,dhdy,frac,salpha,
     *       calpha,harea,kdist
