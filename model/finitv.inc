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
      real bway(0:idim1,0:jdim1,0:kdim1)
      real dhdy(0:idim1,0:jdim1)
      real dhdx(0:idim1,0:jdim1)
      real calpha(0:idim1,0:jdim1)
      real salpha(0:idim1,0:jdim1)
      integer kdist(0:idim1,0:jdim1),iwall(jdim,kdim)
      real frac(0:idim1,0:jdim1)
      real harea(0:idim1,0:jdim1)
      real mslope,yinter
      real hx(1:kdim*4+4)
      real hdx(0:idim1,0:jdim1,0:kdim1)
      real disch,plumw,sedinit,dphip

      common/finitv/ax,ay,bway,way,bay,bax,uax,vay,dhdx,dhdy,frac,salpha,
     *       calpha,harea,kdist,iwall,mslope,yinter,hx,hdx,disch,dphip,
     *       plumw,sedinit
