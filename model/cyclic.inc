c
c   $Id: cyclic.inc,v 1.1.1.1 1998/07/20 22:03:47 eds Exp $
c
c     cyclic common block
c
      real randn(0:idim1,0:jdim1),ranamp,heatf,latent,ustr,vstr
      real lheat(nsize)
      real ust(nsize),vst(nsize),us(0:kdim1),swheat(nsize)
      real vs(0:kdim1),heat(nsize),wavel(nsize),waveh(nsize)
      real angle(nsize),swheatf,rainf,rain(nsize)
      real rk(idim,jdim)
      real rad(kdim)

      real bvf,dsdzlow,dtdzlow,ustk,vstk,wavelen,waveamp,waveang
      real bvft,dsdztop,dtdztop
      real hf_tmask(0:idim1,0:jdim1),lf_tmask(0:idim1,0:jdim1)
      real str_tmask(0:idim1,0:jdim1)
      integer tt(nsize),ncool,nf,nrndom,timsiz,nstage
      integer str_tbc(0:idim1,0:jdim1)
      logical bcon,heaton,uston,vston,wavelon,wavehon,angleon
      logical hfton,lfton,strton,swheaton,rainon,lheaton
c
      common/radia/ randn,ncool,nf,ranamp,nrndom,heatf,latent,swheatf,
     *   tt,ust,vst,wavel,waveh,angle,us,vs,heat,timsiz,bcon,swheat,
     *   bvf,rk,nstage,dtdzlow,ustr,vstr,ustk,vstk,wavelen,waveamp,
     *   waveang,heaton,uston,vston,wavelon,wavehon,angleon,swheaton,
     *   rain,rainf,rainon,lheat,lheaton
     *   ,rad,dsdzlow,dtdztop,dsdztop,bvft
c
      common/bcmask/ hf_tmask,lf_tmask,str_tmask,str_tbc,
     *   hfton,lfton,strton
