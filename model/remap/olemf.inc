c
c   $Id: olemf.inc.in,v 1.3 1998/10/21 22:02:13 dwd Exp $
c
c     u     = horizontal velocity field 1
c     v     = horizontal velocity field 2
c     w     = vertical velocity field
c     p     = pert. pressure field
c     th    = temperature field
c     sal   = salinity field
c
      real u(0:idim1,0:jdim1,0:kdim1),um(0:idim1,0:jdim1,0:kdim1)
      real umm(0:idim1,0:jdim1,0:kdim1)
      real v(0:idim1,0:jdim1,0:kdim1),vm(0:idim1,0:jdim1,0:kdim1)
      real vmm(0:idim1,0:jdim1,0:kdim1)
      real w(0:idim1,0:jdim1,0:kdim1),wm(0:idim1,0:jdim1,0:kdim1)
      real wmm(0:idim1,0:jdim1,0:kdim1)
c
      real vf(0:idim1,0:jdim1,0:kdim1),wf(0:idim1,0:jdim1,0:kdim1)
      real uf(0:idim1,0:jdim1,0:kdim1)
      real pf(0:idim1,0:jdim1,0:kdim1)
      real rhopf(0:idim1,0:jdim1,0:kdim1)
      real dvadvv(0:idim1,0:jdim1,0:kdim1)
      real dcor(0:idim1,0:jdim1,0:kdim1)
      real dpres(0:idim1,0:jdim1,0:kdim1)
      real sgf(0:idim1,0:jdim1,0:kdim1)
c

      real umf(0:idimd1,0:jdimd1,0:kdimd1)
      real vmf(0:idimd1,0:jdimd1,0:kdimd1)
      real wmf(0:idimd1,0:jdimd1,0:kdimd1)
      real thf(0:idimd1,0:jdimd1,0:kdimd1)
      real ppf(0:idimd1,0:jdimd1,0:kdimd1)
      real salf(0:idimd1,0:jdimd1,0:kdimd1)
c
      real p(0:idim1,0:jdim1,0:kdim1)
      real pm(0:idim1,0:jdim1,0:kdim1)
      real pmm(0:idim1,0:jdim1,0:kdim1)
      real th(0:idim1,0:jdim1,0:kdim1),thp(0:idim1,0:jdim1,0:kdim1)
#if defined(SNGLSALT)
      real sal(0:idim1,0:jdim1,0:kdim1),salp(0:idim1,0:jdim1,0:kdim1)
#else
      real*8 sal(0:idim1,0:jdim1,0:kdim1),salp(0:idim1,0:jdim1,0:kdim1)
#endif
      real tke(0:idim1,0:jdim1,0:kdim1),tkep(0:idim1,0:jdim1,0:kdim1)
#if defined(TRACER)
      real trc(0:idim1,0:jdim1,0:kdim1,trcdim)
      real trcp(0:idim1,0:jdim1,0:kdim1)
#endif
      real rhop(0:idim1,0:jdim1,0:kdim1),rhob(0:kdim1),rhobm1(0:kdim1)
      real rho0
      real fluxvd(0:idim1,0:kdim1),fluxv0(0:idim1,0:kdim1)
      real fluxv1(0:idim1,0:kdim1)
      real hbot(0:idim1,0:jdim1)
      real uistar(0:idim1,0:jdim1)
      real uistart(0:idim1,0:jdim1)
#if defined(ANALY)
      real*8 p_ave(kdim1),a_ave(kdim1),b_ave(kdim1),sg_ave(kdim1)
      real*8 sd_ave(kdim1),uf_ave(kdim1)
      real u_bar(0:kdim1),v_bar(0:kdim1)
      real*8 vf_ave(kdim1),wf_ave(kdim1),u2_ave(kdim1),v2_ave(kdim1)
      real*8 w2_ave(kdim1),tke_ave(kdim1),disp_ave(kdim1)
      real*8 wt_ave(kdim1)
      real*8 dpea(kdim1),dpesg(kdim1),pen(kdim1),peo(kdim1)
      real*8 vpupum(kdim1), vpupvm(kdim1),upwpum(kdim1)
      real*8 vpwpvm(kdim1), wrhop(kdim1)
#endif
      real*8 ustara
      real salbar(0:kdim1),thbar(0:kdim1)
c
      common/comuvw/ u,um,umm,w,wm,wmm,v,vm,vmm
      common/comdub/umf,vmf,wmf,thf,ppf,salf
      common/comfilter/uf,vf,wf,pf,rhopf
      common/balance/dvadvv,dcor,dpres,sgf
      common/comt/  p,pm,pmm,th,thp,sal,salp,salbar,thbar 
#if defined(TRACER)
      common/comtrc/ trc,trcp
#endif
      common/comq/ tke,tkep
      common/comrho/ rhop,rhob,rhobm1,rho0
      common/comflx/fluxvd,fluxv0,fluxv1
      common/combot/ustara, hbot,uistar,uistart
#if defined(ANALY)
      common/analy/p_ave,a_ave,b_ave,sg_ave,sd_ave,uf_ave,vf_ave,
     $             wf_ave,u2_ave,v2_ave,w2_ave,
     $             dpesg,dpea,pen,peo,tke_ave,disp_ave,wt_ave,
     $             u_bar,v_bar,vpupum,vpupvm,upwpum,vpwpvm,wrhop
#endif
