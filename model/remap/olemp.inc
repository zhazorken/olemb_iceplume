c
c   $Id: olemp.inc,v 1.4 1998/11/04 00:25:33 dwd Exp $
c
      real tdz,tdy,tdx,dz,dy,dx,edz,edy,edx,odz,ody,odx,fdz,fdy,fdx
      real tdxs,tdys,tdzs
      real ds
      real delt,dtdx,dtdy,dtdz,start,stop
      real xlat,tomegaz,tomegay
      real grav,omega,cp,lv
      real mdvals,bcvals,bcfrce,sgsvals,ttime,icvals
      real filtfac
      real maxu,maxv,maxw,maxt
      real thtamin,thtamax,salmin,salmax,centhta,censal
      real iczi,icramp
      real uave,vave
      real zicriteria
      real raytau(0:kdim1),taumax
      real vgeo,ugeo
      logical isrestart,advsal,filtvert,doanal,rotate,useomegy
      logical appden,linden,potden,uselatent,havefloats
      logical useswheat,usestress,userain
      logical bcuseran,icuseran,usetopref,useray
      logical enddumpon
      integer timstp,filtord,filtskp
c      integer imax1,jmax1,kmax1,imax,jmax,kmax
      integer punit,iunit,ounit,dunit,lunit,runit,funit
      integer mdflgs,bcflgs,sgsflgs,icflgs
      integer notimes,nofld,ofield,nostart,nostop,nodt
      integer ndiag,diloc,djloc,dkloc
      integer nslab,siloc,sjloc,skloc,noslb,oslab
      real sxloc,syloc,szloc
      integer slbtimes,slbtype
      integer ntitle,stime,ntime,nrstrt,ntrace,nstart,relstep
      integer icrseed
      integer zimethod
      character bcfstr*80,bcvstr*80,bcfrcstr*80,bctime*30,bcver*80
      character bccom1*120,bccom2*120,bccom3*120
      character sgsfstr*80,sgsvstr*80
      character sgscom1*120,sgscom2*120,sgscom3*120
      character mver*80,title*80,trfile*80,trccom*120
      character icfstr*80,icvstr*80
c
      common/comdx/tdz,tdy,tdx,dz,dy,dx,edz,edy,edx,odz,ody,
     *     odx,fdz,fdy,fdx,tdzs,tdys,tdxs
      common/comds/ds
c      common/commx/imax,imax1,jmax,jmax1,kmax,kmax1
      common/comtx/delt,dtdx,dtdy,dtdz,start,stop
      common/com3d/xlat,tomegaz,tomegay
      common/timcom/stime,ntime,ttime,nrstrt,timstp,nstart,
     *     enddumpon,relstep
      common/const/grav,omega,cp,lv
      common/modflg/isrestart,advsal,filtord,filtskp,filtvert,
     *     doanal,rotate,useomegy,uselatent,havefloats,
     *     bcuseran,icuseran,usetopref,useswheat,usestress,
     *     userain,useray
      common/modval/filtfac
      common/maxadv/maxu,maxv,maxw,maxt
      common/mdunit/punit,iunit,ounit,dunit,lunit,runit,funit
      common/mdcom/mdflgs(20),mdvals(20)
      common/bcinp/bcfstr,bcvstr,bcfrcstr,bctime,bcver,
     *     bccom1,bccom2,bccom3
      common/bccom/bcflgs(20),bcvals(20),bcfrce(20)
      common/icinp/icfstr,icvstr
      common/iccom/icflgs(20),icvals(20)
      common/icval/iczi,icramp,icrseed
      common/sgsinp/sgsfstr,sgsvstr,sgscom1,sgscom2,sgscom3
      common/sgscom/sgsflgs(20),sgsvals(20)
      common/thtasal/thtamin,thtamax,salmin,salmax,
     *     appden,linden,potden,centhta,censal
      common/modvr/mver
      common/rwparmc/title(5)
      common/outcom/notimes,nofld(5),ofield(5,varoutdim),
     *     nostart(5),nostop(5),nodt(5)
      common/digcom/ndiag,diloc,djloc,dkloc
      common/slbcom/slbtimes,nslab(5),siloc(5),sjloc(5),skloc(5),
     *     noslb(5),slbtype(5),oslab(5,varoutdim),
     *     sxloc(5),syloc(5),szloc(5)
      common/ttlcom/ntitle
      common/trcinp/trfile(trcdim),trccom(trcdim)
      common/trccom/ntrace
      common/reference/uave,vave
      common/zitest/zimethod,zicriteria
      common/topdamp/raytau,taumax
      common/geost/ugeo,vgeo
