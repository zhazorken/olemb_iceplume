      subroutine sfcflx(tstep)
      use ocean_column
c
c   $Id: sfcflx.F,v 1.4 1998/10/21 21:56:42 dwd Exp $
c
c This subroutine computes the image points for scalars with a
c heat flux imposed at the model top
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c #include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "communicate.inc"
#include "finitv.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real rwtop,rwbot,ad,bb
      real znot,cd,ah,b,ch,rib,fh,vkrmn,rrn,znoth
      real psit,zet,psi0,fac,fact,cwan,am,ccmm,alzaz0,alzaz0h
      real psiq,fq,aq,heatl,qvsea,facq,temp
      real cn,cnh,fm,uistartemp,visk,cvisk,usqr,ypos
      real sfcheat,latheat,stress
      real uu,tair,qq,pp,z,zol,tskin,zomu,solar,qsfc,alpha
      real zomo,wu,wt,wq,uwind,vwind,zoq,zot,ztot
      real*8 uave8,vave8,tave8,qave8,sst8
      real*8 taveloc,uaveloc,vaveloc,qaveloc,sstloc
      real*8 shfloc,lhfloc,ustrloc
      real*8 tavebl(6),qavebl(6),tavebl8(6),qavebl8(6)
      real ubar,vbar,sstbar
      real ttbar,qbar,ctime
      real xportshf,xportlhf,rhow
      real spd,uflx,vflx,sflx,htflx,lhflx,swrad,orain,delto,hpl
      real uplw,deltsh
      real uc(1:kodim),vc(1:kodim),tc(1:kodim),sc(1:kodim)
      integer i,j,k,iswitch,tstep,jj,iopt_zo,iopt_sf,ierr,ns
c
c top and bottom boundary
c I have inserted placemark for future stress b.c. calculation
c
c define constants for MOS bottom fluxes
c
      aw = 17.27
      bw = 35.5
      b = 5.0
      ccmm = 7.5
      ch = 5.0
      visk = 1.45e-5
      cvisk = 0.11
      vkrmn = 0.4
      odz = 1./dzv(1)
c
      k=0
c
c compute the average surface winds
c
       shfloc = 0.0
       lhfloc = 0.0
       ustrloc = 0.0
       ztot = 0.0
       do k=1,10
         ztot = ztot+dzv(k)
       enddo
       do k=1,6
         tavebl(k) = 0.0
         qavebl(k) = 0.0
       enddo
#if defined(AVESST)
       sstloc =  0.0
       do j=1,jdim
         do i=1,idim
           sstloc = sstloc+tocn(i,j,1)+273.0
c           sstloc = sstloc+tsfc2(i,j)
         enddo
       enddo
#if defined(USEMPI)
       call MPI_ALLREDUCE(sstloc,sst8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)

c
        sstbar = sst8/dble(itotal*jtotal)
#endif
#endif
c       if(iope)write(*,*)"sstbar ",sstbar
#if defined(AVEFLUX)
       uaveloc = 0.0
       vaveloc = 0.0
       taveloc = 0.0
       qaveloc = 0.0
       ubar = 0.0
       vbar = 0.0
       ttbar = 0.0
       qbar = 0.0

       do j=1,jdim
          do i=1,idim
            uaveloc = uaveloc+um(i,j,1)
            vaveloc = vaveloc+vm(i,j,1)
            qaveloc = qaveloc+qv(i,j,1)
            taveloc = taveloc+th(i,j,1)
          enddo
        enddo
#if defined(USEMPI)
c
c call mpi to collect all sums
c
        call MPI_ALLREDUCE(uaveloc,uave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
        call MPI_ALLREDUCE(vaveloc,vave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
        call MPI_ALLREDUCE(taveloc,tave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
       call MPI_ALLREDUCE(qaveloc,qave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)

c
        ubar = uave8/dble(itotal*jtotal)
        vbar = vave8/dble(itotal*jtotal)
        qbar = qave8/dble(itotal*jtotal)
        ttbar = tave8/dble(itotal*jtotal)
#else
        ubar = uaveloc/dble(idim*jdim)
        vbar = vaveloc/dble(idim*jdim)
        ttbar = taveloc/dble(idim*jtotal)
        qbar = qaveloc/dble(idim*jtotal)
#endif
#endif
c
c move using Langrangian mean velocity
c
c       if(relstep.eq.1) then
c         xsst = 0.0
c         ysst = 0.0
c       endif
cc
c       xsst = xsst+uave*dt
c       ysst = ysst+vave*dt
c       if (abs(xsst).gt.dx) then
c          xsst = 0.0
c          do k=1,kodim
c           do j=1,jdim
c            do i=1,idim
c              tocn(i,j,k) = tocn(i-1,j,k)               
c              socn(i,j,k) = socn(i-1,j,k)               
c              uocn(i,j,k) = uocn(i-1,j,k)               
c              vocn(i,j,k) = socn(i-1,j,k)               
c            enddo
c           enddo
c          enddo
c       endif
c       if (abs(ysst).gt.dx) then
c          ysst = 0.0
c          do k=1,kodim
c           do j=1,jdim
c            do i=1,idim
c              tocn(i,j,k) = tocn(i,j-1,k)
c              socn(i,j,k) = socn(i,j-1,k)
c              uocn(i,j,k) = uocn(i,j-1,k)
c              vocn(i,j,k) = socn(i,j-1,k)
c            enddo
c           enddo
c          enddo
c       endif


          

c 
c bottom boundary heatflux condition
c
c hardwired for a warm land cold ocean
c
      do j=1,jdim
        do i=1,idim
c
c compute surface heat flux based on surface temperature
c
c              tsfc = 291.
c
c              znot = 0.25
c
c sea znot
c 
c              tsfc = 283.
c              tsfc = 296.
c
c coamps runs
c
#if defined(DYCOM)
               tsfc = 292.5
#else
c               tsfc = 293.
c               tsfc = 281.
#endif
c
c              tsfc = 290.
c              tsfc = 292.
c              tsfc = 300.
c               tsfc = 296.
#if defined(OPENBC)

              jj = jblkstart(my_pe+1)+j
              znot = 0.00001
c
c set up linear dt/dx of 10 degrees over 120 km after recirc .
c first case is cool to warm
c
#endif
c              tsfc = 305.
c              tsfc = 286.0
              ctime=tstep*delt
c              if(ctime.gt.300.and.ctime.lt.2100) then
c                ypos = ugeo(1)*(ctime-300.)
c                tsfc = 286.0+ypos*1.0/3000.
c                if(tsfc.gt.287) tsfc = 287.
c              elseif(ctime.gt.2100) then
c                 tsfc = 287.
c              endif
c              tsfc = 305.
c
c compute ocean column every 10 time steps
c
#if defined(TSFCDEF)
              if(ctime.gt.7200) tsfc = tsfc2(i,j)
#endif
              if(tstep.eq.1) then
                  tsfc2(i,j) = tocn(i,j,1)+273.
              endif
              if(mod(ctime,delt).eq.0.or.relstep.eq.1) then
              uplw = 0.0
#if defined(MLOCEAN)
                uplw = 5.67e-8*tsfc2(i,j)**4
c                   tsfc2(i,j) = tsfc2(i,j)+(glw(i,j)+0.5*swdn(i,j)
c     *              -botheat(i,j)
c     *              -uplw-botlhf(i,j))*0.33*delt/(1025.*4000.) 
c              else
c                tsfc2(i,j) = 305.5
c              endif
              
              do k=1,kodim
               uc(k) = uocn(i,j,k)
               vc(k) = vocn(i,j,k)
               tc(k) = tocn(i,j,k)
               sc(k) = socn(i,j,k)
              enddo
              hpl = hocn(i,j)
              htflx = botheat(i,j)+uplw-glw(i,j)
              lhflx = botlhf(i,j)
              swrad = swdn(i,j)
              orain = raint(i,j)/1000.
              raint(i,j) = 0.0
              spd = sqrt((um(i,j,1)+uave)**2+(vm(i,j,1)+vave)**2)
c              spd = sqrt(um(i,j,1)**2+vm(i,j,1)**2)
              uflx = 0.0012*uistar(i,j)*uistar(i,j)*um(i,j,1)/spd
              vflx = 0.0012*uistar(i,j)*uistar(i,j)*vm(i,j,1)/spd
              uflx = uflx+0.00001
              vflx = vflx+0.00001
              if(relstep.eq.1) then
                htflx = 0.0
              endif
c
              if(htflx.gt.1000) then
                write(*,*)"big ht flux ",htflx,uplw,glw(i,j),i,j,my_pe
              endif
              if(my_pe.eq.129.and.j.eq.53.and.i.eq.30) then
                write(*,*)"htflx,lhflx,swrad ",htflx,lhflx,swrad
                write(*,*)"hpl,uflx,vflx ",hpl,uflx,vflx,orain,delt
                write(*,*)"tc(1),uc(1-2)",tc(1),uc(1),uc(2)
                write(*,*)"sc(1),sc(2) ",sc(1),sc(2)
              endif
              if(tc(1).gt.38.or.tc(1).lt.10) then
                write(*,*)"bogus temp,delt,hpl,orain",tc(1),hpl,orain
                write(*,*)"htflx,lhflx,swrad ",htflx,lhflx,swrad
                write(*,*)"uflx,vflx ",uflx,vflx
                write(*,*)"tc(1)",tc(1),i,j,my_pe
                call MPI_ABORT(comm,1,ierr)
c                  write(*,*)"kt ",(kt(k),k=1,kodim)
c                  write(*,*)"uc ",(uc(k),k=1,kodim)
c                  write(*,*)"ghatt ",(ghatt(k),k=1,kodim)
              endif
c
              call ocn_col(uc,vc,tc,sc,uflx,vflx,htflx,lhflx,swrad,
     *                     orain,delt,hpl,1,i,j,tomegaz)
              
c             if(iope.and.j.eq.10.and.i.eq.10) then
c                write(*,*)"tc(1),uc(1-2)",tc(1),uc(1),uc(2)
c                write(*,*)"km,kt,ghatu ",(km(k),kt(k),ghatu(k),k=1,4)
c              endif

c              if(sc(1).gt.36.or.sc(1).lt.2.0) then
c                write(*,*)"bogus sal,delt,hpl,orain",delt,hpl,orain
c                write(*,*)"htflx,lhflx,swrad ",htflx,lhflx,swrad
c                write(*,*)"uflx,vflx ",uflx,vflx
c                write(*,*)"sc(1)",sc(1),sc(2),i,j,my_pe
cc                  write(*,*)"kt ",(kt(k),k=1,kodim)
cc                  write(*,*)"uc ",(uc(k),k=1,kodim)
cc                  write(*,*)"ghatt ",(ghatt(k),k=1,kodim)
c              endif
c              if(my_pe.eq.129) then
c                if(i.eq.30.and.j.eq.53) then
c                  do k=1,kodim
c                   write(*,*)"t,s,ks,kt ",tc(k),sc(k),ks(k),kt(k)
c                  enddo
c                endif
c              endif

              do k=1,kodim
               uocn(i,j,k) = uc(k)
               vocn(i,j,k) = vc(k)
               tocn(i,j,k) = tc(k)
               socn(i,j,k) = sc(k)
              enddo
              hocn(i,j) = hpl
              
              lwupf = lwupf+uplw
#endif
              endif
              tsfc2(i,j) = tc(1)+273.15
#if defined(BIGDOMAIN)
c              if(n*delt.gt.43200) then
c               tsfc = tsfc2(i,j)
c              else
c               tsfc = 303.0
c              endif
#endif
              tsfcave = tsfcave+tsfc2(i,j)
#if defined(AVESST)
              tsfc2(i,j) = sstbar
              tsfc = sstbar
#else

              tsfc = tsfc2(i,j)
#endif
c
c              tsfc = 287.
c              jj = blockstart(my_pe+1)+j
c              if(tstep*delt.lt.3600) then
c                if(jj.gt.450.and.jj.lt.1050)then
c                  if(i.gt.225.and.i.lt.425) then
c                    tsfc = 303.5
c                  endif
c                endif
c              endif
c              if(jj.ge.(jdim*8)+1) then
c                ypos = (jj-jdim*8-1)*dy
c                if(ypos.lt.120000) then
c                   tsfc = 281.0+ 10.*ypos/120000.
c                   tsfc = 281.0+ 8.*ypos/80000.
c                else
c                   tsfc = 289.
c                   tsfc = 291.
c                endif
c              else
c                tsfc = 281.0
c              endif
c
c this case is warm to cold, land to water
c land = 303, water = 290
c
c              if(jj.ge.(jdim*6)+1) then
cc                ypos = (jj-jdim*6-1)*dy
cc                if(ypos.lt.10000) then
cc                   tsfc = 293.0+ 5.*ypos/10000.
cc                else
c                   tsfc = 290.
cc                endif
c              else
c                tsfc = 301.0
c              endif
c
c               tsfc = 293.
c
c
         uwind=um(i,j,1)+uave
         vwind=vm(i,j,1)+vave
         tair=th(i,j,1)/((100000./(pbar(1)+p(i,j,1)))**0.286)
         if(iope.and.tair.gt.320.) then
            write(*,*)"tair,tsfc ",tair,tsfc2(i,j)
         endif
#if defined(AVEFLUX)
c         uwind = ubar+uave
c         vwind = vbar+vave
c         tair = ttbar/((100000./pbar(1))**0.286)
#endif
         uu=sqrt(uwind*uwind+vwind*vwind)
c
         qq = 0.0
#if defined(MOIST)
         qq=qv(i,j,1)
#if defined(AVEFLUX)
c         qq = qbar
#endif
#endif
         pp = (p(i,j,k)+pbar(1))/100.
         solar=0.
         tskin=tsfc
         if(tskin.lt.1.0) then
          write(*,*)"tskin too small ",tskin,i,j,iblkstart(my_pe+1)
          call MPI_ABORT(comm,1,ierr)
         endif
         zomu=0.
#if defined(MOIST)
         qvsea = (380.*exp(aw*(tsfc-273.16)/
     *          (tsfc-bw)))/(pbar(1)+p(i,j,1))
#else
         qvsea = 0.0
#endif
         qsfc=qvsea
         alpha=0.011
c         iopt_zo=2
c         iopt_zo=7
         iopt_zo=9
c         iopt_sf=1
         iopt_sf=2
         z = 0.5*dzv(1)
#if defined(LAND)
         qsfc = (380.*exp(aw*(280.-273.16)/
     *          (tsfc-bw)))/(pbar(1)+p(i,j,1))
         tskin = 315.
         iopt_zo=6
         iopt_sf=1
#endif
         call mosmodel (z,uu,tair,qq,pp, solar,tskin,zomu,
     &      qsfc,alpha,iopt_zo,iopt_sf, zol,zomo, wu,wt,wq,
     &      stress,sfcheat,latheat,zot,zoq)
         latheata = latheata+latheat
         sfcheata = sfcheata+sfcheat
         lhfloc = lhfloc+latheat
         shfloc = shfloc+sfcheat
         ustrloc = stress+ustrloc
#if defined(AVEFLUX)
         do k=1,6
           tavebl(k) = tavebl(k)+th(i,j,k)     
           qavebl(k) = qavebl(k)+qv(i,j,k)     
         enddo
#endif
c
         if(iope.and.i.eq.2.and.j.eq.2) then
           write(14,*)"t,rho,pp,q,qsfc,tau,u,zt,zq",tair,rhob(1),pp,qq,
     *                 qsfc,stress,uu,zot,zoq,ztot
         write(*,*)"Heat Fluxes,tsfc,htf,lh...randn ",
     *              sfcheat,latheat,-glw(i,j),
     *              uplw,swdn(i,j),tocn(i,j,1),htflx,lhflx,randn(i,j)
         endif
#if !defined(AVESTRESS)
         uistar(i,j) = sqrt(stress/rhob(1))
         ustrsum(i,j) = uistar(i,j)
#endif
c
         xportlhf = 0.
c         if(iblkstart(my_pe).lt.itotal/2) then
c           if(jblkstart(my_pe).lt.jtotal/2) then
c              xportlhf = 450.
c           endif
c         endif
         xportshf = 0.
c         xportshf = 50.
#if defined(LAND)
         xportlhf = 0.
#endif
 610      format(8(1x,e8.2))
c          if(i.eq.1.and.j.eq.1 ) then
c            write(12,610)um(i,j,1),vm(i,j,1),uistar(1,1),latheat,
c     *                 sfcheat,p(i,j,1),th(i,j,1),qv(i,j,1)
c          endif
#if defined(MOIST)
c
c turn off heat flux
c
c          sfcheat = 0.0
          rhow = 0.622*rhob(1)/(0.622+qv(i,j,1))
          if(tstep*delt.gt.120.) then
c            if(latheat.lt.0.0) then
c              qv(i,j,1)=qv(i,j,1)
c     *        +0.5*delt*randn(i,j)*latheat*odz/(lv*rhob(1))
c            else
c
          botlhf(i,j) = latheat*randn(i,j)
#if !defined(AVEFLUXLH)
              qv(i,j,1)=qv(i,j,1)
     *           +delt*(latheat)/(dzv(1)*lv*rhow)
c              do k=1,10
c                rhow = 0.622*rhob(k)/(0.622+qv(i,j,k))
c                qv(i,j,k) = qv(i,j,k)
c     *           +delt*(xportlhf)/(ztot*lv*rhow)
c              enddo
c              if(iope.and.i.eq.2.and.j.eq.2) then
c                write(*,*)"lath,qv",latheat,qv(i,j,1), delt,odz,lv,
c     *          rhob(1),xportlhf
c              endif
              botlhf(i,j) = latheat*randn(i,j)
#endif
c            endif
            wqvs_ave(1) = wqvs_ave(1)+
     *               (latheat+xportlhf)/(lv*rhob(1))
          endif
#endif
#if !defined(AVEFLUXSH)
          th(i,j,1) = th(i,j,1)+delt*randn(i,j)*
     *           (xportshf+sfcheat)/(cp*rhow*dzv(1))
#if defined(BOTTOM)
          botheat(i,j) = sfcheat*randn(i,j)
#endif
#endif
          wts_ave(1) = wts_ave(1)+
     *            (xportshf+sfcheat)*randn(i,j)/(cp*rhow)
c
        enddo
      enddo
#if defined(AVESTRESS)
      call MPI_ALLREDUCE(ustrloc,vave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
      stress = vave8/dble(itotal*jtotal)
      do i=1,idim
        do j=1,jdim
           uistar(i,j) = sqrt(stress/rhob(1))
           ustrsum(i,j) = uistar(i,j)
        enddo
      enddo
          
#endif
c
#if defined(AVEFLUXLH)
#if defined(USEMPI)
c
c call mpi to collect all sums
c
      call MPI_ALLREDUCE(lhfloc,uave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
c      call MPI_ALLREDUCE(shfloc,vave8,1,
c     &        MPI_REAL8,MPI_SUM,comm, ierr)

c
c      sfcheat = vave8/dble(itotal*jtotal)
      latheat = uave8/dble(itotal*jtotal)
#else
c      sfcheat = shfloc/dble(idim*jdim)
      latheat = lhfloc/dble(idim*jdim)
#endif

      do i=1,idim
        do j=1,jdim
           rhow = 0.622*rhob(1)/(0.622+qv(i,j,1))
c          do k=1,10
            qv(i,j,1)=qv(i,j,1)
     *        +delt*latheat/(dzv(1)*lv*rhow)
c            th(i,j,1) = th(i,j,1)+delt*randn(i,j)*
c     *           (xportshf+sfcheat)/(cp*rhow*dzv(1))
c          enddo

c          botlhf(i,j) = latheat*randn(i,j)
c           
c          botheat(i,j) = wm(i,j,2)*(th(i,j,1)-ttbar)*cp*rhob(1)
        enddo
      enddo
#endif
c
#if defined(AVEFLUXSH)
#if defined(USEMPI)
c
c call mpi to collect all sums
c
      call MPI_ALLREDUCE(shfloc,vave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)

c
      sfcheat = vave8/dble(itotal*jtotal)
#else
      sfcheat = shfloc/dble(idim*jdim)
#endif
      do i=1,idim
        do j=1,jdim
           rhow = 0.622*rhob(1)/(0.622+qv(i,j,1))
c          do k=1,10
            th(i,j,1) = th(i,j,1)+delt*randn(i,j)*
     *           (xportshf+sfcheat)/(cp*rhow*dzv(1))
c          enddo

#if defined(BOTTOM)
          botheat(i,j) = sfcheat*randn(i,j)
#endif
c           
c          botheat(i,j) = wm(i,j,2)*(th(i,j,1)-ttbar)*cp*rhob(1)
        enddo
      enddo

#endif
#if defined(AVEQV)
       qbar = 0.0

       do j=1,jdim
          do i=1,idim
            qaveloc = qaveloc+qv(i,j,1)
          enddo
        enddo
c
c call mpi to collect all sums
c
       call MPI_ALLREDUCE(qaveloc,qave8,1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)

c
       qbar = qave8/dble(itotal*jtotal)
       do j=1,jdim
          do i=1,idim
            qv(i,j,1) = qbar
          enddo
        enddo
#endif

      call imag(th)
      call imag(qv)
      return
      end
