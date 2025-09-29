      subroutine dag_write(n,nerr)
c
c   $Id: dag_write.F,v 1.13 1999/09/20 23:07:02 dwd Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
#include "cyclic.inc"
#include "finitv.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif
#include "communicate.inc"
c 
#include "moddef.inc"
c 
      integer time,strtt(4),countt(4),countu(4),countw(4)
      integer n,k,nerr,rcode,ix,kx,jx,i,j
#if defined(USEMPI)
      integer ierr
      real*8 pew
#endif
      logical first
      real zi,q,qs,beta,ws,ts,sigthta,s
      real ustar,tau,rainfall,buf(11)
      real la, S0,tave1(kdim),kmave(kdim1),w2ave(kdim1)
      real tave(kdim1),salave(kdim1),rhopave(kdim1),cave(kdim1)
      real uaved(0:kdim1),vaved(0:kdim1)
      real sig0,sig1,sig2,zz,urt,urb
      real*8 ep,ept,pe,pet
      real rr(kdim),z,pes,pets,z1
      real pews,rave
      real*8 uwave(kdim1),vwave(kdim1),udif,vdif
#if defined(BACKPE)
      real rhotot(0:idim1,0:jtotal1,0:kdim1)
      real rhoin(idim*jtotal*kdim)
      integer indr(idim*jtotal*kdim),nn
#endif
      real rtemp(0:idim1,0:jdim1,0:kdim1)
#if defined(cray)
      real darea(kdim1),area(kdim1),dave
      real ave(kdim1)
#if defined(USEMPI)
      real sd(kdim),sd1(kdim1)
#endif
#else
      real*8 darea(kdim1),area(kdim1),dave,ave1,areai,areaid
      real*8 lafact(kdim1)
      real*8 ave(kdim1), ave2(kdim1)
#if defined(USEMPI)
      real*8 sd(kdim),sd1(kdim1)
#endif
#endif
      real sarea,tarea
      real ttop,tbot
      integer indx(kdim)
c
      data countt/1,1,1,1/
      data countu/1,1,1,1/
      data countw/1,1,1,1/
      data strtt/1,1,1,1/
      data first/.true./
c
      countu(3)=kdim
      countw(3)=kdim1
c

#if defined(ICE)
      sarea = 0
      do i=1,idim
        do j=1,jdim
          if(hice(i,j).gt.0.0) then
           sarea = sarea+1
          endif
        enddo
      enddo
      do k=1,kdim
        lafact(k) = 0.0
        do j=1,jdim
          do i=1,idim
            if(k.le.kdist(i,j).or.hice(i,j).eq.0.0) then
               lafact(k) = lafact(k)+1.
            endif
          enddo
        enddo
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(sarea,tarea,1,
     &        MPI_REAL,MPI_SUM,mastertask,
     &        comm, ierr)
      areai = tarea
      call MPI_REDUCE(lafact,area,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
#else
      areai = sarea
      do k=1,kdim1
        area(k) = lafact(k)
      enddo
#endif
#else
      areai = 1.
      do k=1,kdim
        area(k) = float(idim*jtotal)
      enddo
#endif
      if(areai.eq.0.0) areai = 1.0
      areaid = areai*dble(ndiag)
c
      area(kdim1) = area(kdim)
      do k=1,kdim1
        if(area(k).eq.0.0) area(k) = 1.0
        darea(k) = area(k)*dble(ndiag)
      enddo

c      area=idim*jtotal
c      darea = area*dble(ndiag)
c
c write output fields
c 
      if(iope) then
        time=n*delt
        if(first)then
           first=.false.
           tindext = 1
           strtt(4) = 1
           call ncvpt(timid, ttid, 1, 1, time,rcode)
        else
           tindext = tindext + 1
           strtt(4)=tindext
           call ncvpt(timid, ttid, tindext, 1, time, rcode)
        endif
      endif
      ix=diloc
      jx=djloc
      kx=dkloc
c
c u (actually um)
c
#if defined(USEMPI)
      call collect_point(ix,jx,kx,buf)
#else
      buf(1) = um(ix,jx,kx)
      buf(2) = vm(ix,jx,kx)
      buf(3) = wm(ix,jx,kx)
      buf(4) = p(ix,jx,kx)
      buf(5) = th(ix,jx,kx)
      buf(6) = sal(ix,jx,kx)
      buf(7) = tke(ix,jx,kx)
      buf(8) = kheat(ix,jx,kx)
      buf(9) = kmom(ix,jx,kx)
      buf(10) = rhop(ix,jx,kx)
#if defined(TRACER)
      buf(11) = trc(ix,jx,kx,1)
#endif
#endif
      if(iope) then
        call ncvpt(timid,utid,strtt,countt,buf(1),rcode)
c
c v (actually vm)
c
        call ncvpt(timid,vtid,strtt,countt,buf(2),rcode)
c
c w (actually wm)
c
        call ncvpt(timid,wtid,strtt,countt,buf(3),rcode) 
c
c p_p
c
        call ncvpt(timid,ptid,strtt,countt,buf(4),rcode)
c
c t
c
        call ncvpt(timid,thtid,strtt,countt,buf(5),rcode)
c
c s
c
        call ncvpt(timid,saltid,strtt,countt,buf(6),rcode)
c
c tke*tke
c
        call ncvpt(timid,tketid,strtt,countt,buf(7),rcode)
c
c kh
c
        call ncvpt(timid,khtid,strtt,countt,buf(8),rcode)
c
c km
c
        call ncvpt(timid,kmtid,strtt,countt,buf(9),rcode)
c
c rho_p
c
        call ncvpt(timid,rho_ptid,strtt,countt,buf(10),rcode)
c
c trcr (only output 1 at this time, needs modification if more than
c one tracer is examined)
c
        do 20 k=1,ntrace
           call ncvpt(timid,trcrtid(k),strtt,countt,buf(11),rcode)
 20     continue
      endif
c
c zi, have to do this over all processers because of averaging
c
      call zmixed(zi)
      if(iope) then
        call ncvpt(timid,zitid,strtt,countt,zi,rcode)
c
c compute w_star, u_star, t_star, q, rainfall
c these are from j=1 on the mastertask

        q=0.0
        j=1
        do 30 i=1,idim
           q=q+heatf*randn(i,j)
 30     continue
        q=q/float(idim)
c
        qs=-q/(cp*rho0)
c
        s=sal(1,j,kdim)
        beta=-(sigthta(s,th(1,j,kdim)+0.001,0.0) - 
     *       sigthta(s,th(1,j,kdim),0.0))/(0.001*rho0)
c
        if(beta.le.0 .or. qs.le.0)then
           ws=1.
           ts=1.
        else
           ws=(beta*grav*zi*qs)**0.3333333
           ts=qs/ws
        endif
        tau = sqrt(ustr**2 + vstr**2)
      endif
c
c u_star
c
c
c ustar
c
#if defined(USEMPI)
      call MPI_REDUCE(ustara,ave1,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        ustar = ave1/darea(kdim)
      endif
#else
      ustar = ustara/darea(kdim)
#endif
      if(iope) then
        call ncvpt(timid,u_startid,strtt,countt,ustar,rcode)
      endif
      ustara = 0.0

c
c compute S0 and La (without km term)
c
      if(iope) then
#if defined(STOKES)
c
c        S0 = 0.5*((3.14159*waveamp/wavelen)**2)*sqrt(1.561*wavelen)
        S0 = sqrt(us(kdim1)**2 + vs(kdim1)**2+0.0001)
        if(tau .eq. 0.0) then
           La = 0.0
        else
           La = (6.283185/(wavelen*ustar))**1.5/sqrt(S0/ustar)
        endif
#endif
c
        if(rainon)rainfall = 3600.0*rainf*1000.0
c
c w_star
c
        call ncvpt(timid,w_startid,strtt,countt,ws,rcode)
c
c t_star
c
        call ncvpt(timid,t_startid,strtt,countt,ts,rcode)
c
c q
c
        call ncvpt(timid,qtid,strtt,countt,q,rcode)
c
c hf_top
c
        if(heaton) call ncvpt(timid,hf_toptid,strtt,countt,
     *     heatf,rcode)
c
c ustr_t
c
        if(uston) call ncvpt(timid,ustr_ttid,strtt,countt,
     *     ustr,rcode)
c
c vstr_t
c
        if(vston) call ncvpt(timid,vstr_ttid,strtt,countt,
     *     vstr,rcode)
c
c swf_top
c
        if(swheaton) call ncvpt(timid,swf_toptid,strtt,countt,
     *     swheatf,rcode)
c
c rain
c
        if(rainon) call ncvpt(timid,raintid,strtt,countt,
     *     rainfall,rcode)
c
c lhf_top
c
        if(lheaton) call ncvpt(timid,lhf_toptid,strtt,countt,
     *     latent,rcode)
c
#if defined(STOKES)
c
c wave_l
c
        if(wavelon) call ncvpt(timid,wave_ltid,strtt,countt,
     *     wavelen,rcode)
c
c wave_h
c
        if(wavehon) call ncvpt(timid,wave_htid,strtt,countt,
     *     waveamp,rcode)
c
c w_angle
c
        if(angleon) call ncvpt(timid,w_angletid,strtt,countt,
     *     waveang,rcode)
#endif
      endif
c
c t_bar
c
      call intg(th,thbar)
      if(iope) then
        call ncvpt(timid,t_bartid,strtt,countt,thbar,rcode)
      endif
c
c s_bar
c
#if defined(SNGLSALT)
      call intg(sal,salbar)
#else
      call intgdp(sal,salbar)
#endif
      if(iope) then
        call ncvpt(timid,s_bartid,strtt,countt,salbar,rcode)
      endif
c
c write the horizontal average fields
c
c
c t_ave
c
      call horzavg(th,tave)
c
      if(iope) then
        call ncvpt(timid,t_avetid,strtt,countu,tave,rcode)
      endif
#if defined(ICE)
#if !defined(WARMICE)
c
c cfr_ave
c
      call horzavg(cfr,cave)
c
      if(iope) then
        call ncvpt(timid,cfratid,strtt,countu,cave,rcode)
      endif
c
c cfr flux wc_ave
c
      do k=1,kdim
        dave=0.0
          do j=1,jdim
            do i=1,idim
               dave=dave+thp(i,j,k)*(cfr(i,j,k)-cave(k))
            enddo
          enddo
          ave(k)=dave
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
c
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/area(k)
        enddo
        call ncvpt(timid,wctid,strtt,countu,tave1,rcode)
      endif
#else
      do k=1,kdim
         tave1(k)= ave(k)/area(k)
      enddo
      call ncvpt(timid,wctid,strtt,countu,tave1,rcode)
#endif
#endif
#endif
c
c t2_ave
c
      do k=1,kdim
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.le.kdist(i,j)) then
                 dave=dave+(th(i,j,k)-tave(k))**2
               endif
            enddo
         enddo
         ave(k)=dave
      enddo
c
#if defined(USEMPI) 
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k)= ave(k)/area(k)
      enddo
#endif

      if(iope) then
        call ncvpt(timid,t2_avetid,strtt,countu,tave1,rcode)
      endif
c
c s_ave
c
#if defined(SNGLSALT)
      call horzavg(sal,salave)
#else
      call horzavgdp(sal,salave)
#endif
c
      if(iope) then
        call ncvpt(timid,s_avetid,strtt,countu,salave,rcode)
      endif
c
c s2_ave
c
      do k=1,kdim
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.le.kdist(i,j)) then
                 dave=dave+(sal(i,j,k)-salave(k))**2
               endif
            enddo
         enddo
         ave(k)=dave
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,s2_avetid,strtt,countu,tave1,rcode)
      endif
c
c tke_ave
c
      do k=1,kdim
         ave(k) = 0.5*(u2_ave(k)+v2_ave(k)+w2_ave(k))
      enddo
c
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*sd(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = ave(k)/darea(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,tke_avetid,strtt,countu,tave1,rcode)
      endif
c
c km_ave 
c
      call horzavg(kmom,kmave)
c
      if(iope) then
        call ncvpt(timid,km_avetid,strtt,countu,kmave,rcode)
c
#if defined(STOKES)
c
c S_0
c
        call ncvpt(timid,S_0tid,strtt,countt,S0,rcode)
c
c L_a, add kmom term
c
        La = La*ave(kdim)**1.5
        call ncvpt(timid,L_atid,strtt,countt,La,rcode)
#endif
      endif
c
c w2_ave
c
c      do k=1,kdim1
c         dave=0.0
c         do j=1,jdim
c            do i=1,idim
c               dave=dave+wm(i,j,k)*wm(i,j,k)
c            enddo
c         enddo
c         ave(k)=dave
c      enddo
      w2_ave(kdim1) = 0.0
#if defined(USEMPI)
      call MPI_REDUCE(w2_ave,sd1,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          w2ave(k) = sd1(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim1
         w2ave(k)= w2_ave(k)/darea(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,w2_avetid,strtt,countw,w2ave,rcode)
      endif
c
c KE_ave
c
c
c compute wm*,  place in thp as a temporary array
c
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               thp(i,j,k)=0.5*(wm(i,j,k)+wm(i,j,k+1))
            enddo
         enddo
      enddo
c
      do k=1,kdim
         dave=0.0
         do j=1,jdim
            do i=1,idim
c               if(k.le.kdist(i,j)) then
                 dave=dave+um(i,j,k)*um(i,j,k)+
     *              vm(i,j,k)*vm(i,j,k)+0.5*(wm(i,j,k)*wm(i,j,k)+
     *              wm(i,j,k+1)*wm(i,j,k+1))
c               endif
            enddo
         enddo
         ave(k)=0.5*dave
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,KE_avetid,strtt,countu,tave1,rcode)
      endif

c
c wt_ave
c
c top grid point is either the external heat flux or the heatflux at
c the ice/ocean interface
c
c      do k=2,kdim1
c         dave=0.0
c         do j=1,jdim
c            do i=1,idim
c              if(k.eq.kdim1) then
c#if defined(ICE)
c                dave = 0.0
c#else
c                dave = dave-heatf*randn(i,j)/(cp*rho0)
c#endif
c              else
cc               dave=dave+thp(i,j,k)*(th(i,j,k)-tave(k))
c#if defined(ICE)
c               if(k.le.kdist(i,j)) then
c#else
c               if(k.gt.kdist(i,j) then
c#endif
c                 dave=dave+wm(i,j,k)*((th(i,j,k)-tave(k))
c     &                   +(th(i,j,k-1)-tave(k-1)))*0.5
c               endif
c              endif
c            enddo
c         enddo
c#if defined(ICE)
cc         if(k.eq.kdim1) then
cc           dave = -hicebc/dble(ndiag)
cc         endif
c#endif
c         ave(k)=dave
c      enddo
c      ave(1) = ave(2)
c#if defined(USEMPI)
c      call MPI_REDUCE(ave,sd1,kdim1,
c     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
c     &        comm, ierr)
c      if(iope) then
c        do k=1,kdim1
c          tave(k) = sd1(k)/area(k)
c        enddo
c      endif
c#else
c      do k=1,kdim1
c         tave(k)=ave(k)/area(k)
c      enddo
c#endif
c
#if defined(USEMPI)
      call MPI_REDUCE(wt_ave,ave2,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
#else
      do k=1,kdim
        ave2(k) = wt_ave(k)
      enddo
#endif
      if(iope) then
        do k=1,kdim
          tave(k) = ave2(k)/darea(k)
        enddo
        call ncvpt(timid,wt_avetid,strtt,countu,tave,rcode)
      endif

c      if(iope) then
c        call ncvpt(timid,wt_avetid,strtt,countw,tave,rcode)
c      endif
c
c ws_ave
c
      do k=2,kdim1
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.eq.kdim1) then
                 dave = dave+latent*sal(i,j,k)/(1000.*lv)
     &                  + rainf*sal(i,j,k)
               else
c               dave=dave+thp(i,j,k)*(sal(i,j,k)-salave(k))
                 dave=dave+wm(i,j,k)*((sal(i,j,k)-salave(k))
     &                   +(sal(i,j,k-1)-salave(k-1)))*0.5
               endif
            enddo
         enddo
         ave(k)=dave
      enddo
      ave(1) = ave(2)
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd1,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          tave(k) = sd1(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim1
         tave(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,ws_avetid,strtt,countw,tave,rcode)
      endif
#if !defined(BACKPE)
c
c wp_ave at k=bottom
c
         dave=0.0
         do j=1,jdim
            do i=1,idim
               dave=dave+thp(i,j,1)*p(i,j,1)
            enddo
         enddo
#if defined(USEMPI)
      call MPI_REDUCE(dave,pew,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          pews = pew/area(k)
        enddo
      endif
#else
      do k=1,kdim1
        pews = dave/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,pewtid,strtt,countt,pews,rcode)
      endif
#endif 
c
c wrho_ave
c
      call horzavg(rhop,rhopave)
c
c      do k=1,kdim
c         dave=0.0
c         do j=1,jdim
c            do i=1,idim
c		dave=dave+((rhop(i,j,k+1)+rhop(i,j,k)
c     *          )*wm(i,j,k+1)-(rhop(i,j,k)+
c     *          rhop(i,j,k-1))*wm(i,j,k))*
c     *         (float(k-1)-0.5)*dz*0.5*odz
c            enddo
c         enddo
c         ave(k)=dave
c      enddo
c#if defined(USEMPI)
c      call MPI_REDUCE(ave,sd,kdim,
c     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
c     &        comm, ierr)
c#endif
c
#if defined(USEMPI)
      call MPI_REDUCE(dpea,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = dpea(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,wrho_avetid,strtt,countu,tave1,rcode)
      endif
c
c do dpesg
c
#if defined(USEMPI)
      call MPI_REDUCE(dpesg,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope)then
         do k=1,kdim
            tave1(k) = ave(k)/darea(k)
         enddo
      endif
#else
      do k=1,kdim
         tave1(k) = dpesg(k)/darea(k)
      enddo
#endif
      if(iope) then
         call ncvpt(timid,sig_avetid,strtt,countu,tave1,rcode)
      endif

c
c rhop_ave
c
      if(iope) then
c
        call ncvpt(timid,rhop_avetid,strtt,countu,rhopave,rcode)
      endif
c
c u_ave
c
      call horzavg(um,tave)
      do k=1,kdim
        uaved(k) = tave(k)
      enddo
      uaved(0) = 2.*uaved(1) - uaved(2)
      uaved(kdim1) = 2*uaved(kdim)-uaved(kdim-1)
c
      if(iope) then
        call ncvpt(timid,u_avetid,strtt,countu,tave,rcode)
      endif
c
c u2_ave
c
c      do k=1,kdim
c         dave=0.0
c         do j=1,jdim
c            do i=1,idim
c               dave=dave+(um(i,j,k)-tave(k))**2
c            enddo
c         enddo
c         ave(k)=dave
c      enddo
#if defined(USEMPI)
      call MPI_REDUCE(u2_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave(k) = u2_ave(k)/darea(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,u2_avetid,strtt,countu,tave,rcode)
      endif
c
c v_ave
c
      call horzavg(vm,tave)
      do k=1,kdim 
        vaved(k) = tave(k) 
      enddo
      vaved(0) = 2.*vaved(1) - vaved(2)
      vaved(kdim1) = 2*vaved(kdim)-vaved(kdim-1)
c
      if(iope) then
        call ncvpt(timid,v_avetid,strtt,countu,tave,rcode)
      endif
c
c v2_ave
c
c      do k=1,kdim
c         dave=0.0
c         do j=1,jdim
c            do i=1,idim
c               dave=dave+(vm(i,j,k)-tave(k))**2
c            enddo
c         enddo
c         ave(k)=dave
c      enddo
#if defined(USEMPI)
      call MPI_REDUCE(v2_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave(k) = v2_ave(k)/darea(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,v2_avetid,strtt,countu,tave,rcode)
      endif
c
c hf_ave
c
      do k=2,kdim1
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.eq.kdim1) then
                dave = 0.0
               else if(k.le.kdist(i,j)) then
                dave=dave-((kheat(i,j,k)+kheat(i,j,k-1)+2.*khb)
     *                 *(th(i,j,k)-th(i,j,k-1)))*0.5
               endif
            enddo
         enddo
         ave(k)=odz*dave
      enddo
      ave(1) = ave(2)
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd1,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          tave(k) = sd1(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim1
         tave(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,hf_avetid,strtt,countw,tave,rcode)
      endif
c
c rhof_ave
c
      do k=1,kdim
         dave=0.0
         do j=1,jdim
            do i=1,idim
               z1 = dz*(k-1)+0.5*dz
               zz = dz*(kdim-k)+0.5*dz
c	       sig1 = rhob(k)+rhop(i,j,k)
c	       sig2 = rhob(k+1)+rhop(i,j,k+1)
c	       sig0 = rhob(k-1)+rhop(i,j,k-1)
c	       if(k.eq.kdim) sig2 = sig1+(sig1-sig0)
c	       if(k.eq.1) sig0 = sig1
               s = sal(i,j,k)
	       sig1 = sigthta(s,th(i,j,k),zz)
               s = sal(i,j,k+1)
	       sig2 = sigthta(s,th(i,j,k+1),zz+dz)
               s = sal(i,j,k-1)
               sig0 = sigthta(s,th(i,j,k-1),zz-dz)
	       urt = (kheat(i,j,k+1)+kheat(i,j,k))*(sig2-sig1)*0.5
	       urb = (kheat(i,j,k)+kheat(i,j,k-1))*(sig1-sig0)*0.5
c	       if(k.eq.kdim) urt = 0.0
               dave=
     *           dave+(urt-urb)*odz*odz*z1
            enddo
         enddo
         ave(k)=dave
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave(k) = sd(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim
         tave(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,rhof_avetid,strtt,countu,tave,rcode)
      endif
c
c sf_ave
c
      do k=2,kdim1
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.eq.kdim1) then
                 dave = 0.0
               else
                if(k.le.kdist(i,j)) then
                  dave=dave-((kheat(i,j,k)+kheat(i,j,k-1)+2.*kms)
     *                 *(sal(i,j,k)-sal(i,j,k-1)))*0.5
                endif
               endif
            enddo
         enddo
         ave(k)=odz*dave
      enddo
      ave(1) = ave(2)
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd1,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          tave(k) = sd1(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim1
         tave(k) = ave(k)/area(k)
      enddo
#endif
c
c compute temperature distribution
c
      do k=1,kdim1
        sd1(k) =0.0
      enddo
      do j=1,jdim
        do i=1,idim
          do k=1,kdim
c            ttop = 16.75+float(k)*0.5/float(kdim)
            ttop = 16.85+float(k)*0.3/float(kdim)
c            tbot = 16.75+float(k-1)*0.5/float(kdim)
            tbot = 16.85+float(k-1)*0.3/float(kdim)
            if(th(i,j,k).gt.tbot.and.th(i,j,k).lt.ttop) then
              sd1(k) = sd1(k)+1
            endif
          enddo
        enddo
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd1,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim1
          tave(k) = sd1(k)
        enddo
      endif
#else
      do k=1,kdim1
         tave(k) = ave(k)
      enddo
#endif

            
      if(iope) then
        call ncvpt(timid,sf_avetid,strtt,countw,tave,rcode)
      endif
c
c p2_ave
c
      do k=1,kdim
         dave=0.0
         do j=1,jdim
            do i=1,idim
               if(k.le.kdist(i,j))then
                 dave=dave+p(i,j,k)**2
               endif
            enddo
         enddo
         ave(k)=dave
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/area(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = ave(k)/area(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,p2_avetid,strtt,countu,tave1,rcode)
      endif
c
#if defined(ICE) || defined(MPOND)
c
c hbc
c
#if defined(MPOND)
#if defined(USEMPI)
      call MPI_REDUCE(bottoma,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        MPI_COMM_WORLD, ierr)
      bottoma = dave
      if(bottoma.eq.0.0) then
        bottoma = 1.0
        if(iope) write(*,*)"No bottom Area defined"
      endif
      if(iope) write(*,*)"bottoma ",bottoma
#else
      if(bottoma.eq.0.0) then
        bottoma = 1.0
        write(*,*)"No bottom area defined"
      endif
#endif
      areaid = bottoma
#endif
#if defined(USEMPI)
      call MPI_REDUCE(hicebc,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        rave = dave/areaid
      endif
#else
      rave = hicebc/areaid
#endif
      hicebc = 0.0
      if(iope) then
        call ncvpt(timid,hbctid,strtt,countt,rave,rcode)
      endif
c
c sbc
c
#if defined(USEMPI)
      call MPI_REDUCE(sicebc,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        rave = dave/areaid
      endif
#else
      rave = sicebc/areaid
#endif
      if(iope) then
        call ncvpt(timid,sbctid,strtt,countt,rave,rcode)
      endif
      sicebc = 0.0
c
c
c bvst
c
#if defined(USEMPI)
      call MPI_REDUCE(bvst,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        rave = dave/areaid
      endif
#else
      rave = bvst/areaid
#endif
      if(iope) then
        call ncvpt(timid,bvsttid,strtt,countt,rave,rcode)
      endif
      bvst = 0.0
c
c wall average salt flux
c
#if defined(USEMPI)
      call MPI_REDUCE(wallarea,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        MPI_COMM_WORLD, ierr)
      wallarea = dave
      if(wallarea.eq.0.0) then
        wallarea = 1.0
        if(iope) write(*,*)"No wall Area defined"
      endif
      if(iope) write(*,*)"wallarea ",wallarea
#else
      if(wallarea.eq.0.0) then
        wallarea = 1.0
        write(*,*)"No Wall area defined"
      endif
#endif
#if defined(USEMPI)
      call MPI_REDUCE(wsf,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        MPI_COMM_WORLD, ierr)
      if(iope) then
        rave = dave/wallarea
      endif
#else
      rave = wsf/wallarea
#endif
      if(iope) then
        call ncvpt(timid,wsatid,strtt,countt,rave,rcode)
      endif
      wsf = 0.0
c
c wall v_star
c
#if defined(USEMPI)
      call MPI_REDUCE(wvst,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        MPI_COMM_WORLD, ierr)
      if(iope) then
        rave = dave/wallarea
      endif
#else
      rave = wvst/wallarea
#endif
      if(iope) then
        call ncvpt(timid,wvsttid,strtt,countt,rave,rcode)
      endif
      wvst = 0.0

c
c wall average heat flux
c
#if defined(USEMPI)
      call MPI_REDUCE(whf,dave,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        MPI_COMM_WORLD, ierr)
      if(iope) then
        rave = dave/wallarea
      endif
#else
      rave = whf/wallarea
#endif
      if(iope) then
        call ncvpt(timid,wtatid,strtt,countt,rave,rcode)
      endif
      whf = 0.0
      wallarea = 0.0
      bottoma = 0.0
#endif
#if !defined(WARMICE)
c
c cheat
c
#if defined(USEMPI)
      call MPI_REDUCE(cheat,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
        tave1(k) = cheat(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,chtid,strtt,countu,tave1,rcode)
      endif
c
c sflux
c
#if defined(USEMPI)
      call MPI_REDUCE(sflux,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
        tave1(k) = sflux(k)/darea(k)
      enddo
#endif
      if(iope)then
        call ncvpt(timid,sftid,strtt,countu,tave1,rcode)
      endif

#endif
#endif
#if defined(ANALY)
c
c uf_ave
c
#if defined(USEMPI)
      call MPI_REDUCE(uf_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = uf_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,ufilt_avetid,strtt,countu,tave1,rcode)
      endif
c
c vf_ave
c
#if defined(USEMPI)
      call MPI_REDUCE(vf_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim 
          tave1(k) = ave(k)/darea(k) 
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = vf_ave(k)/darea(k)
      enddo
#endif   
      if(iope) then 
        call ncvpt(timid,vfilt_avetid,strtt,countu,tave1,rcode)
      endif
c
c wf_ave
c
#if defined(USEMPI)
      call MPI_REDUCE(wf_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim 
          tave1(k) = ave(k)/darea(k) 
        enddo 
      endif
#else
      do k=1,kdim
         tave1(k)=wf_ave(k)/darea(k)
      enddo
#endif   
      if(iope) then 
        call ncvpt(timid,wfilt_avetid,strtt,countu,tave1,rcode)
      endif
c
c disipation rate
c
#if defined(USEMPI)
      call MPI_REDUCE(disp_ave,sd,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = sd(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = disp_ave(k)/darea(k)
      enddo
#endif
c
      if(iope) then
        call ncvpt(timid,disptid,strtt,countu,tave1,rcode)
      endif
c
c pressure term
c
#if defined(USEMPI)
      call MPI_REDUCE(p_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim 
          tave1(k) = ave(k)/darea(k)
        enddo 
      endif
#else
      do k=1,kdim
         tave1(k) = p_ave(k)/darea(k)
      enddo
#endif   
      if(iope) then 
        call ncvpt(timid,p_avetid,strtt,countu,tave1,rcode)
      endif
c
c total advection term (advection plus shear production)
c
#if defined(USEMPI)
      call MPI_REDUCE(a_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,a_avetid,strtt,countu,tave1,rcode)
      endif
c
c buoyancy term
c
#if defined(USEMPI)
      call MPI_REDUCE(b_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = b_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,b_avetid,strtt,countu,tave1,rcode)
      endif
c
c subgrid term (dissipation plus subgrid transport)
c
#if defined(USEMPI)
      call MPI_REDUCE(sg_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = sg_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,sg_avetid,strtt,countu,tave1,rcode)
      endif
c
c stokes term (if used)
c
#if defined(STOKES)
#if defined(USEMPI)
      call MPI_REDUCE(sd_ave,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = sd_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,sd_avetid,strtt,countu,tave1,rcode)
      endif
#endif
c
c uvdu for filtered fields
c
#if defined(USEMPI)
      call MPI_REDUCE(vpupum,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,uvdutid,strtt,countu,tave1,rcode)
      endif
c
c uvdv for filtered fields
c
#if defined(USEMPI)
      call MPI_REDUCE(vpupvm,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,uvdvtid,strtt,countu,tave1,rcode)
      endif
c
c uwdu for filtered fields
c
#if defined(USEMPI)
      call MPI_REDUCE(upwpum,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,uwdutid,strtt,countu,tave1,rcode)
      endif
c
c vwdv for filtered fields
c
#if defined(USEMPI)
      call MPI_REDUCE(vpwpvm,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,vwdvtid,strtt,countu,tave1,rcode)
      endif
c
c wrhop for filtered fields
c
#if defined(USEMPI)
      call MPI_REDUCE(wrhop,ave,kdim,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &        comm, ierr)
      if(iope) then
        do k=1,kdim
          tave1(k) = float(ndiag)*ave(k)/darea(k)
        enddo
      endif
#else
      do k=1,kdim
         tave1(k) = a_ave(k)/darea(k)
      enddo
#endif
      if(iope) then
        call ncvpt(timid,wrhotid,strtt,countu,tave1,rcode)
      endif





      do k=1,kdim
        uf_ave(k)=0.0
        vf_ave(k)=0.0
        wf_ave(k)=0.0
        a_ave(k)=0.0
        b_ave(k)=0.0
        p_ave(k)=0.0
        sg_ave(k)=0.0
        sd_ave(k)=0.0
        u2_ave(k) = 0.0
        v2_ave(k) = 0.0
        w2_ave(k) = 0.0
        disp_ave(k) = 0.0
        tke_ave(k) = 0.0
        dpea(k) = 0.0
        dpesg(k) = 0.0
#if defined(ICE)
        sflux(k) = 0.0
        cheat(k) = 0.0
#endif
        wt_ave(k) = 0.0
        vpupum(k) = 0.0
        vpupvm(k) = 0.0
        upwpum(k) = 0.0
        vpwpvm(k) = 0.0
        wrhop(k) = 0.0
      enddo
c
c
c next compute the shear production term separately
c
c     u'w' du/dz + v'w' dv/dz
c
c
      do k=1,kdim
         uwave(k) = 0.0
         vwave(k) = 0.0
         ave(k) = 0.0
         do j=1,jdim
            do i=1,idim
              if(k.le.kdist(i,j)) then
                udif = um(i,j,k)-uaved(k)+um(i,j,k-1)-uaved(k-1)
                vdif = vm(i,j,k)-vaved(k)+vm(i,j,k-1)-vaved(k-1)
                uwave(k) = uwave(k)
     *          + 0.25*udif*(wm(i,j,k)+wm(i-1,j,k))
                vwave(k) = vwave(k)
     *          + 0.25*vdif*(wm(i,j,k)+wm(i,j-1,k))
              endif
            enddo
         enddo
      enddo
#if defined(USEMPI)
      call MPI_REDUCE(uwave,ave,kdim,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &            comm, ierr)
#else
      do k=1,kdim
        ave(k) = uwave(k)
      enddo
#endif

      if(iope) then
         do k=1,kdim
           uwave(k) = ave(k)/area(k)
         enddo
      endif
#if defined(USEMPI)
      call MPI_REDUCE(vwave,ave,kdim, 
     &            MPI_DOUBLE_PRECISION,MPI_SUM,mastertask, 
     &            comm, ierr) 
#else
      do k=1,kdim
        ave(k) = vwave(k)
      enddo
#endif
      if(iope) then
        do k=1,kdim
          vwave(k) = ave(k)/area(k)
        enddo
      endif
c
c compute the full term
c
      if(iope) then
        uwave(kdim1) = 0.0
        vwave(kdim1) = 0.0
c
        do k=1,kdim
          ave(k) =   0.5*((uwave(k+1)+uwave(k))*
     *               (uaved(k+1)-uaved(k-1))*tdz
     *               + (vwave(k+1)+vwave(k))*
     *               (vaved(k+1)-vaved(k-1))*tdz)
        enddo
c
c output u'w' and v'w'
c
        do k=1,kdim
          tave(k) = uwave(k)
        enddo
        call ncvpt(timid,uw_avetid,strtt,countu,tave,rcode)
        do k=1,kdim 
          tave(k) = vwave(k) 
        enddo
        call ncvpt(timid,vw_avetid,strtt,countu,tave,rcode)
c
c output shear production term
c
        do k=1,kdim
          tave(k) = ave(k)
        enddo
        call ncvpt(timid,sp_avetid,strtt,countu,tave,rcode)
      endif
c
c next, do potential energy calculations
c
c first compute potential energy of flow with and without reordered 
c density profiles (column reordering ala Thorpe)
c
      ept = 0.0
      ep = 0.0
      do j=1,jdim
        do i=1,idim
          do k=1,kdim
            z = dz*(k-1)+0.5*dz
            zz = dz*(kdim-k)+0.5*dz
            s = sal(i,j,k)
            rr(k) = sigthta(s,th(i,j,k),zz)
            ep = ep+dble(rr(k)*z)
            rr(k) = -1.*rr(k)
          enddo
c
c next call indexing program to get reordered index
c
          call indexx(kdim,rr,indx)
c
c calculate p.e. using the reordered density
c
          do k=1,kdim
            z = dz*(k-1)+0.5*dz
            ept = ept-dble(rr(indx(k))*z)
          enddo
        enddo
      enddo
c
c next collect data from processors
c
#if defined(USEMPI)
      call MPI_REDUCE(ept,pet,1,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &            comm, ierr)

      call MPI_REDUCE(ep,pe,1,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
     &            comm, ierr)
#else
      pet = ept
      pe = ep
#endif
      if(iope) then
        pet = pet/dble(idim*jtotal*kdim)
        pe = pe/dble(idim*jtotal*kdim)
c        pets = pet/dble(sigthta(censal,centhta,0.0))
c        pes = pe/dble(sigthta(censal,centhta,0.0))
        pets = pet
        pes = pe
        call ncvpt(timid,pettid,strtt,countt,pets,rcode)
c        
        call ncvpt(timid,petid,strtt,countt,pes,rcode)
      endif
c
c next, compute base level p.e. by sorting 3-d density field
c
c first do node level sorting, then collect at mastertask and
c finish sort.
c
      do i=0,idim1
        do j=0,jdim1
          do k=1,kdim
c            rtemp(i,j,k) = -1.*sigthta(sal(i,j,k),th(i,j,k),0.0)
             rtemp(i,j,k) = rhop(i,j,k)+rhob(k)-1024.
          enddo
	  rtemp(i,j,0) = rtemp(i,j,1)
	  rtemp(i,j,kdim1) = rtemp(i,j,kdim)
        enddo
      enddo
c
c output the horizontal average of sigma theta
c
c      call horzavg(rtemp,tave)
c      do k=1,kdim
c        tave(k) = tave(k)*((k-1)*dz+0.5*dz)
c      enddo
c
c set top and bottom boundaries to zero, these points will get sorted 
c to the top of the index list since all defined densities are negative
c
#if defined(BACKPE)
      call collect_3d(rtemp,rhotot)
c
      nn = 1
      do k=1,kdim
        do j=1,jtotal
          do i=1,idim
            rhoin(nn) = rhotot(i,j,k)
            nn = nn+1
          enddo
        enddo
      enddo
      if(iope) then
        call indexx(kdim*idim*jtotal,rhoin,indr)
c
c start index count at 1 since all b.c. points are at the top of the
c index
c
        nn=1
        pew = 0.0
        do k=1,kdim
          z = dz*(k-1)+0.5*dz
          do j=1,jtotal
            do i=1,idim
              pew = pew-dble(z*rhoin(indr(nn)))
              nn = nn+1
            enddo
          enddo
        enddo
        pew = pew/dble(idim*jtotal*kdim)
        pews = pew/dble(sigthta(censal,centhta,0.0))
c
c output to dag file
        call ncvpt(timid,pewtid,strtt,countt,pews,rcode)
c
      endif
#endif
c
#endif

c sync the netCDF file
c
      if(iope) then
        call ncsnc(timid, rcode)
      endif
      return
      end
      subroutine horzavg(phi,avg)
c
c compute horizontal average of field phi and return the column
c 
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif
#include "finitv.inc"
#include "communicate.inc"

      real phi(0:idim1,0:jdim1,0:kdim1)
      real avg(kdim1)
      real*8 locavg(kdim1)
      integer i,j,k
      real ijsum(kdim1),lijsum(kdim1)
#if defined(USEMPI)
      real*8 lavg(kdim1)
      integer ierr
#endif
      do k=1,kdim1
         locavg(k) = 0.0
         lijsum(k) = 0
         do j=1,jdim
           do i=1,idim
#if defined(ICE)
             if(k.le.kdist(i,j).or.hice(i,j).eq.0.0) then
               locavg(k) = locavg(k) + phi(i,j,k)
               lijsum(k) = lijsum(k)+1
             endif
#endif
#if defined(HBOT)
             if(k.ge.kdist(i,j)) then
               locavg(k) = locavg(k) + phi(i,j,k)
               lijsum(k) = lijsum(k)+1
             endif
#endif
#if !defined(HBOT) && !defined(ICE)
               locavg(k) = locavg(k) + phi(i,j,k)
               lijsum(k) = lijsum(k)+1
#endif
           enddo
         enddo
      enddo
#if defined(USEMPI)
      call MPI_ALLREDUCE(locavg,lavg,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,
     &        comm, ierr)
      call MPI_ALLREDUCE(lijsum,ijsum,kdim1,
     &        MPI_REAL,MPI_SUM,
     &        comm, ierr)
      do k=1,kdim1
         if(ijsum(k).eq.0.0) ijsum(k) = 1.0
         avg(k) = lavg(k)/dble(ijsum(k))
      enddo
#else 
      do k=1,kdim1
         if(ijsum(k).eq.0.0) ijsum(k) = 1.0
#if defined(cray)
         avg(k) = locavg(k)/lijsum(k)
#else
         avg(k) = locavg(k)/dble(lijsum(k))
#endif
      enddo
#endif
      return
      end
      subroutine horzavgdp(phi,avg)
c
c compute horizontal average of field phi and return the column
c 
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "communicate.inc"

      real*8 phi(0:idim1,0:jdim1,0:kdim1)
      real avg(kdim1)
      real*8 locavg(kdim1)
      integer i,j,k
#if defined(USEMPI)
      real*8 lavg(kdim1)
      integer ierr
#endif
      do k=1,kdim1
         locavg(k) = 0.0
         do j=1,jdim
            do i=1,idim
               locavg(k) = locavg(k) + phi(i,j,k)
            enddo
         enddo
      enddo
#if defined(USEMPI)
      call MPI_ALLREDUCE(locavg,lavg,kdim1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,
     &        comm, ierr)
      do k=1,kdim1
         avg(k) = lavg(k)/dble(float(idim*jtotal))
      enddo
#else 
      do k=1,kdim1
#if defined(cray)
         avg(k) = locavg(k)/float(idim*jtotal)
#else
         avg(k) = locavg(k)/dble(float(idim*jtotal))
#endif
      enddo
#endif
      return
      end
