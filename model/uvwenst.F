      Subroutine uvwflux (dofilt,filtu,filtv,filtw,ru,rv,rw,n,nerr)
c
c  $Id: uvwflux.F,v 1.8 1999/08/16 20:22:32 eds Exp $
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "olemtim.inc"
#include "cyclic.inc"
#include "moddef.inc"
#include "communicate.inc"
#include "finitv.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif
c 
c  process all the data slab by slab
c
c
      logical dofilt
      integer j,n,nerr,i,k,kk
#if defined(USEMPI)
      integer mpierr
      real*8 u_bar8(kdim),v_bar8(kdim)
#endif
      real uadvu,wadvu,uadvw,wadvw,uadvv,wadvv,vadvu
      real ua,va,wa
      real filtu(0:idim1,0:jdim1,0:kdim1),filtw(0:idim1,0:jdim1,0:kdim1)
      real filtv(0:idim1,0:jdim1,0:kdim1)
      real dbx,dby,dbz,spd,junk
c
      real ru(0:idim1,0:jdim1,0:kdim1),rv(0:idim1,0:jdim1,0:kdim1)
      real rw(0:idim1,0:jdim1,0:kdim1)
      real dragub(idim,jdim)
      real draguh(idim,jdim)
      real dragu(idim,jdim)
      real dragvb(idim,jdim)
      real dragvh(idim,jdim)
      real dragv(idim,jdim)
c
c
c temporary arrays for fluxes
c
      real et1(0:idim1,0:jdim1,0:kdim1),et2(0:idim1,0:jdim1,0:kdim1)
      real et3(0:idim1,0:jdim1,0:kdim1),kin(0:idim1,0:jdim1,0:kdim1)
c
      real sgsru,sgsrw,sgsrv
      real b1,b2,b3,kmo
      real vadvv,vadvw
      real sd,cd,znot
      real uws(0:idim1,0:jdim1,0:kdim1)
      real vws(0:idim1,0:jdim1,0:kdim1)
#if defined(ANALY)
      real up(0:idim1,0:jdim1,0:kdim1)
      real vp(0:idim1,0:jdim1,0:kdim1)
      real uvs(0:idim1,0:jdim1,0:kdim1)
      real tkew(0:idim1,0:jdim1,0:kdim1)
      real disp(0:idim1,0:jdim1,0:kdim1)
      real w_ave(kdim1)
      real bw_ave(kdim1),sgw_ave(kdim1),sdw_ave(kdim1)
      real*8 ubr(kdim),vb(kdim)
      real*8 denom
#endif
#if defined(ICE) || defined(HBOT)
      real deptho
      integer olev
#endif
c
c
c setup prime u and v for analysis of tke budget
c
#if defined(ANALY)
      denom = dble(itotal*jtotal)
c      if(mod(n,ndiag).eq.0) then
        do k=1,kdim
          ubr(k) = 0.0
          do j=1,jdim
            do i=1,idim
              ubr(k)=um(i,j,k) + ubr(k)
            enddo
          enddo
#if defined(USEMPI)
          call MPI_ALLREDUCE(ubr(k),u_bar8(k),1,
     &         MPI_REAL8,MPI_SUM,
     &         comm, mpierr)
          u_bar(k) = u_bar8(k)/denom
#else
          u_bar(k) = ubr(k)/denom
#endif
c
          do j=1,jdim
            do i=1,idim
              up(i,j,k)= um(i,j,k) - u_bar(k)
c              if(k.le.kdist(i,j)) then
                u2_ave(k) = u2_ave(k)+up(i,j,k)**2
c              endif
            enddo
          enddo
c
          vb(k) = 0.0
          do j=1,jdim
            do i=1,idim
              vb(k)=vm(i,j,k)+vb(k)
            enddo
          enddo
#if defined(USEMPI)
          call MPI_ALLREDUCE(vb(k),v_bar8(k),1, 
     &         MPI_REAL8,MPI_SUM,
     &         comm, mpierr)
          v_bar(k) = v_bar8(k)/denom
#else
          v_bar(k) = vb(k)/denom
#endif
          do j=1,jdim
            do i=1,idim
              vp(i,j,k)= vm(i,j,k) - v_bar(k)
c              if(k.le.kdist(i,j)) then
                v2_ave(k) = v2_ave(k)+vp(i,j,k)**2
                w2_ave(k) = w2_ave(k)+wm(i,j,k)**2
c              endif
            enddo
          enddo
        enddo
c      endif
#else
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            up(i,j,k) = um(i,j,k)+uave
            vp(i,j,k) = vm(i,j,k)+vave
          enddo
        enddo
      enddo
#endif
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=0,idim1
c            wm(i,j,k) = wm(i,j,k)*way(i,j,k)
c          enddo
c        enddo
c      enddo

c
c first compute the fluxes for uv,uw, and vw
c
c add in ice motion uice
c
      do k=1,kdim1
        do j=1,jdim
          do i=1,idim
c
            odx = 2./(dx(i-1)+dx(i))
c
            et3(i,j,k) = -0.25*
     *                   ((um(i,j,k)-um(i,j-1,k))*ody*vay(i,j-1,k) -
     *                   (vm(i,j,k)-vm(i-1,j,k))*odx*vay(i-1,j,k)
     *                    - tomegaz)
            et3(i,j,k) = et3(i,j,k)*vay(i,j,k)
c
            et2(i,j,k) = 0.25*
     *                   ((um(i,j,k)-um(i,j,k-1))*odz*vay(i,j,k-1) -
     *                   (wm(i,j,k)-wm(i-1,j,k))*odx*vay(i-1,j,k)
     *                    + tomegay)
            et2(i,j,k) = et2(i,j,k)*vay(i,j,k)
c
            et1(i,j,k) = 0.25*((vm(i,j,k)-vm(i,j,k-1))*odz*vay(i,j,k-1)
     *                   -(wm(i,j,k)-wm(i,j-1,k))*ody*vay(i,j-1,k))
c
            et1(i,j,k) = et1(i,j,k)*vay(i,j,k)
c
c            kin(i,j,k) = 0.125*((um(i+1,j,k)+um(i,j,k))**2
cc     *                   + (vay(i,j+1,k)*vm(i,j+1,k)
cc     *                   +vay(i,j,k)*vm(i,j,k))**2 +
c     *                   + (vm(i,j+1,k)
c     *                   +vm(i,j,k))**2 +
c     *                   (wm(i,j,k+1)+wm(i,j,k))**2)
c     *                   +us(k)**2 + vs(k)**2
c            if(k.eq.kdist(i,j)) then
c              kin(i,j,k) = kin(i,j,k)*(1.-frac(i,j))
c            endif

c
c compute the subgrid scale terms
c
            kmo = 0.25*
     *        (kmom(i,j,k)+kmom(i-1,j,k)+kmom(i,j,k-1)+kmom(i-1,j,k-1))
c            kmo = kmo+kmb

            uws(i,j,k) = kmo*((umm(i,j,k)-umm(i,j,k-1))*odz*
     *                   vay(i,j,k)*vay(i,j,k-1) + 
     *                   vay(i,j,k)*vay(i-1,j,k)*
     *                   (wmm(i,j,k)-wmm(i-1,j,k))*odx)
            kmo = 0.25*
     *        (kmom(i,j,k)+kmom(i,j-1,k)+kmom(i,j,k-1)+kmom(i,j-1,k-1))
c            kmo = kmo+kmb

            vws(i,j,k) = kmo*((vmm(i,j,k)-vmm(i,j,k-1))*odz*
     *                   vay(i,j,k)*vay(i,j,k-1) + 
     *                   (wmm(i,j,k)-wmm(i,j-1,k))*ody*
     *                   vay(i,j,k)*vay(i,j-1,k))
            kmo =  0.25*
     *        (kmom(i,j,k)+kmom(i-1,j,k)+kmom(i,j-1,k)+kmom(i-1,j-1,k))
c            kmo = kmo+kmb
c
            uvs(i,j,k) =kmo*(vay(i-1,j,k)*(vmm(i,j,k)-vmm(i-1,j,k))*odx+
     *                   vay(i,j-1,k)*(umm(i,j,k)-umm(i,j-1,k))*ody)*
     *                   vay(i,j,k)
          enddo
        enddo
      enddo
      do j=1,jdim
        do i=1,idim
c          uws(i,j,kdim1) = 0.0
c          vws(i,j,kdim1) = 0.0
          kin(i,j,kdim1) = kin(i,j,kdim)
c          kin(i,j,0) = kin(i,j,1) - (kin(i,j,2) - kin(i,j,1))
          kin(i,j,0) = kin(i,j,1)
c
c upper b.c. for momentum wind stress input
c
#if !defined(RADTOP)
          if(strton) then
           dragu(i,j)=1.*str_tmask(i,j)*ustr*randn(i,j)/rhob(kdim)
           dragv(i,j)=1.*str_tmask(i,j)*vstr*randn(i,j)/rhob(kdim)
          else
            dragu(i,j) = ustr/rhob(kdim)
            dragv(i,j) = vstr/rhob(kdim)
          endif
#if defined(HBOT)
          uistart(i,j) = sqrt(abs(dragu(i,j))+abs(dragv(i,j)))
#else
          uistart(i,j) = sqrt(abs(dragu(i,j))+abs(dragv(i,j)))
          uistar(i,j) = sqrt(abs(dragu(i,j))+abs(dragv(i,j)))
          ustara = ustara + sqrt(abs(dragu(i,j))+abs(dragv(i,j)))
#endif
#endif

          uvs(i,j,kdim1) = uvs(i,j,kdim-1)
          uws(i,j,0) = 0.0
          uws(i,j,1) = 0.0
          vws(i,j,0) = 0.0
          uws(i,j,kdim1) = 0.0
          vws(i,j,kdim1) = 0.0
          vws(i,j,1) = 0.0
          et1(i,j,kdim1) = 0.0
          et2(i,j,kdim1) = 0.0
c #if !defined(RADIATE)
c          et1(i,j,1) = 0.0
c          et2(i,j,1) = 0.0
c          et3(i,j,1) = 0.0
c #endif
c
c
          tkew(i,j,kdim1) = 0.0
        enddo
      enddo
#if defined(HBOT)
c      vice = 0.
c      uice = 0.
c      write(lunit,*)"Ice Motion ",uice,vice
      deptho = 0.5*dz
      cd = 0.16/(log(0.5*dx(1)/0.006))**2
      do i=0,idim
        do j=0,jdim
          olev = kdist(i,j)
c          if(olev.lt.kdim) then
            do k=olev+1,kdim1
                  uws(i,j,k) = 0.0
                  uws(i+1,j,k) = 0.0
                  vws(i,j,k) = 0.0
                  vws(i,j+1,k) = 0.0
                  uvs(i,j,k) = 0.0
                  uvs(i+1,j,k) = 0.0
                  uvs(i,j+1,k) = 0.0
                  uvs(i+1,j+1,k) = 0.0
c                  um(i,j,k) = 0.0
c                  um(i+1,j,k) = 0.0
c                  umm(i,j,k) = 0.0
c                  umm(i+1,j,k) = 0.0
c                  vm(i,j,k) = 0.0
c                  vm(i,j+1,k) = 0.0
c                  vmm(i,j,k) = 0.0
c                  vmm(i,j+1,k) = 0.0
c                  wm(i,j,k) = 0.0
c                  wmm(i,j,k) = 0.0
                  et1(i,j,k) = 0.0
                  et2(i,j,k) = 0.0
                  kin(i,j,k) = 0.0
                  kmom(i,j,k) = 0.0
                  kheat(i,j,k) = 0.0
             enddo
             et2(i,j,kdist(i,j)) = 0.0
        enddo
      enddo
      cd = 0.16/(log(deptho/0.02))**2
c      cd = 0.0
      do j=1,jdim
        do i=1,idim
c          if(hice(i,j).gt.0.01) then
c            olev = kdist(i,j)
            olev = kdist(i,j)
            if(olev.lt.kdim-4) then
              draguh(i,j) = -cd*(um(i,j,olev)+uave)*
     *             abs(um(i,j,olev)+uave)
              dragvh(i,j) = -cd*(vm(i,j,olev)+vave)*
     *             abs(vm(i,j,olev)+vave)
c            dpea(1) = dpea(1) + dragv(i,j)
c            dragv(i,j) = 0.0
c            dragu(i,j) = 0.0
              uistarh(i,j) = sqrt(abs(dragu(i,j))+abs(dragv(i,j)))
c            ustara = ustara + uistar(i,j)
c          endif
           else
            dragvh(i,j) = 0.0
            draguh(i,j) = 0.0
            uistarh(i,j) = 0.0
           endif
        enddo
      enddo
      call imag2d(uistar)
      call imag2d(uistarh)
#endif
c

#if defined(BOTTOM)
      znot = 0.02
      do i=1,idim
        do j=1,jdim
          cd = 0.16/(log(0.5*dz/znot))**2
          
          dragub(i,j) = -1.*cd*(um(i,j,1)+uave-ugeo)*
     *          abs(um(i,j,1)+uave-ugeo)
          dragvb(i,j) = -1.*cd*(vm(i,j,1)+vave-vgeo(1))*
     *          abs(vm(i,j,1)+vave-vgeo(1))
          dpea(1) = dpea(1) + dragvb(i,j)
          uistar(i,j) = sqrt(abs(dragub(i,j))+abs(dragvb(i,j)))
          ustara = ustara + sqrt(abs(dragub(i,j))+abs(dragvb(i,j)))
        enddo
      enddo
#endif
#if defined(HBOTOLD)
      znot = 0.006
      do i=1,idim
        do j=1,jdim
          olev = kdist(i,j)
          deptho = 0.5*dz
          if(olev.gt.0) then
            do k=0,olev-1
                  uws(i,j,k) = 0.0
                  uws(i+1,j,k) = 0.0
                  vws(i,j,k) = 0.0
                  vws(i,j+1,k) = 0.0
                  uvs(i,j,k) = 0.0
                  uvs(i+1,j,k) = 0.0
                  uvs(i,j+1,k) = 0.0
                  uvs(i+1,j+1,k) = 0.0
                  um(i+1,j,k) = 0.0
                  um(i,j,k) = 0.0
                  umm(i,j,k) = 0.0
                  umm(i+1,j,k) = 0.0
                  vm(i,j,k) = 0.0
                  vm(i,j+1,k) = 0.0
                  vmm(i,j,k) = 0.0
                  vmm(i,j+1,k) = 0.0
                  wm(i,j,k) = 0.0
                  wmm(i,j,k) = 0.0
                  et1(i,j,k) = 0.0
                  et2(i,j,k) = 0.0
                  kmom(i,j,k) = 0.0
                  kheat(i,j,k) = 0.0
            enddo
            uws(i,j,olev) = 0.0
            uws(i+1,j,olev) = 0.0
            vws(i,j+1,olev) = 0.0
            vws(i,j,olev) = 0.0
            wm(i,j,olev) = 0.0
            wmm(i,j,olev) = 0.0
            et1(i,j,olev) = 0.0
            et2(i,j,olev) = 0.0
          endif
          cd = 0.16/(log(deptho/znot))**2
          dragub(i,j) = -1.*randn(i,j)*cd*(um(i,j,olev)+uave)*
     *             abs(um(i,j,olev)+uave)
          dragvb(i,j)=-1.*randn(i,j)*cd*(vm(i,j,olev)+vave)*
     *             abs(vm(i,j,olev)+uave)

c          dragvb(i,j) = 0.0
c          dragub(i,j) = 0.0
          dpea(1) = dpea(1) + dragvb(i,j)
          uistar(i,j) = sqrt(abs(dragub(i,j))+abs(dragvb(i,j)))
          ustara = ustara + sqrt(abs(dragub(i,j))+abs(dragvb(i,j)))
        enddo
      enddo
#endif
c
c
c set boundary image points
c
      call imag(et1)
      call imag(kin)
      call imag2d(uistar)
      call imag2d(uistart)
      call imag(uvs)
      call imag(uws)
      call imag(vws)
#if defined(OPENBC)
c      if(iblkstart(my_pe+1).eq.1) then
c      if(dxx*0.5/0.005 .gt. 1.0) then
c       cd = 0.16/(log(0.5*dxx/0.005))**2
c      else
c       cd = 0.0
c      endif
c
c  3-eq cd
c
       cd =.015
c       cd =0.0
c
c      write(*,*)"cd dxx",my_pe,dxx,cd
      do j=0,jdim1
        do i=1,idim
c          i=1
c          do k=1,kdist(i,j)
          do k=1,kdim
c
c look for boundary corners for inflo areas
c
             if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
               uws(i,j,k) = 0.0
               uvs(i,j,k) = 0.0
c               et3(i,j,k) = 0.0
c               et2(i,j,k) = 0.0
c               et1(i,j,k) = 0.0
               kin(0,j,k) = kin(1,j,k)
               dragw(j,k) = -cd*(wm(i,j,k))*
     *             abs(wm(i,j,k))
               dragvw(j,k) = -cd*(vm(i,j,k))*
     *             abs(vm(i,j,k))

             else
               dragw(j,k) = 0.0
               dragvw(j,k) = 0.0
             endif
            enddo
         enddo
      enddo
c      endif
#endif
      call imag(et3)
      call imag(et2)

c
c calculate PV
c
      do k=1,kdim
       do j=1,jdim
         do i=1,idim
c           dbz = ((+1)+rhop(i,j,k+1))-(rhob(k-1)+rhop(i,j,k-1)))
c     *            *0.5*odz
           dbz = -(th(i,j,k+1)-th(i,j,k-1))*0.5*odz*2e-4*rhob(k)
           dbx = (rhop(i+1,j,k)-rhop(i-1,j,k))*0.5*odz
     *           +0.0001*0.0002*rhob(k)
           dby = (rhop(i,j+1,k)-rhop(i,j-1,k))*0.5*odz
c
c           tke(i,j,k) = -9.81*(-(et1(i,j,k)+et1(i,j+1,k)+et1(i,j,k+1)
c     *       +et1(i,j+1,k+1))*dbx +
c     *       (et2(i,j,k)+et2(i+1,j,k)+et2(i,j,k+1)+et2(i+1,j,k+1))
c     *       *dby +
c     *       (et3(i,j,k)+et3(i+1,j,k)+et3(i,j+1,k)+et3(i+1,j+1,k))
c     *       *dbz)/rhob(k)+tke(i,j,k)
          enddo
        enddo
      enddo
      call imag(tke)

      do k=1,kdim
#if defined(ANALY)
c        a_ave(k) = 0.0
c        sd_ave(k) = 0.0
c        sg_ave(k) = 0.0
c        b_ave(k) = 0.0
c        uf_ave(k) = 0.0
c        vf_ave(k) = 0.0
c        wf_ave(k) = 0.0
        w_ave(k) = 0.0
        bw_ave(k) = 0.0
        sgw_ave(k) = 0.0
        sdw_ave(k) = 0.0
#endif
        do j=1,jdim
          do i=1,idim
       
c
            odx = 2./(dx(i-1)+dx(i))
c
c compute d(kinetic)/dx
c
c            uadvu = (kin(i,j,k)-kin(i-1,j,k))*odx
            uadvu = 0.0
c
c
c
            vadvu = (vm(i-1,j+1,k)+vm(i,j+1,k))*et3(i,j+1,k)+
     *              (vm(i-1,j,k)+vm(i,j,k))*et3(i,j,k)

            if(usetopref) then
              vadvu = vadvu + vave*tomegaz
            endif
c
c compute  ( d(wu)/dz + 2*wa*omega cos(lat) )
c

            if(k.ne.kdim1) then
              wadvu = (wm(i-1,j,k+1)+wm(i,j,k+1))*et2(i,j,k+1)+
     *              (wm(i-1,j,k)+wm(i,j,k))*et2(i,j,k)
            else
              wadvu = 0.0
            endif
c
c
c compute d/dx(Km d/dx(u)) and then other subgrid terms
c
c
            b1 = vay(i,j,k)*(vay(i+1,j,k)*(kmb+kmom(i,j,k))*
     *           (umm(i+1,j,k)-umm(i,j,k))/dx(i)-
     *           vay(i-1,j,k)*(kmb+kmom(i-1,j,k))*
     *           (umm(i,j,k)-umm(i-1,j,k))/dx(i-1))
c
            b2 = (uvs(i,j+1,k)-uvs(i,j,k))*ody
            if(k.ne.kdim1) then
              b3 = vay(i,j,k)*vay(i,j,k+1)*(uws(i,j,k+1)-uws(i,j,k))*odz
            else
              b3 = 0.0
            endif
c
#if defined(STOKES)
c
c compute stokes drift effects
c
c            sd = -0.5*vs(k)*((um(i,j+1,k)-um(i,j-1,k))*ody -
c     *           (vm(i,j+1,k)+vm(i,j,k)-vm(i-1,j+1,k)-vm(i-1,j,k))*odx)

            sd = 2.0*vs(k)*(et3(i,j+1,k)+et3(i,j,k))
#if defined(ICE)
            sd = sd*str_tmask(i,j)
#endif

#endif
c
c compute the tke dissipation rate
c
c            if(mod(n,ndiag).eq.0) then
             if(k.ne.kdim1) then
              disp(i,j,k) = (kmom(i,j,k)*(um(i+1,j,k)-um(i,j,k))*
     *                   (um(i+1,j,k)-um(i,j,k))/dx(i) +
     *                   kmom(i-1,j,k)*(um(i,j,k)-um(i-1,j,k))*
     *                   (um(i,j,k)-um(i-1,j,k))/dx(i-1)) +
     *                   uvs(i,j+1,k)*(um(i,j+1,k)-um(i,j,k))*tdy +
     *                   uvs(i,j,k)*(um(i,j,k)-um(i,j-1,k))*tdy +
     *                   uws(i,j,k+1)*(um(i,j,k+1)-um(i,j,k))*tdz +
     *                   uws(i,j,k)*(um(i,j,k)-um(i,j,k-1))*tdz
             endif
c            endif
c
c
            sgsru = 2.*b1+b2+b3
c
c
c compute the resolved tke budget term
c
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then
#if defined(STOKES)
              sd_ave(k) = sd_ave(k) + (up(i,j,k) * sd)
#endif
c              if(k.le.kdist(i,j)) then
                sg_ave(k) = sg_ave(k) + (sgsru * up(i,j,k))	  
c              endif
c            endif
#endif
c
c compute ru
c
            if(dofilt) then
              ru(i,j,k)=-vay(i,j,k)*filtu(i,j,k)/delt + uadvu - vadvu +
     *           wadvu - sgsru - sd
     *                    + raytau(k)*up(i,j,k)
c     *           - slopu*grav*rhop(i,j,k)/rhob(k)
#if defined(ANALY)
                 uf_ave(k) = uf_ave(k)+up(i,j,k)*filtu(i,j,k)/
     *                       delt
#endif
            else
              ru(i,j,k) = uadvu - vadvu + wadvu - sgsru - sd
     *                    + raytau(k)*up(i,j,k)
c     *                    - slopu*grav*rhop(i,j,k)/rhob(k)
            endif
#if defined(ANALY)
c            if(mod(n,ndiag).eq.0) then
#if !defined(RADIATE)
c              uf_ave(k) = uf_ave(k) + raytau(k)*up(i,j,k)**2
#endif
c              if(k.le.kdist(i,j)) then
                a_ave(k) = a_ave(k) - up(i,j,k)*(uadvu + wadvu -
     *          vadvu)
c              endif
c            endif
#endif
c
c compute ( d(uv)/dx + 2*ua*omega sine(lat) )
c
c
            uadvv = (um(i+1,j-1,k)+um(i+1,j,k))*et3(i+1,j,k)+
     *              (um(i,j-1,k)+um(i,j,k))*et3(i,j,k)
c
            if(usetopref) then
              uadvv = uadvv+uave*tomegaz
            endif
c
c compute v d/dy(v)
c
c            vadvv = (kin(i,j,k)-kin(i,j-1,k))*ody
            vadvv = 0.0
c
c compute  d(wv)/dz
c
            
            if(k.ne.kdim1) then
              wadvv = (wm(i,j-1,k+1)+wm(i,j,k+1))*et1(i,j,k+1) +
     *              (wm(i,j-1,k)+wm(i,j,k))*et1(i,j,k)
            else
              wadvv = 0.0
            endif

c
c compute d/dy(Km d/dy(v))
c
            b2 = (vay(i,j+1,k)*(kmb+kmom(i,j,k))*
     *           (vmm(i,j+1,k)-vmm(i,j,k))-
     *           vay(i,j-1,k)*(kmb+kmom(i,j-1,k))*
     *           (vmm(i,j,k)-vmm(i,j-1,k)))
     *           *ody*ody
c
            b1 = (uvs(i+1,j,k)-uvs(i,j,k))/dx(i)
            if(k.ne.kdim1) then
              b3 = (vws(i,j,k+1)-vws(i,j,k))*odz
            else
              b3 = 0.0
            endif
c
            sgsrv = b1+2.*b2+b3
c 
c compute the tke dissipation rate term
c 
c            if(mod(n,ndiag).eq.0) then
             if(k.ne.kdim1) then
              disp(i,j,k)=disp(i,j,k)+
     *                   (kmom(i,j,k)*(vmm(i,j+1,k)-vmm(i,j,k))*
     *                   (vm(i,j+1,k)-vm(i,j,k))+kmom(i,j-1,k)*
     *                   (vmm(i,j,k)-vmm(i,j-1,k))*
     *                   (vm(i,j,k)-vm(i,j-1,k)))*ody*ody +
     *                   uvs(i+1,j,k)*(vm(i+1,j,k)-vm(i,j,k))*odx*odx +
     *                   uvs(i,j,k)*(vm(i,j,k)-vm(i-1,j,k))*0.5*odx +
     *                   vws(i,j,k+1)*(vm(i,j,k+1)-vm(i,j,k))*tdz+
     *                   vws(i,j,k)*(vm(i,j,k)-vm(i,j,k-1))*tdz
             endif
c            endif
c

c
#if defined(STOKES)
c
c compute stokes drift effects
c
c            sd = 0.5*us(k)*((vm(i+1,j,k)-vm(i-1,j,k))*odx - 
c     *           (um(i+1,j,k)+um(i,j,k)-um(i+1,j-1,k)-um(i,j-1,k))*ody)

            sd = 2.0*us(k)*(et3(i+1,j,k)+et3(i,j,k))
#if defined(ICE)
	    sd = sd*str_tmask(i,j)
#endif
c
#endif
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then
#if defined(STOKES)
              sd_ave(k) = sd_ave(k) - (vp(i,j,k) * sd)
#endif
c              if(k.le.kdist(i,j))then
                sg_ave(k) = sg_ave(k) + (vp(i,j,k) * sgsrv)
c              endif
c            endif
#endif
c
c compute rv
c
            if(dofilt)then
              rv(i,j,k) = -ay(i,j,k)*filtv(i,j,k)/delt + uadvv + wadvv
     *           -sgsrv +vadvv + sd 
     *           + raytau(k)*vp(i,j,k)
#if defined(ANALY) 
              vf_ave(k) = vf_ave(k)+vp(i,j,k)*filtv(i,j,k)/ 
     *                       delt 
#endif 
            else
              rv(i,j,k) = uadvv + wadvv - sgsrv + vadvv + sd
     *                    + raytau(k)*vp(i,j,k)
            endif
c
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then  
#if !defined(RADIATE)
c              vf_ave(k) = vf_ave(k)+raytau(k)*vp(i,j,k)**2
#endif
c              if(k.le.kdist(i,j))then
                a_ave(k) = a_ave(k) -   
     *                   vp(i,j,k)*(uadvv
     *                    + vadvv + wadvv)
c              endif
#if defined(HBOT) 
c              if(k.ge.kdist(i,j)) then
                dpesg(k) = dpesg(k) - ay(i,j,k)*(vadvv + uadvv + wadvv)
                uf_ave(k) = uf_ave(k) - ay(i,j,k)*vadvv
                dpea(k) = dpea(k) + 
     *              ay(i,j,k)*(sgsrv + filtv(i,j,k)/delt)
c              endif
#endif
#if defined(ICE)
              if(k.le.kdist(i,j)) then
                dpesg(k) = dpesg(k) - vadvv - uadvv - wadvv
                uf_ave(k) = uf_ave(k) - vadvv
                dpea(k) = dpea(k) + sgsrv
              endif
#endif

c            endif
#endif
c
c
c compute ( d(uw)/dx - 2*ua*omega cos(lat) )
c
c
            uadvw = (um(i+1,j,k-1)+um(i+1,j,k))*et2(i+1,j,k)+
     *              (um(i,j,k-1)+um(i,j,k))*et2(i,j,k)
c
            if(usetopref) then
               uadvw = uadvw + uave*tomegay
            endif
c
c compute v d/dy(w)
c
           vadvw = (vm(i,j+1,k-1)+vm(i,j+1,k))*et1(i,j+1,k) +
     *              (vm(i,j,k-1)+vm(i,j,k))*et1(i,j,k)
c
c compute w d/dz(w)
c
c           if(k.eq.kdim1) then
c             wadvw = 0.0
c           else
#if defined(HBOT)
c             if(k.eq.kdist(i,j)) then
c              wadvw = (kin(i,j,k)-kin(i,j,k-1))*odz/(1.-frac(i,j))
c             else
c              wadvw = (kin(i,j,k)-kin(i,j,k-1))*odz
c             endif
#else
c            wadvw = (kin(i,j,k)-kin(i,j,k-1))*odz
#endif
c           endif
           wadvw = 0.0
c
c compute sgsrw
c
c compute d/dz(Km d/dz(w)) and other flux like terms
c
            if(k.eq.kdim1) then
              b3 = 0.0
            else
              b3u = vay(i,j,k+1)*vay(i,j,k)*
     *              (kmb+kmom(i,j,k))*(wmm(i,j,k+1)-wmm(i,j,k))
              b3d = vay(i,j,k)*vay(i,j,k-1)*
     *              (kmb+kmom(i,j,k-1))*(wmm(i,j,k)-wmm(i,j,k-1))
              b3 = (b3u-b3d)*odz*odz
            endif
c
            b1 = (uws(i+1,j,k)-uws(i,j,k))/dx(i)
            b2 = (vws(i,j+1,k)-vws(i,j,k))*ody
c
            sgsrw = b1+b2+2.*b3
c  
c compute the tke dissipation rate term 
c  
            tkew(i,j,k) =((kmb+kmom(i,j,k))*(wmm(i,j,k+1)-wmm(i,j,k))*
     *                   (wm(i,j,k+1)-wm(i,j,k))+(kmb+kmom(i,j,k-1))*
     *                   (wmm(i,j,k)-wmm(i,j,k-1))*
     *                   (wm(i,j,k)-wm(i,j,k-1)))*odz*odz +
     *                   uws(i+1,j,k)*(wm(i+1,j,k)-wm(i,j,k))*0.5*odx +
     *                   uws(i,j,k)*(wm(i,j,k)-wm(i-1,j,k))*0.5*odx +
     *                   vws(i,j+1,k)*(wm(i,j+1,k)-wm(i,j,k))*tdy+
     *                   vws(i,j,k)*(wm(i,j,k)-wm(i,j-1,k))*tdy
c 


c
#if defined(STOKES)
c
            if(k.ne.1) then
c              sd = (us(k)+us(k-1))*0.25*((um(i+1,j,k)+um(i,j,k) -
c     *           um(i+1,j,k-1)-um(i,j,k-1))*odz -
c     *           (wm(i+1,j,k)-wm(i-1,j,k))*odx)
c     *           +(vs(k-1)+vs(k))*0.25*((vm(i,j+1,k)+vm(i,j,k)-
c     *           vm(i,j+1,k-1)-vm(i,j,k-1))*odz -
c     *           (wm(i,j+1,k)-wm(i,j-1,k))*ody)
              sd = (us(k)+us(k-1))*(et2(i+1,j,k)+et2(i,j,k))+
     *             (vs(k-1)+vs(k))*(et1(i,j+1,k)+et1(i,j,k))
            else
              sd = 0.0
            endif
c
#if defined(ICE)
            sd = sd*str_tmask(i,j)
#endif
c
#endif
c
#if defined(ANALY)
c        if (mod(n,ndiag).eq.0) then
#if defined(STOKES)
          sdw_ave(k) = sdw_ave(k) + (wm(i,j,k) * sd)
#endif
          sgw_ave(k) = sgw_ave(k) + (sgsrw * wm(i,j,k))	  
c        endif
#endif
c
c compute rw

         rhot1 = rhobm1(k) + rhop(i,j,k)
         rhot2 = rhobm1(k-1) + rhop(i,j,k-1)
         rhot = (rhot1+rhot2)*0.5
c shut off sediment buoyancy effect
         rhot = rhot+(1.0-rhot/2650.)*tke(i,j,k)
         rhopt = rhot-0.5*(rhobm1(k)+rhobm1(k-1))
c         rhopt = (rhop(i,j,k)+ rhop(i,j,k-1))*0.5

            if(dofilt)then
              rw(i,j,k) = -vay(i+1,j,k)*filtw(i,j,k)/delt - uadvw+wadvw
     *           + grav*rhopt/rhob(k)
     *           - sgsrw - vadvw -sd
     *           + raytau(k)*wm(i,j,k)
c     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rhob(k)*0.5
#if defined(ICE)
#if !defined(WARMICE)
     *           + grav*cfr(i,j,k)*(920.-rho0)/rho0
#endif
#endif
#if defined(ANALY) 
              wf_ave(k) = wf_ave(k)+wm(i,j,k)*filtw(i,j,k)/ 
     *                       delt 
#endif 
            else
              rw(i,j,k) = -uadvw + wadvw
     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rhob(k)*0.5
     *           - sgsrw - vadvw -sd
     *           + raytau(k)*wm(i,j,k)
#if defined(ICE)
#if !defined(WARMICE)
     *           + grav*cfr(i,j,k)*(920.-rho0)/rho0
#endif
#endif

            endif
c
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then
#if !defined(RADIATE)
c              wf_ave(k) = wf_ave(k) + raytau(k)*wm(i,j,k)**2
#endif
              w_ave(k) = w_ave(k) - wm(i,j,k)*(-uadvw - vadvw + wadvw)
c              if(k.le.kdist(i,j)) then
              bw_ave(k) = bw_ave(k) - wm(i,j,k)*grav*
     *                   (rhop(i,j,k)+rhop(i,j,k-1))/rhob(k)*0.5
c              endif
#if defined(ICE)
#if !defined(WARMICE)
     *                   - wm(i,j,k)*grav*cfr(i,j,k)*(920.-rho0)/rho0
#endif
#endif
c            endif
#endif
          enddo
        enddo
      enddo
c 
c if the fields are larger than the cutoffs stop the program
c 
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
c
            if (abs(wm(i,j,k)).gt.maxw) then
               write (*,*)'w overflow = ',wm(i,j,k),th(i,j,k),sal(i,j,k)
               write (*,*) my_pe,'at i,j,k ',i,j,k
               write(*,*)i+iblkstart(my_pe+1),
     *               j+jblkstart(my_pe+1),k
               write (*,*)'th,sal boundary ',th(1,j,k),sal(1,j,k)
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return 
            end if
            if (abs(um(i,j,k)).gt.maxu) then
              write (*,*) 'u overflow = ',um(i,j,k),th(i,j,k),sal(i,j,k)
              write(*,*) 'rhop ',(rhop(i,j,kk),kk=1,k)
              write (*,*) my_pe,'at i,j,k ',i+iblkstart(my_pe+1),
     *               j+jblkstart(my_pe+1),k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return 
            end if
            if (abs(vm(i,j,k)).gt.maxv) then
              write (*,*) 'v overflow = ',vm(i,j,k),th(i,j,k),sal(i,j,k)
              write(*,*) 'rhop ',(rhop(i,j,kk),kk=1,k)
              write (*,*) my_pe,' at i,j,k ',i+iblkstart(my_pe+1),
     *               j+jblkstart(my_pe+1),k
              nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            end if
c
c test temperature, salinity range
c
            if(th(i,j,k).lt.thtamin) then
               write(*,*) 'theta less than thtamin on process',my_pe
               write(*,*) th(i,j,k),thtamin,vm(i,j,k),wm(i,j,k)
               write (*,*) 'at i,j,k ',i+iblkstart(my_pe+1),j
     *           +jblkstart(my_pe+1),k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            else if(th(i,j,k).gt.thtamax) then
               write(*,*) 'theta greater than thtamax on ',my_pe
               write(*,*) th(i,j,k),thtamax
c               write (*,*) 'at i,j,k ',i,j,k
               write (*,*) 'at i,j,k ',i+iblkstart(my_pe+1),j
     *           +jblkstart(my_pe+1),k,ax(i,j,k),th(i-1,j,k),
     *            ax(i+1,j,k),um(i,j,k),um(i+1,j,k),wm(i,j,k),
     *            wm(i,j,k+1),vm(i,j,k),vm(i,j+1,k)
               call slb_write(n,nerr)

               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            endif
            if(sal(i,j,k).lt.salmin) then
               write(*,*) 'sal less than salmin on',my_pe
               write(*,*) sal(i,j,k),salmin,vm(i,j,k),kdist(i,j)
c               write (*,*) 'at i,j,k ',i,j,k
              write (*,*) 'at i,j,k ',i+iblkstart(my_pe+1),j
     *           +jblkstart(my_pe+1),k,um(i,j,k),um(i-1,j,k),um(i+1,j,k)
              write(*,*) sal(i,j,k+1),sal(i+1,j,k),sal(i,j,k-1)

               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
c               return
            else if(sal(i,j,k).gt.salmax) then
               write(*,*) 'sal greater than salmax on',my_pe
               write(*,*)sal(i,j,k),salmax,um(i,j,k),um(i+1,j,k),
     *             kdist(i,j),sal(i,j,k-1),vm(i,j,k),wm(i,j,k+1),
     *             wm(i,j,k),sal(i,j,k+1),way(i,j,k),way(i,j,k+1),
     *             sal(i+1,j,k),vm(i,j+1,k),sal(i-1,j,k),ay(i,j,k),
     *             ay(i,j+1,k)
               write (*,*) 'at i,j,k ',i+iblkstart(my_pe+1),
     *             j+jblkstart(my_pe+1),k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
#if defined(OPENSPONGE2)
c
c sponge away vertical motion
c
                fmax = itotal-itotal/5
                if(iblkstart(my_pe+1)+i.gt. fmax) then
                  nsp = iblkstart(my_pe+1)+i-fmax
                  rw(i,j,k) = rw(i,j,k)
     *              +w(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                  ru(i,j,k) = ru(i,j,k)
     *              +up(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                  rv(i,j,k) = rv(i,j,k)
     *              +vp(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                endif

#endif

            endif
          enddo
        enddo
      enddo
c 
c
c do image points for ru and rw, then top and bot for rw
c
      k = kdim1
      do j=1,jdim
        do i=1,idim
c
c
c            uadvw = (um(i+1,j,k-1)+um(i+1,j,k))*et2(i+1,j,k) +
c     *              (um(i,j,k-1)+um(i,j,k))*et2(i,j,k)
            uadvw = 2.*(u_bar(k-1))*et2(i+1,j,k) +
     *              2.*(u_bar(k-1))*et2(i,j,k)
c
            if(usetopref) then
               uadvw = uadvw + uave*tomegay
            endif
c
c compute v d/dy(w)
c
           vadvw = (vm(i,j+1,k-1)+vm(i,j+1,k))*et1(i,j+1,k) +
     *              (vm(i,j,k-1)+vm(i,j,k))*et1(i,j,k)
           vadvw = 2.*(v_bar(k-1)*et1(i,j+1,k)+v_bar(k-1)*et1(i,j,k))
c
c compute w d/dz(w)
c
           wadvw = 0.0
c
c compute sgsrw
c
c compute d/dz(Km d/dz(w)) and other flux like terms
c
            odx = 1./dx(i)
c
            b1 = (uws(i+1,j,k)-uws(i,j,k))*odx
            b2 = (vws(i,j+1,k)-vws(i,j,k))*ody
c
            sgsrw = b1+b2
c
c stokes
c
#if defined(STOKES)
            sd = us(kdim1)*0.5*(um(i+1,j,kdim)+um(i,j,kdim) -
     *         um(i+1,j,kdim-1)-um(i,j,kdim-1))*odz
     *        +vs(kdim1)*0.5*(vm(i,j+1,kdim)+vm(i,j,kdim)-
     *         vm(i,j+1,kdim-1)-vm(i,j,kdim-1))*odz
#endif

c
c
c add du/dz term if needed later
c
c compute rw
c
            if(dofilt)then
              rw(i,j,k) = -filtw(i,j,k)/delt - uadvw + wadvw
     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rhob(k)*0.5
     *           - sgsrw - vadvw -sd
#if !defined(RADIATE)
c     *           + raytau(k)*wm(i,j,k)
#endif
#if defined(ANALY)
              wf_ave(k) = wf_ave(k)+wm(i,j,k)*filtw(i,j,k)/
     *                       delt
#endif
            else
              rw(i,j,k) = -uadvw + wadvw
     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rhob(k)*0.5
     *           - sgsrw - vadvw -sd
#if !defined(RADIATE)
c     *           + raytau(k)*wm(i,j,k)
#endif
            endif
c
c rw at the bottom
c
c
c         rw(i,j,kdim1) = rw(i,j,kdim)
         rw(i,j,1) = -0.5*(um(i,j,1)+um(i+1,j,1))*tomegay +
     *        grav*rhop(i,j,1)/rhob(1)
          
#if defined(HBOT)
         if(kdist(i,j).lt.kdim-4) then
           ru(i,j,kdist(i,j)) = ru(i,j,kdist(i,j))- draguh(i,j)*odz
           rv(i,j,kdist(i,j)) = rv(i,j,kdist(i,j))- dragvh(i,j)*odz
         endif
#endif
#if defined(BOTTOM)
         ru(i,j,1) = ru(i,j,1)- dragub(i,j)*odz
         rv(i,j,1) = rv(i,j,1)- dragvb(i,j)*odz
#endif
#if !defined(RADTOP) && !defined(ICE) 
         ru(i,j,kdim) = ru(i,j,kdim)- dragu(i,j)*odz
         rv(i,j,kdim) = rv(i,j,kdim)- dragv(i,j)*odz
#endif
#if defined(ICE)
         ru(i,j,kdist(i,j)) = ru(i,j,kdist(i,j))- dragu(i,j)*odz
         rv(i,j,kdist(i,j)) = rv(i,j,kdist(i,j))- dragv(i,j)*odz
#endif
        enddo
      enddo
c
c add ice edge roughness
c
#if defined(OPENBC)
c      if(iblkstart(my_pe+1).eq.1) then
        do i=1,idim
          do j=1,jdim
          do k=1,kdim
            if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
              rw(i,j,k) = rw(i,j,k)-dragw(j,k)*odz
              rv(i,j,k) = rv(i,j,k)-dragvw(j,k)*odz
              ru(i,j,k) = 0.0
            endif
          enddo
          enddo
        enddo
c      endif
#endif
c
#if defined(HBOT)
      do i=1,idim
        do j=1,jdim

          olev = kdist(i,j)
          if(olev.lt.kdim)then
            do k=olev+1,kdim1
              rw(i,j,k) = 0.0
              ru(i,j,k) = 0.0
              ru(i+1,j,k) = 0.0
              rv(i,j,k) = 0.0
              rv(i,j+1,k) = 0.0
            enddo
          endif
        enddo
      enddo
#endif
#if defined(HBOTOLD)
      do i=1,idim
        do j=1,jdim
          olev = kdist(i,j)
          if(olev.gt.0) then
            do k=0,olev-1
                  ru(i,j,k) = 0.0
                  ru(i+1,j,k) = 0.0
                  rv(i,j,k) = 0.0
                  rv(i,j+1,k) = 0.0
                  rw(i,j,k) = 0.0
            enddo
            rw(i,j,olev) = 0.0
          endif
        enddo
      enddo
#endif

c
#if !defined(USEMPI)
      do k=1,kdim
        do j=1,jdim
c 
c left boundary image point
c 
         ru(0,j,k) = ru(idim,j,k)
         rv(0,j,k) = rv(idim,j,k)
         rw(0,j,k) = rw(idim,j,k)
c 
c right boundary image point
c 
         ru(idim1,j,k) = ru(1,j,k)
         rv(idim1,j,k) = rv(1,j,k)
         rw(idim1,j,k) = rw(1,j,k)
c
        enddo
      enddo
      do j=1,jdim
        rw(0,j,kdim1)=rw(idim,j,kdim1)
        rw(idim1,j,kdim1)=rw(1,j,kdim1)
c

        if(nerr.eq.1) then
           write(*,*)'Error in uvwslb '
  	    return
        endif
      enddo
#endif

#if defined(ANALY)

      w_ave(kdim1) = 0.0
      sgw_ave(kdim1) = 0.0
      sdw_ave(kdim1) = 0.0
      bw_ave(kdim1) = 0.0
      do k=1,kdim
          a_ave(k) = a_ave(k) + 0.5*(w_ave(k+1)+w_ave(k))
          b_ave(k) = b_ave(k) + 0.5*(bw_ave(k+1)+bw_ave(k))
          sg_ave(k) = sg_ave(k) + 0.5*(sgw_ave(k+1)+sgw_ave(k))
          sd_ave(k) = sd_ave(k) + 0.5*(sdw_ave(k+1)+sdw_ave(k))
      enddo
c
c next compute the  total dissipation rate
c
c      if(mod(n,ndiag).eq.0) then
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              disp(i,j,k) = disp(i,j,k)+0.5*(tkew(i,j,k+1)+tkew(i,j,k))
c              tke(i,j,k) = disp(i,j,k)
              disp_ave(k) = disp_ave(k)-disp(i,j,k)
c              dcor(i,j,k) = disp(i,j,k)
            enddo
          enddo
        enddo
c      endif
#endif
      return
      end
