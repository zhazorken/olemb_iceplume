      subroutine uvwflux (dofilt,filtu,filtv,filtw,ru,rv,rw,n,nerr)
c
c  $Id: uvwflux.F,v 1.9 1999/09/08 22:34:06 eds Exp $
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
#if defined(ICE)
#include "ice.inc"
#endif
c 
c  process all the data slab by slab
c
c
      logical dofilt
      integer j,n,nerr,i,k
#if defined(USEMPI)
      integer mpierr
      real*8 u_bar8(kdim),v_bar8(kdim)
      real ub2dl(idim),vb2dl(idim),wb2dl(idim)
      real  ub2d(idim),vb2d(idim),wb2d(idim)

#endif
      real uadvu,wadvu,uadvw,wadvw,uadvv,wadvv,vadvu
      real ua,va,wa
      real filtu(0:idim1,0:jdim1,0:kdim1),filtw(0:idim1,0:jdim1,0:kdim1)
      real filtv(0:idim1,0:jdim1,0:kdim1)
c
      real ru(0:idim1,0:jdim1,0:kdim1),rv(0:idim1,0:jdim1,0:kdim1)
      real rw(0:idim1,0:jdim1,0:kdim1)
c
c temporary arrays for fluxes
c
      real uw(0:idim1,0:jdim1,0:kdim1),uv(0:idim1,0:jdim1,0:kdim1)
      real vw(0:idim1,0:jdim1,0:kdim1)
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
      real deptho
      integer olev
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
     &         MPI_COMM_WORLD, mpierr)
          u_bar(k) = u_bar8(k)/denom
#else     
          u_bar(k) = ubr(k)/denom
#endif    
c
          do j=1,jdim
            do i=1,idim
              up(i,j,k)= um(i,j,k) - u_bar(k)
              u2_ave(k) = u2_ave(k)+up(i,j,k)**2
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
     &         MPI_COMM_WORLD, mpierr)
          v_bar(k) = v_bar8(k)/denom
#else     
          v_bar(k) = vb(k)/denom
#endif    
          do j=1,jdim
            do i=1,idim
              vp(i,j,k)= vm(i,j,k) - v_bar(k)
              v2_ave(k) = v2_ave(k)+vp(i,j,k)**2
              w2_ave(k) = w2_ave(k)+wm(i,j,k)**2
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

      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            wm(i,j,k) = wm(i,j,k)*way(i,j,k)
          enddo
        enddo
      enddo

c
c first compute the fluxes for uv,uw, and vw
c
c add in ice motion uice
c
      do k=1,kdim1
        do j=1,jdim
          do i=1,idim
            uv(i,j,k) = 0.25*(um(i,j,k)+um(i,j-1,k))*
     *                  (vm(i,j,k)+vm(i-1,j,k))
            uw(i,j,k) = 0.25*(um(i,j,k)+um(i,j,k-1))*
     *                  (wm(i,j,k)+wm(i-1,j,k))
            vw(i,j,k) = 0.25*(vm(i,j,k)+vm(i,j,k-1))*
     *                  (wm(i,j,k)+wm(i,j-1,k))
c
c compute the subgrid scale terms
c
            kmo = 0.25*
     *        (kmom(i,j,k)+kmom(i-1,j,k)+kmom(i,j,k-1)+kmom(i-1,j,k-1))

            uws(i,j,k) = kmo*((um(i,j,k)-um(i,j,k-1))*odz + (wm(i,j,k)
     *                   -wm(i-1,j,k))*odx)
            kmo = 0.25*
     *        (kmom(i,j,k)+kmom(i,j-1,k)+kmom(i,j,k-1)+kmom(i,j-1,k-1))

            vws(i,j,k) = kmo*((vm(i,j,k)-vm(i,j,k-1))*odz + (wm(i,j,k)
     *                   -wm(i,j-1,k))*ody)
            kmo =  0.25*
     *        (kmom(i,j,k)+kmom(i-1,j,k)+kmom(i,j-1,k)+kmom(i-1,j-1,k))
c
            uvs(i,j,k) = kmo*((vm(i,j,k)-vm(i-1,j,k))*odx + (um(i,j,k)
     *                   -um(i,j-1,k))*ody)
          enddo
        enddo
      enddo
      do j=1,jdim
        do i=1,idim
c          uws(i,j,kdim1) = 0.0
c          vws(i,j,kdim1) = 0.0
c
c upper b.c. for momentum wind stress input
c
#if !defined(RADTOP)
          if(strton) then
            uw(i,j,kdim1)=-1.*str_tmask(i,j)*ustr*randn(i,j)/rhob(kdim)
            vw(i,j,kdim1)=-1.*str_tmask(i,j)*vstr*randn(i,j)/rhob(kdim)
          else
            uw(i,j,kdim1) = -1.*ustr*randn(i,j)/rhob(kdim)
            vw(i,j,kdim1) = -1.*vstr*randn(i,j)/rhob(kdim)
          endif
          uistar(i,j) = sqrt(abs(uw(i,j,kdim1))+abs(vw(i,j,kdim1)))
          ustara = ustara + uistar(i,j)
#else
c          uw(i,j,kdim1) = uw(i,j,kdim)
c          vw(i,j,kdim1) = vw(i,j,kdim)
#endif
          uvs(i,j,kdim1) = uvs(i,j,kdim)
          uws(i,j,1) = 0.0
          vws(i,j,1) = 0.0
          uws(i,j,0) = uws(i,j,1)
          vws(i,j,0) = uws(i,j,1)
          uw(i,j,0) = uw(i,j,1)
          vw(i,j,0) = vw(i,j,1)
c
c
          tkew(i,j,kdim1) = 0.0
        enddo
      enddo
#if defined(ICE)
      do i=1,idim
        do j=1,jdim
          kdist = int(hice(i,j)/dz)
          frac = hice(i,j)/dz - kdist
#if defined(WARMICE)
          if(frac.gt.0.8) then
#else
          if(frac.gt.0.4) then
#endif
            olev = kdim-kdist-1
            frac = 0.0
          else
            olev = kdim-kdist
          endif
#if defined(WARMICE)
          deptho = 0.5*dz + 0.5*frac*dz
#else
          deptho =  (kdim-olev)*dz + 0.5*dz - hice(i,j)
#endif
c
c try simple boundary box
c
          deptho = 0.5*dz
          olev = kdim-kdist

          if(olev.le.kdim) then
            do k=olev+1,kdim1
                  uw(i,j,k) = 0.0
                  vw(i,j,k) = 0.0
                  uws(i,j,k) = 0.0
                  vws(i,j,k) = 0.0
            enddo
          endif
          if(hice(i,j).gt.0.0) then
            cd = 0.16/(log(deptho/0.03))**2
                uw(i,j,olev+1) = randn(i,j)*cd*(um(i,j,olev)-uice+uave)*
     *             abs(um(i,j,olev)-uice+uave)
                vw(i,j,olev+1) = randn(i,j)*cd*(vm(i,j,olev)+vave)*
     *             abs(vm(i,j,olev)+vave)
          endif
          uistar(i,j) = sqrt(abs(uw(i,j,olev+1))+abs(vw(i,j,olev+1)))
        enddo
      enddo
#endif
#if defined(BOTTOM)
      znot = 0.03
      deptho = 0.5*dz
      do i=1,idim
        do j=1,jdim
          cd = 0.16/(log(deptho/znot))**2
              uw(i,j,1) = -1.*randn(i,j)*cd*(um(i,j,1)+uave)*
     *           abs(um(i,j,1)+uave)
              vw(i,j,1) = -1.*randn(i,j)*cd*(vm(i,j,1)+vave)*
     *           abs(vm(i,j,1)+vave)
        enddo
      enddo
#endif
#if defined(HBOT)
      znot = 0.03
      do i=1,idim
        do j=1,jdim
          kdist = int(hbot(i,j)/dz)
          frac = hbot(i,j)/dz - kdist
          if(frac.gt.0.5) then
            olev = kdist+2
            frac = 0.0
          else
            olev = kdist+1
          endif
          deptho = (1.-frac)*0.5*dz
          deptho = 0.5*dz
          if(olev.gt.0) then
            do k=0,olev-1
                  uw(i,j,k) = 0.0
                  vw(i,j,k) = 0.0
                  uws(i,j,k) = 0.0
                  vws(i,j,k) = 0.0
                  um(i,j,k) = 0.0
                  umm(i,j,k) = 0.0
                  vm(i,j,k) = 0.0
                  vmm(i,j,k) = 0.0
                  wm(i,j,k) = 0.0
                  wmm(i,j,k) = 0.0
                  wm(i,j,k+1) = 0.0
                  wmm(i,j,k+1) = 0.0
            enddo
          endif
          if(hbot(i,j).gt.0.0) then
            cd = 0.16/(log(deptho/znot))**2
                uw(i,j,olev) = -1.*randn(i,j)*cd*(um(i,j,olev)+uave)*
     *             abs(um(i,j,olev)+uave)
                vw(i,j,olev) = -1.*randn(i,j)*cd*(vm(i,j,olev)+vave)*
     *             abs(vm(i,j,olev)+vave)
          endif
c          uistar(i,j) = sqrt(abs(uw(i,j,olev))+abs(vw(i,j,olev)))
        enddo
      enddo
#endif
#if defined(OPENBC)
      cd = 0.16/(log(0.5*dz/0.06))**2
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
c
c look for boundary corners for inflo areas
c
                  if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c                  if(iblkstart(my_pe+1).eq.1 .and. i.le.1) then
                    spd = wm(i,j,k)**2+um(i,j,k)**2+vm(i,j,k)**2
                    uws(i,j,k) = 0.0
                    uvs(i,j,k) = 0.0
                    uw(i,j,k) = 0.0
                    vw(i,j,k) = 0.0
                    uv(i,j,k) = 0.0
                    dragw(j,k) = cd*spd
c                    dragw(j,k) = cd*(wm(i,j,k)*abs(wm(i,j,k)))
                  else
                    dragw(j,k) = 0.0
                  endif
            enddo
        enddo
      enddo
#endif

c
c
c set boundary image points
c
      call imag(uv)
      call imag(uw)
      call imag(vw)
      call imag(uvs)
      call imag(uws)
      call imag(vws)

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
c compute d(u^2)/dx
c
            uadvu =  ((um(i,j,k)+um(i+1,j,k))**2 - 
     *               (um(i,j,k)+um(i-1,j,k))**2)*0.25*odx
c
c compute ( d(vu)/dy - 2 v omega sin(lat) )
c
            va = 0.25*(vm(i,j+1,k)+vm(i-1,j+1,k)+vm(i,j,k)+vm(i-1,j,k))
c
#if defined(STOKESF)
            vadvu = (uv(i,j+1,k)-uv(i,j,k))*ody
     *           - va*tomegaz - vs(k)*tomegaz
c     *           - va*tomegaz - vs(k)*1.028e-4
#else
            vadvu = (uv(i,j+1,k)-uv(i,j,k))*ody
     *           - va*tomegaz
#endif
            if(usetopref) then
              vadvu = vadvu - vave*tomegaz
            endif
c
c compute  ( d(wu)/dz + 2*wa*omega cos(lat) )
c
            wa = (wm(i,j,k)+wm(i,j,k+1)+wm(i-1,j,k+1)+wm(i-1,j,k))*0.25

            wadvu = (uw(i,j,k+1)-uw(i,j,k))*odz+
     *            wa*tomegay
c
c
c compute d/dx(Km d/dx(u)) and then other subgrid terms
c
c
            b1 = (kmom(i,j,k)*(um(i+1,j,k)-um(i,j,k))-
     *           kmom(i-1,j,k)*(um(i,j,k)-um(i-1,j,k)))*odx*odx
c
            b2 = (uvs(i,j+1,k)-uvs(i,j,k))*ody
            b3 = (uws(i,j,k+1)-uws(i,j,k))*odz
c
#if defined(STOKES)
c
c compute stokes drift effects
c
            sd = -0.5*vs(k)*((um(i,j+1,k)-um(i,j-1,k))*ody -
     *           (vm(i,j+1,k)+vm(i,j,k)-vm(i-1,j+1,k)-vm(i-1,j,k))*odx)
#if defined(ICE)
            sd = sd*str_tmask(i,j)
#endif
 
#endif

c
c compute the tke dissipation rate
c
c            if(mod(n,ndiag).eq.0) then
              disp(i,j,k) = (kmom(i,j,k)*(um(i+1,j,k)-um(i,j,k))*
     *                   (um(i+1,j,k)-um(i,j,k)) +
     *                   kmom(i-1,j,k)*(um(i,j,k)-um(i-1,j,k))*
     *                   (um(i,j,k)-um(i-1,j,k)))*odx*odx +
     *                   uvs(i,j+1,k)*(um(i,j+1,k)-um(i,j,k))*tdy +
     *                   uvs(i,j,k)*(um(i,j,k)-um(i,j-1,k))*tdy +
     *                   uws(i,j,k+1)*(um(i,j,k+1)-um(i,j,k))*tdz +
     *                   uws(i,j,k)*(um(i,j,k)-um(i,j,k-1))*tdz
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
              sg_ave(k) = sg_ave(k) + (sgsru * up(i,j,k))	  
c            endif
#endif
c
c compute ru
c
            if(dofilt) then
              ru(i,j,k) = -filtu(i,j,k)/delt + uadvu + vadvu + 
     *           wadvu - sgsru - sd
#if defined(ANALY)
                 uf_ave(k) = uf_ave(k)+up(i,j,k)*filtu(i,j,k)/
     *                       delt
#endif
            else
              ru(i,j,k) = uadvu + vadvu + wadvu - sgsru - sd
            endif
c            uf_ave(k) = uf_ave(k)+wadvu
c            vf_ave(k) = vf_ave(k)+ru(i,j,k)
c            wf_ave(k) = wf_ave(k)+sgsru
#if defined(ANALY)
c            if(mod(n,ndiag).eq.0) then
c              sd_ave(k) = sd_ave(k)+uw(i,j,k)
              a_ave(k) = a_ave(k) - up(i,j,k)*(uadvu + wadvu+
     *        (uv(i,j+1,k)-uv(i,j,k))*ody)
c            endif
#endif
c
c compute ( d(uv)/dx + 2*ua*omega sine(lat) )
c
            ua = (um(i+1,j,k)+um(i+1,j-1,k)+um(i,j,k)+um(i,j-1,k))*0.25
c
#if defined(STOKESF)
            uadvv = (uv(i+1,j,k)-uv(i,j,k))*odx +
     *            ua*tomegaz+us(k)*tomegaz
c     *            ua*tomegaz+us(k)*1.028e-4
#else
            uadvv = (uv(i+1,j,k)-uv(i,j,k))*odx +
     *            ua*tomegaz
#endif
c
            if(usetopref) then
              uadvv = uadvv+uave*tomegaz
            endif
c
c compute v d/dy(v)
c
            vadvv = ((vm(i,j+1,k)+vm(i,j,k))**2 -
     *              (vm(i,j,k)+vm(i,j-1,k))**2)*ody*0.25
c
c compute  d(wv)/dz
c
            wadvv = (vw(i,j,k+1)-vw(i,j,k))*odz
c
c compute d/dy(Km d/dy(v))
c
            b2 = (kmom(i,j,k)*(vmm(i,j+1,k)-vmm(i,j,k))-
     *           kmom(i,j-1,k)*(vmm(i,j,k)-vmm(i,j-1,k)))*ody*ody
c
            b1 = (uvs(i+1,j,k)-uvs(i,j,k))*odx
            b3 = (vws(i,j,k+1)-vws(i,j,k))*odz
c
            sgsrv = b1+2.*b2+b3
c 
c compute the tke dissipation rate term
c 
c            if(mod(n,ndiag).eq.0) then
              disp(i,j,k)=disp(i,j,k)+
     *                   (kmom(i,j,k)*(vmm(i,j+1,k)-vmm(i,j,k))*
     *                   (vm(i,j+1,k)-vm(i,j,k))+kmom(i,j-1,k)*
     *                   (vmm(i,j,k)-vmm(i,j-1,k))*
     *                   (vm(i,j,k)-vm(i,j-1,k)))*ody*ody +
     *                   uvs(i+1,j,k)*(vm(i+1,j,k)-vm(i,j,k))*tdx +
     *                   uvs(i,j,k)*(vm(i,j,k)-vm(i-1,j,k))*tdx +
     *                   vws(i,j,k+1)*(vm(i,j,k+1)-vm(i,j,k))*tdz+
     *                   vws(i,j,k)*(vm(i,j,k)-vm(i,j,k-1))*tdz
c            endif
c

c
#if defined(STOKES)
c
c compute stokes drift effects
c
            sd = 0.5*us(k)*((vm(i+1,j,k)-vm(i-1,j,k))*odx - 
     *           (um(i+1,j,k)+um(i,j,k)-um(i+1,j-1,k)-um(i,j-1,k))*ody)
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
              sg_ave(k) = sg_ave(k) + (vp(i,j,k) * sgsrv)
c            endif
#endif
c
c compute rv
c
            if(dofilt)then
              rv(i,j,k) = -filtv(i,j,k)/delt + uadvv + wadvv - sgsrv
     *           + vadvv + sd
#if defined(ANALY) 
              vf_ave(k) = vf_ave(k)+vp(i,j,k)*filtv(i,j,k)/ 
     *                       delt 
#endif 
            else
              rv(i,j,k) = uadvv + wadvv - sgsrv + vadvv + sd
            endif
c
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then  
              a_ave(k) = a_ave(k) -   
     *                   (vp(i,j,k) * ((uv(i+1,j,k)-uv(i,j,k))*odx
     *                    + vadvv + wadvv))
c            endif
#endif
c
c
c compute ( d(uw)/dx - 2*ua*omega cos(lat) )
c

            ua = (um(i,j,k-1)+um(i,j,k)+um(i+1,j,k)+um(i+1,j,k-1))
     *           *0.25
c
#if defined(STOKESF)
            uadvw = (uw(i+1,j,k)-uw(i,j,k))*odx -
     *           ua*tomegay-(us(k)+us(k-1))*0.5*tomegay
#else
            uadvw = (uw(i+1,j,k)-uw(i,j,k))*odx -
     *           ua*tomegay
#endif
c
            if(usetopref) then
               uadvw = uadvw - uave*tomegay
            endif
c
c compute v d/dy(w)
c
            vadvw = (vw(i,j+1,k)-vw(i,j,k))*ody
c
c compute w d/dz(w)
c
            wadvw = ((wm(i,j,k+1)+wm(i,j,k))**2 - 
     *              (wm(i,j,k)+wm(i,j,k-1))**2)*0.25*odz
c
c compute sgsrw
c
c compute d/dz(Km d/dz(w)) and other flux like terms
c
            b3 = (kmom(i,j,k)*(wmm(i,j,k+1)-wmm(i,j,k))
     *           -kmom(i,j,k-1)*(wmm(i,j,k)-wmm(i,j,k-1)))*odz*odz
c
            b1 = (uws(i+1,j,k)-uws(i,j,k))*odx
            b2 = (vws(i,j+1,k)-vws(i,j,k))*ody
c
            sgsrw = b1+b2+2.*b3
c  
c compute the tke dissipation rate term 
c  
            tkew(i,j,k) =(kmom(i,j,k)*(wmm(i,j,k+1)-wmm(i,j,k))*
     *                   (wm(i,j,k+1)-wm(i,j,k))+kmom(i,j,k-1)*
     *                   (wmm(i,j,k)-wmm(i,j,k-1))*
     *                   (wm(i,j,k)-wm(i,j,k-1)))*odz*odz +
     *                   uws(i+1,j,k)*(wm(i+1,j,k)-wm(i,j,k))*tdx +
     *                   uws(i,j,k)*(wm(i,j,k)-wm(i-1,j,k))*tdx +
     *                   vws(i,j+1,k)*(wm(i,j+1,k)-wm(i,j,k))*tdy+
     *                   vws(i,j,k)*(wm(i,j,k)-wm(i,j-1,k))*tdy
c 


c
#if defined(STOKES)
c
            if(k.ne.1) then
              sd = (us(k)+us(k-1))*0.25*((um(i+1,j,k)+um(i,j,k) -
     *           um(i+1,j,k-1)-um(i,j,k-1))*odz -
     *           (wm(i+1,j,k)-wm(i-1,j,k))*odx)
     *           +(vs(k-1)+vs(k))*0.25*((vm(i,j+1,k)+vm(i,j,k)-
     *           vm(i,j+1,k-1)-vm(i,j,k-1))*odz -
     *           (wm(i,j+1,k)-wm(i,j-1,k))*ody)
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
c
            if(dofilt)then
              rw(i,j,k) = -filtw(i,j,k)/delt + uadvw + wadvw
     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rho0*0.5
     *           - sgsrw + vadvw -sd
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
              rw(i,j,k) = uadvw + wadvw
     *           + grav*(rhop(i,j,k)+rhop(i,j,k-1))/rho0*0.5
     *           - sgsrw + vadvw -sd
            endif
c
#if defined(ANALY)
c            if (mod(n,ndiag).eq.0) then
              w_ave(k) = w_ave(k) - wm(i,j,k)*(uadvw + vadvw + wadvw)
              bw_ave(k) = bw_ave(k) - wm(i,j,k)*grav*
     *                   (rhop(i,j,k)+rhop(i,j,k-1))/rho0*0.5
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
               write (*,*) 'w overflow = ',wm(i,j,k)
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return 
            end if
            if (abs(um(i,j,k)).gt.maxu) then
               write (*,*) 'u overflow = ',um(i,j,k),maxu
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return 
            end if
            if (abs(vm(i,j,k)).gt.maxv) then
               write (*,*) 'v overflow = ',vm(i,j,k)
               write (*,*) 'at i,j,k ',i,j,k
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
               write(*,*) th(i,j,k),thtamin
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            else if(th(i,j,k).gt.thtamax) then
               write(*,*) 'theta greater than thtamax on ',my_pe
               write(*,*) th(i,j,k),thtamax
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            endif
            if(sal(i,j,k).lt.salmin) then
               write(*,*) 'sal less than salmin'
               write(*,*) sal(i,j,k),salmin
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            else if(sal(i,j,k).gt.salmax) then
               write(*,*) 'sal greater than salmax'
               write(*,*) sal(i,j,k),salmax
               write (*,*) 'at i,j,k ',i,j,k
               nerr=1
#if defined(USEMPI)
               call MPI_ABORT(comm,1,mpierr)
#endif
               return
            endif
          enddo
        enddo
      enddo
c 
c
c do image points for ru and rw, then top and bot for rw
c
      do j=1,jdim
        do i=1,idim
c
c rw at the top.  w=0 is assumed, du/dz can be non-zero at top.
c
#if !defined(RADTOP)
         rw(i,j,kdim1) = -0.5*(um(i,j,kdim)+um(i+1,j,kdim))*tomegay +
     *        grav*rhop(i,j,kdim)/rho0 -
     *        2.*kmom(i,j,kdim)*wmm(i,j,kdim)*odz*odz
#if defined(STOKES)
     *        -us(kdim1)*0.5*(um(i+1,j,kdim)+um(i,j,kdim) -
     *         um(i+1,j,kdim-1)-um(i,j,kdim-1))*odz
     *        -vs(kdim1)*0.5*(vm(i,j+1,kdim)+vm(i,j,kdim)-
     *         vm(i,j+1,kdim-1)-vm(i,j,kdim-1))*odz
#if defined(STOKESF)
     *        -us(kdim1)*tomegay
#endif
#endif
#endif

c
c add du/dz term if needed later
c
         rw(i,j,1) = -0.5*(um(i,j,1)+um(i+1,j,1))*tomegay +
     *        grav*rhop(i,j,1)/rho0 
c     *        - 2.*kmom(i,j,1)*wmm(i,j,2)*odz*odz

        enddo
      enddo
#if defined(ICE)
      do i=1,idim
        do j=1,jdim
          kdist = int(hice(i,j)/dz)
          frac = hice(i,j)/dz - kdist
#if defined(WARMICE)
          if(frac.gt.0.8) then
#else
          if(frac.gt.0.4) then
#endif
            olev = kdim-kdist-1
          else
            olev = kdim-kdist
          endif
c simple box
          olev = kdim-kdist
          if(olev.lt.kdim)then
            do k=olev+1,kdim1
              rw(i,j,k) = 0.0
              ru(i,j,k) = 0.0
              rv(i,j,k) = 0.0
            enddo
          endif
        enddo
      enddo
#endif
#if defined(HBOT)
      do i=1,idim
        do j=1,jdim
          kdist = int(hbot(i,j)/dz)
          frac = hbot(i,j)/dz - kdist
          if(frac.gt.0.5) then
            olev = kdist+2
            frac = 0.0
          else
            olev = kdist+1
          endif
          if(olev.gt.0) then
            do k=0,olev-1
                  ru(i,j,k) = 0.0
                  rv(i,j,k) = 0.0
                  rw(i,j,k) = 0.0
            enddo
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
              disp_ave(k) = disp_ave(k)-disp(i,j,k)
            enddo
          enddo
        enddo
c      endif
#endif
      return
      end
