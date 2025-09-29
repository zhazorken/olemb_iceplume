       subroutine loop3d (n,nerr)
c 
c   $Id: loop3d.F,v 1.11 1999/09/20 22:48:22 dwd Exp $
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
c #if defined(ICE)
#include "ice.inc"
c #endif
#include "finitv.inc"
#include "communicate.inc"

#if defined(AVEFLUX)
#include "olemflx.inc"
#endif
c 
c     umm    um
c      +     +     +
c     t-1    t    t+1
c 
      integer j,n,nerr,i,k,itrc,jj,irelax,fmax,nsp
      real*8 drho,drhol
#if defined(USEMPI)
      integer mpierr
      real*8 uave8,vave8
#endif
      real*8 uaveloc,vaveloc
      real*8 aforc,aforcl
      real cputime
      real forc(0:idim1,0:jdim1,0:kdim1)
#if defined(cray)
      real wtavg(0:kdim1)
      real wmave(0:kdim1)
      real wtavgloc(0:kdim1)
      real wmavel(0:kdim1)
#else
      real*8 wtavg(0:kdim1)
      real*8 wmave(0:kdim1)
      real*8 wtavgloc(0:kdim1)
      real*8 wmavel(0:kdim1)
#endif
      real filtu(0:idim1,0:jdim1,0:kdim1),filtw(0:idim1,0:jdim1,0:kdim1)
      real filtv(0:idim1,0:jdim1,0:kdim1),rv(0:idim1,0:jdim1,0:kdim1)
      real ru(0:idim1,0:jdim1,0:kdim1),rw(0:idim1,0:jdim1,0:kdim1)
      real spd(0:idim1,0:jdim1,0:kdim1),spd1(0:idim1,0:jdim1,0:kdim1)
      real rum2(0:idim1,0:jdim1,0:kdim1),rwm2(0:idim1,0:jdim1,0:kdim1)
      real rum1(0:idim1,0:jdim1,0:kdim1),rwm1(0:idim1,0:jdim1,0:kdim1)
      real rvm1(0:idim1,0:jdim1,0:kdim1),rvm2(0:idim1,0:jdim1,0:kdim1)
c     real usk(0:idim1,0:jdim1,0:kdim1),vsk(0:idim1,0:jdim1,0:kdim1)
      real t_bar(0:jdim1,0:kdim1),tb
      real ua(kdim),va(kdim)
#if defined(ICE)
#if !defined(WARMICE)
      real cfrp(0:idim1,0:jdim1,0:kdim1)
#endif
#endif
      real xpgr,ypgr,zpgr
      real rut,rvt,rwt
      real dsadx,dsady,dtadx,dtady
      real fac,relax
      real znot,zzz,vland(0:kdim1)
      real lh,cd,m,kmm
      real aa,bb,lamb1,lamb2,f1,f2,tw,sw
      real cw,ci,latf,lamb3,gams,gamt,zdep
      real vedge,wedge,uedge
#if defined(ANALY)
      real wp_ave(kdim1)
      real*8 pel(kdim),zz
#endif
      real tranv
      real dischinit
c
      save rum1,rwm1,rvm1
      save rum2,rwm2,rvm2
      save uedge,vedge,wedge
      save dischinit
#if NEWSOLVE
      complex y1(idim,jdim,kdim)
      real gzb(idim,jdim),gzt(idim,jdim),ddx,ddz
      integer iimax,jjmax,kkmax,kk2max
#else
      real tmparr(10000)
#endif
c      real disch
c      save disch
c
      logical dofilt
c
c tmparr space for 128^3 = 1566
c tmpdir = 30+idim+kdim+2*jdim+max(idim,jdim,kdim)+7*(int((idim+1)/2)
c                 +int((kdim+1)/2))
c
      integer ierr,olev,ncl,nc
      real a(kdim),b(kdim),c(kdim),ft(idim,jdim,kdim)
c
#if defined(PRINTDEBUG)
      real pmax
#endif
c
#include "moddef.inc"
c 
c sediment discharge concentration
      sedinit = 0.0
      dischinit = 75.
c
c plume half width
      plumw = 105
      dphip = 2
      if(iope.and.n.eq.1) then
       write(*,*)"discharge velocity ",dischinit/(dphip*2*plumw*dxx*dxx)
       write(lunit,*)"plum half width and height ",plumw,dphip
      endif
c
c latent heat of fusion
c
c
      lh = 3.34e5/cp

      aa = rhob(1)*khb/(lh*920.);
      bb = rhob(1)*kms/920.;
      lamb1 = -5.73e-2;
c      lamb1 = -5.4e-2;
      lamb2 = 8.32e-2;
      lamb3 = 0.000761
      gamt = 0.022
      gams = 0.00062
      cw = 3974.
      ci = 2009.
      latf = 335000.


      cd = 0.16/(log(0.5/0.001))**2
      znot = 0.001
c      cd = 0.16/(log(0.5/znot))
      m = 0.054
c
c 3-eq cd
c
      cd = 0.0025
c      cd = 0.0075


      dsadx = 0.0/1000.
c      dsady = 0.0
c      dsadx = 0.0
c      dtadx = -0.5/5000.
      dtadx = 0.0
c      dtady = 0.0
#if defined(ICE)
c      vgeo = 0.0
      ugeo = 0.0
#else
c
c increase velocity over a 1.5 hour period
c
c      if(n*delt.lt.3600) then
c        vgeo = (n*delt)*0.16/3600.
c      else
c        vgeo = 0.16
c      endif
      ugeo = 0.0
#endif
#if defined(OPENBC)
      ugeo = 0.0
#endif
c       vgeo = 0.15
#if defined(FRONT2)
c      vgeo(0) = 0.0
c      vz = 0.0
c      vz = -19.62*0.0001
c      vz = +19.62*0.0001
c      vz = 2e-4*9.81*0.0001/tomegaz
c      vz = -2e-4*9.81*0.0004/tomegaz
      vz = -2e-4*9.81*0.0001/tomegaz
c      vz = -2e-4*9.81*0.000216/tomegaz
c      vz = -2e-4*9.81*0.00002/tomegaz
c      vz = 0.0
      vgeo(0) = vz*(-(kdim/2)*dz-0.5*dz)
      
      do k=1,kdim1
c        if(k.gt.kdim-6) then
c          vz = -2e-4*9.81*0.001/tomegaz
c          vgeo(k) = vgeo(k-1)+vz*dz
c        else
c          vz = -2e-4*9.81*0.0001/tomegaz
          vgeo(k) = vgeo(k-1)+vz*dz
c           vgeo(k) = 0.0
c        endif
      enddo
c      if(n.eq.1.and.iope) write(*,*)"vgeo ",(vgeo(k),k=1,kdim1)
#else
      do k=0,kdim1
c        vgeo(k) = 0.4
c        vgeo(k) = 0.188
c        vgeo(k) = 0.625
         vgeo(k) = 0.0
      enddo
      if(n.eq.1.and.iope) write(*,*)"vgeo ",(vgeo(k),k=1,kdim1)
#endif

#if !defined(HTOP) && !defined(ICE)
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            ax(i,j,k) = 1.0
            ay(i,j,k) = 1.0
            vay(i,j,k) = 1.0
            way(i,j,k) = 1.0
            bay(i,j,k) = 1.0
          enddo
        enddo
      enddo
      do j=0,jdim1
        do i=0,idim1
          kdist(i,j) = kdim
        enddo
      enddo
#endif
#if defined(SPONGE)
c      ustar = sqrt(0.007*vgeo*vgeo)
c      znot = 0.1
      do k=1,kdim
c        zzz = k*dz - 0.5*dz
c        vland(k) = 1.6*ustar*log(zzz/znot)
c        if(vland(k).gt.vgeo )then
c          vland(k) = vgeo
c        endif
         vland(k) = vgeo
      enddo
      vland(0) = vland(1)
      vland(kdim1) = vland(kdim)
#endif
      if(n.eq.1) then
c         call readUV('UVinit.nc',ua,va,dz)
c         if(iope) write(*,*)" velocity profiles"
c         if(iope) write(*,'(1x,2f7.2)'),(ua(k),va(k),k=1,kdim)
#if defined(EKMAN)
c      dek = 3.14159*sqrt(2*kmb/abs(tomegaz))
      dek = sqrt(2*kmb/abs(tomegaz))
      dek = 30.0
      vnot = 1.414*vstr/(dek*rhob(kdim)*tomegaz)
      if(iope) write(*,*)"Ekman depth and surface velocity ",dek,vnot
      zz = 0.5*dz
      do k=1,kdim
c       va(kdim-k+1)=vgeo(k)*(1.-cos(zz/dek)*exp(-zz/dek))
       ua(kdim-k+1)=vnot*cos(0.7854-zz/dek)*exp(-zz/dek)
c       ua(kdim-k+1)=-vgeo(k)*sin(zz/dek)*exp(-zz/dek)
       va(kdim-k+1)=vnot*sin(0.7854-zz/dek)*exp(-zz/dek)+vgeo(kdim-k+1)
c       va(kdim-k+1)=vnot*sin(0.7854-zz/dek)*exp(-zz/dek)
       zz = zz+dz
      enddo
      write(*,*)"ua ",(ua(k),k=1,kdim)
      write(*,*)"va ",(va(k),k=1,kdim)
#endif
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            rum1(i,j,k) = 0.0
            rvm1(i,j,k) = 0.0
            rwm1(i,j,k) = 0.0
            rum2(i,j,k) = 0.0
            rvm2(i,j,k) = 0.0
            rwm2(i,j,k) = 0.0
c               p(i,j,k) = 0.0
          enddo
        enddo
      enddo
      do k=1,kdim
        do j=0,jdim1
          do i=0,idim1
c              u(i,j,k) = ua(kdim-k+1)
c              um(i,j,k) = ua(kdim-k+1)
c              umm(i,j,k) = ua(kdim-k+1)
c              v(i,j,k) = va(kdim-k+1)
c              vm(i,j,k) = va(kdim-k+1)
c              vmm(i,j,k) = va(kdim-k+1)
#if defined(EKMAN)
              u(i,j,k) = ua(k)
              um(i,j,k) = ua(k)
              umm(i,j,k) = ua(k)
              v(i,j,k) = va(k)+vgeo(k)
              vm(i,j,k) = va(k)+vgeo(k)
              vmm(i,j,k) = va(k)+vgeo(k)
#else
              u(i,j,k) = 0.0
              um(i,j,k) = 0.0
              umm(i,j,k) = 0.0
#if defined(FRONT2)
               v(i,j,k) = vgeo(k)
               vm(i,j,k) = vgeo(k)
               vmm(i,j,k) = vgeo(k)
#else
c               v(i,j,k) = vgeo(k)
c               vm(i,j,k) = vgeo(k)
c               vmm(i,j,k) = vgeo(k)
              v(i,j,k) = 0.0
              vm(i,j,k) = 0.0
              vmm(i,j,k) = 0.0
#endif
#endif
        enddo
       enddo
      enddo
      if(iope) write(*,*)"call topbot"
      call topbotb(v,MDVVEL,n)
      call topbotb(u,MDUVEL,n)
      call topbotb(vm,MDVVEL,n)
      call topbotb(um,MDUVEL,n)
      call topbotb(vmm,MDVVEL,n)
      call topbotb(umm,MDUVEL,n)
      endif
c
c set vgeo to zero
c
c      do k=1,kdim
c        vgeo(k) = 0.0
c      enddo
      if(iope.and.n.eq.1) write(*,*)"vgeo 2 ",(vgeo(k),k=1,kdim)
c
c see if reference frame is shifted
c
      uaveloc = 0.0
      vaveloc = 0.0
      uave = 0.0
      vave = 0.0
      if(usetopref) then
        do j=1,jdim
          do i=1,idim
            uaveloc = uaveloc+um(i,j,kdim)
            vaveloc = vaveloc+vm(i,j,kdim)
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
c
        uave = uave8/dble(itotal*jtotal)
        vave = vave8/dble(itotal*jtotal)
#else
        uave = uaveloc/dble(itotal*jtotal)
        vave = vaveloc/dble(itotal*jtotal)
#endif
#if defined(WARMICE)
c        if(n*delt.lt.600) then
        if(n*delt.lt.1800) then
c        if(n*delt.lt.5400) then
c           uice = -0.25 - 0.2*(n*delt-3600)/1800
c          vice = -0.06*n*delt/1800.
           uice = -0.05*n*delt/1800.
c          uice = -0.01*n*delt/1800
c          uice = 0.17*n*delt/1800
c           uice = 0.0
        else
          uice = -0.05
c          uice = -0.01
c          vice = -0.45
c          vice = -0.1
c          uice = 0.17
c           vice = -0.06
c           uice = 0.0
        endif
c        uice = 0.0
        vice = 0.0
c        uice = -0.05
        uave = uice
        vave = vice
#endif
#if defined(HBOT)
c        vave = vgeo
#endif
#if defined(BOTTOM)
c        vave = vgeo
#endif
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              umm(i,j,k) = umm(i,j,k)-uave
              um(i,j,k) = um(i,j,k)-uave
              vmm(i,j,k) = vmm(i,j,k)-vave
              vm(i,j,k) = vm(i,j,k)-vave
            enddo
          enddo
        enddo
#if defined(ICE)
        do k=0,kdim
          do j=0,jdim
            do i=0,idim
              if(k.gt.kdist(i,j)) then
                umm(i,j,k) = 0.0
                um(i,j,k) = 0.0
                umm(i+1,j,k) = 0.0
                um(i+1,j,k) = 0.0
                vm(i,j,k) = 0.0
                vmm(i,j,k) = 0.0
                vm(i,j+1,k) = 0.0
                vmm(i,j+1,k) = 0.0
              endif
            enddo
          enddo
        enddo
        call imag(um)
        call imag(umm)
        call imag(vm)
        call imag(vmm)
#endif
      endif
c
c
#if defined(STOKES1)
c
c define stokes drift array using random perturbation
c maintain mass balance
c
c     call ranstok(usk,vsk)
             
      if(n.eq.1) then
         do k=0,kdim1
            do j=0,jdim1
               do i=0,idim1
                  um(i,j,k) = um(i,j,k)+us(k)*randn(i,j)
                  vm(i,j,k) = vm(i,j,k)+vs(k)*randn(i,j)
               enddo
            enddo
         enddo
      endif
#endif
c call filter first
c 
c 
      if (mod(n,filtskp).eq.0) then
c 
c filter the u and w fields, and store coefficients in filtu and filtw.
c 
         call filter(umm,filtu)
         call filter(wmm,filtw)
         call filter(vmm,filtv)
         dofilt=.true.
      else
         dofilt=.false.
c         do k=0,kdim1
c            do j=0,jdim1
c               do i=0,idim1
c                  filtu(i,j,k)=0.0
c                  filtw(i,j,k)=0.0
c                  filtv(i,j,k)=0.0
c               enddo
c            enddo
c         enddo
      end if
c      dofilt = .false.
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=1,idim1
c            filtw(i,j,k)=0.0
c          enddo
c        enddo
c      enddo
c 
c
c debugging stuff
c
c
	if(nerr.eq.1) then
		write(*,*)'Error before uvw in loop3d '
		return
	endif
      if(iope) write(*,*)"call uvw"
c
c
c call uvw to get the advection and diffusion of momentum
c
#if defined(UVWCENTER)
      call uvwcenter(dofilt,filtu,filtv,filtw,ru,rv,rw,n,nerr)
#endif
#if defined(UVWFLUX) 
      call uvwflux(dofilt,filtu,filtv,filtw,ru,rv,rw,n,nerr)

#endif
#if defined(UVWENST)
      call uvwenst(dofilt,filtu,filtv,filtw,ru,rv,rw,n,nerr)

#endif
      if(iope) write(*,*)"uvw done"

c
c need rv boundary updates to get pressure
c
      call imag(rv)
c
#if defined(HBOT)
c      do i=0,idim
c        do j=0,jdim
c          olev = kdist(i,j)
c          if(olev.gt.0) then
c            do k=0,olev-1
c                  rv(i,j,k) = 0.0
c                  rv(i,j+1,k) = 0.0
c            enddo
c          endif
c        enddo
c      enddo
#endif
c
c do the idim1 image point for ru
c
      call imag(ru)
c
      if(relstep.eq.1) then
       do k=0,kdim1
         do j=0,jdim1
            do i=0,idim1
              rum1(i,j,k) = ru(i,j,k)
              rvm1(i,j,k) = rv(i,j,k)
              rwm1(i,j,k) = rw(i,j,k)
              rum2(i,j,k) = ru(i,j,k)
              rvm2(i,j,k) = rv(i,j,k)
              rwm2(i,j,k) = rw(i,j,k)
            enddo
          enddo
       enddo
      endif
c      do k=1,kdim
c        do j=1,jdim
c          do i=1,idim
c            ru(i,j,k) = ru(i,j,k)-
c     *             kmb*((umm(i+1,j,k)+umm(i-1,j,k)-2.*umm(i,j,k))/
c     *             (dx(i)*dx(i))+(umm(i,j-1,k)+
c     *             umm(i,j+1,k)+umm(i,j,k+1)+umm(i,j,k-1)-4.*umm(i,j,k))
c     *             *tdxs)
c            rv(i,j,k) = rv(i,j,k)-
c     *             kmb*((vmm(i+1,j,k)+vmm(i-1,j,k)-2.*vmm(i,j,k))/
c     *             (dx(i)*dx(i))+(vmm(i,j-1,k)+
c     *             vmm(i,j+1,k)+vmm(i,j,k+1)+vmm(i,j,k-1)-4.*vmm(i,j,k))
c     *             *tdxs)
c            rw(i,j,k) = rw(i,j,k)-
c     *             kmb*((wmm(i+1,j,k)+wmm(i-1,j,k)-2.*wmm(i,j,k))/
c     *             (dx(i)*dx(i))+(wmm(i,j-1,k)+
c     *             wmm(i,j+1,k)+wmm(i,j,k+1)+wmm(i,j,k-1)-4.*wmm(i,j,k))
c     *             *tdxs)
c
c          enddo
c        enddo
c      enddo



c
c Adams-Bashforth
c
      do k=0,kdim1
         do j=0,jdim1
            do i=0,idim1
            rut = (23.*ru(i,j,k)-16.*rum1(i,j,k)+5.
     *            *rum2(i,j,k))/12.
            rvt = (23.*rv(i,j,k)-16.*rvm1(i,j,k)+5.
     *            *rvm2(i,j,k))/12.
            rwt = (23.*rw(i,j,k)-16.*rwm1(i,j,k)+5.
     *            *rwm2(i,j,k))/12.
c

            rum2(i,j,k) = rum1(i,j,k)
            rvm2(i,j,k) = rvm1(i,j,k)
            rwm2(i,j,k) = rwm1(i,j,k)
            rum1(i,j,k) = ru(i,j,k)
            rvm1(i,j,k) = rv(i,j,k)
            rwm1(i,j,k) = rw(i,j,k)
            ru(i,j,k) = rut*bax(i,j,k)
            rv(i,j,k) = rvt*vay(i,j,k)
            rw(i,j,k) = rwt*vay(i,j,k)
c
          enddo
        enddo
      enddo
c
c perform calculations that require 3-D operations
c
#if !defined(NOWTAVG)
#if !defined(HTOP)
c
      do k=0,kdim1
         wtavgloc(k) = 0.0
         wtavg(k) = 0.0
         wmavel(k) = 0.0
         wmave(k) = 0.0
         nc = 0
         ncl = 0
         do j=1,jdim
            do i=1,idim
#if defined(ICE)
               if(k.le.kdist(i,j)) then
#endif
#if defined(HBOT)
               if(k.gt.kdist(i,j)) then
#endif
                 wtavgloc(k) = wtavgloc(k) + rw(i,j,k)
                 wmavel(k) = wmavel(k) + wm(i,j,k)
                 ncl = ncl+1
#if defined(ICE)||defined(HBOT)
               endif
#endif
            enddo
         enddo
#if defined(USEMPI)
         call MPI_ALLREDUCE(wtavgloc(k),wtavg(k),1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
         call MPI_ALLREDUCE(wmavel(k),wmave(k),1,
     &        MPI_REAL8,MPI_SUM,comm, ierr)
         call MPI_ALLREDUCE(ncl,nc,1,
     &        MPI_INTEGER,MPI_SUM,comm, ierr)

#else
         wtavg(k) = wtavgloc(k)
         wmave(k) = wmavel(k)
         nc = ncl
#endif

         if(nc.gt.0) then
         wtavg(k) = wtavg(k)/dble(float(nc))
         wtavg(k) = 0.0
         wmave(k) = wmave(k)/dble(float(nc))
         endif
      enddo
      wtavg(kdim1) = 0.0

#else
      do k=0,kdim1
c       wtavg(k) = 0.0
c      wmave(k) = 0.0
      enddo
#endif
#if !defined(RADTOP)
      wtavg(kdim1) = 0.0
#endif
#if !defined(RADIATE)
      wtavg(1) = 0.0
#endif
      call deltime(cputime)
      veltim=veltim+cputime
c
c #if defined(HBOT)
c      wmave(0) = wmave(1)
       do j=0,jdim1
         do i=0,idim1
#if defined(ICE)
          do k=0,kdist(i,j)
#endif
#if defined(HBOT)
          do k=kdist(i,j),kdim1
#endif
#if !defined(HBOT)&&!defined(ICE)
          do k=0,kdim1
#endif
c            wm(i,j,k) = wm(i,j,k) - wmave(k)
c            wmm(i,j,k) = wmm(i,j,k) - wmave(k)
          enddo
        enddo
      enddo
c #endif
#else
      do k=0,kdim1
        wtavg(k) = 0.0
      enddo
#endif
      do k=0,kdim1
        wtavg(k) = 0.0
      enddo

c
c *******************************************************
c *******************************************************
c
c WARNING:  The following code for the pressure solver
c           assumes that dx=dy=dz.  ds=(dx=dy=dz)
c
c *******************************************************
c *******************************************************
c
      fmax = itotal-200
c
      if(iope) write(*,*)"compute velocities"

      do k=1,kdim
         do j=1,jdim
            do i=1,idim

c                xpgr=(p(i,j,k)-p(i-1,j,k))*odx
                xpgr = 0.0
                u(i,j,k) = umm(i,j,k)-delt*(ru(i,j,k)+xpgr)*bax(i,j,k)
     *                  -delt*vgeo(k)*tomegaz
c     *                  -0.2*delt*(hdx(i,j,k)+hdx(i-1,j,k))*umm(i,j,k)/40.
c
c                ypgr=(p(i,j,k)-p(i,j-1,k))*ody
                ypgr = 0.0

                v(i,j,k) = vmm(i,j,k)-delt*(rv(i,j,k)+ypgr)*ay(i,j,k)
c     *                  -0.2*delt*hdx(i,j,k)*vmm(i,j,k)/20.
c     *                  +delt*ugeo*tomegaz
c
c                zpgr=(p(i,j,k)-p(i,j,k-1))*odz
                zpgr = 0.0
                w(i,j,k)=wmm(i,j,k)-delt*(rw(i,j,k)
     *                  -wtavg(k)+zpgr)*vay(i,j,k)
c     *                  -0.2*delt*(hdx(i,j,k)+hdx(i,j,k-1))*wmm(i,j,k)/40.
#if defined(OPENSPONGE2)
c
c sponge away vertical motion
c
                if(iblkstart(my_pe+1)+i.gt. fmax) then
                  nsp = iblkstart(my_pe+1)+i-fmax
                  w(i,j,k) = w(i,j,k)
     *              -delt*w(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                endif
                if(iblkstart(my_pe+1)+i.gt. fmax) then
                  nsp = iblkstart(my_pe+1)+i-fmax
                  u(i,j,k) = u(i,j,k)
     *              -delt*up(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                endif
                if(iblkstart(my_pe+1)+i.gt. fmax) then
                  nsp = iblkstart(my_pe+1)+i-fmax
                  v(i,j,k) = v(i,j,k)
     *              -delt*vp(i,j,k)*(sin(1.5708*nsp/fmax)**2)/10.
                endif

#endif
c
c                dpres(i,j,k) = dpres(i,j,k)+(p(i+1,j,k)+p(i-1,j,k)+
c     *                  p(i,j+1,k)+p(i,j-1,k)+
c     *                  -4.*p(i,j,k))*tdxs

            enddo
         enddo
      enddo
      call topbotb(w,MDWVEL,n)
      call topbotb(v,MDVVEL,n)
      call topbotb(u,MDUVEL,n)
      call imag(v)
      call imag(u)
#if defined(ICE)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               if(k.gt.kdist(i,j))then
                 w(i,j,k) = 0.0
                 u(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 v(i,j+1,k) = 0.0
               endif
            enddo
         enddo
      enddo
      call imag(u)
      call imag(v)
#endif
#if defined(HTOP)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               u(i,j,k) = u(i,j,k)*vay(i,j,k)
               u(i+1,j,k) = u(i+1,j,k)*vay(i,j,k)
               v(i,j,k) = v(i,j,k)*vay(i,j,k)
               v(i,j+1,k) = v(i,j+1,k)*vay(i,j,k)
               w(i,j,k) = w(i,j,k)*vay(i,j,k)
               w(i,j,k+1) = w(i,j,k+1)*vay(i,j,k)
c               if(ay(i,j,k).eq.0.0)then
c                 w(i,j,k) = 0.0
c                 u(i,j,k) = 0.0
c                 u(i+1,j,k) = 0.0
c                 v(i,j,k) = 0.0
c                 v(i,j+1,k) = 0.0
c               endif
            enddo
         enddo
      enddo
      do j=1,jdim
        do i=1,idim
          w(i,j,kdim1) = 0.0
          w(i,j,0) = 0.0
        enddo
      enddo
  
#if defined(DISCHARGE)
c      disch = 150.
      do k=1,kdim
         do j=0,jdim
            do i=1,idim

           jj = jblkstart(my_pe+1)+j
           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c           if(bway(i,j,k).eq.1) then
           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
c           if(jj.gt.jtotal/2-50.and.jj.lt.jtotal/2+50) then
c           if(jj.gt.jtotal/2-30.and.jj.lt.jtotal/2+30) then
c           if(jj.gt.jtotal/2+200.and.jj.lt.jtotal/2+300) then

c           if(jj.gt.jtotal/2-30.and.jj.lt.jtotal/2+30) then
             if(k.lt.dphip+1) then
c             if(k.lt.5) then
c               u(i-1,j,k) = disch/(4*100.)
c               u(i,j,k) = disch/(4*60.*dxx*dxx)
c               um(i,j,k) = disch/(4*60.*dxx*dxx)
c               umm(i,j,k) = disch/(4*60.*dxx*dxx)
c               u(i,j,k) = disch/(4*100.*dxx*dxx)
c               um(i,j,k) = disch/(4*100.*dxx*dxx)
c               umm(i,j,k) = disch/(4*100.*dxx*dxx)
               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               ax(i,j,k) = 1.0
               ay(i,j,k) = 1.0
               ax(i+1,j,k) = 1.0
               ay(i,j+1,k) = 1.0
               th(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
c
c set sediment concentration in discharge water kg/m^3
c
               tke(i-1,j,k) = sedinit

             endif
           endif
           endif
           enddo
         enddo
      enddo
#endif

      call imag(u)
      call imag(th)
      call imag(sal)
      call imag(ax)
      call imag(v)
      call imag(w)
#endif
#if defined(HTOPO)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               if(hdx(i,j,k).gt.10.0)then
                 w(i,j,k) = 0.0
                 w(i,j,k+1) = 0.0
                 u(i,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 v(i,j,k) = 0.0
                 v(i,j+1,k) = 0.0
               endif
            enddo
         enddo
      enddo
      do j=1,jdim
        do i=1,idim
          w(i,j,kdim1) = 0.0
          w(i,j,0) = 0.0
        enddo
      enddo
      call imag(u)
      call imag(v)
#endif


#if defined(HBOT)
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
               if(k.gt.kdist(i,j))then
                 w(i,j,k) = 0.0
                 u(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 v(i,j+1,k) = 0.0
               endif
            enddo
         enddo
      enddo
      do j=1,jdim
        do i=1,idim
          w(i,j,kdist(i,j)+1) = 0.0
        enddo
      enddo
      call imag(u)
      call imag(v)
#endif
      call imag(w)
#if defined(OPENBC)
      call obcmom(rv,n)
#endif
c
c
      call topbotb(w,MDWVEL,n)
      call topbotb(v,MDVVEL,n)
      call topbotb(u,MDUVEL,n)

c
c next diffusion with 1/2 time step
c

#if defined(DIFFMELT)
      if(iope) write(*,*)"adi routine"
      call adi(u,n)
      call adi(v,n)
      call adi(w,n)

c      do ndif=1,4
c      call adi(u,n+ndif-1)
c      call adi(v,n+ndif-1)
c      call adi(w,n+ndif-1)
      call topbotb(w,MDWVEL,n)
      call topbotb(v,MDVVEL,n)
      call topbotb(u,MDUVEL,n)
      call obcmom(rv,n)
c     enddo
#endif
#if defined(CHANNEL)
        if(jblkstart(my_pe+1).eq.1) then
          do i=0,idim1
            do k=0,kdim1
              v(i,0,k) = 0.0
              v(i,1,k) = 0.0
              ay(i,0,k) = 0.0
              ay(i,1,k) = 0.0
            enddo
          enddo
        endif
        if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
          do i=0,idim1
            do k=0,kdim1
c              v(i,jdim,k) = 0.0
              v(i,jdim1,k) = 0.0
              ay(i,jdim1,k) = 0.0
            enddo
          enddo
        endif
#endif


      aforcl = 0.0
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
#if defined(CONRES)
c#if defined(HBOT) || defined(ICE)
               forc(i,j,k) = 
     +         - 1.*((ax(i+1,j,k)*u(i+1,j,k)-ax(i,j,k)*u(i,j,k))/dx(i)
     +         + ((ay(i,j+1,k)*v(i,j+1,k)-ay(i,j,k)*v(i,j,k))
     +         + (way(i,j,k+1)*w(i,j,k+1)-way(i,j,k)*w(i,j,k)))
     +            *ody)/delt
c
c make sure forcing is zero in discharge flow
c
c               forc(i,j,k) = forc(i,j,k)*bax(i+1,j,k)
c#else
c               forc(i,j,k) =
c     +       - 1.*((u(i+1,j,k)-u(i,j,k))/dx(i)
c     +       + ((v(i,j+1,k)-v(i,j,k))
c     +       + (way(i,j,k+1)*w(i,j,k+1)-way(i,j,k)*w(i,j,k)))*ody)/delt
c#endif
#else
               forc(i,j,k) = (ru(i+1,j,k)-ru(i,j,k)
     +              + ay(i,j+1,k)*rv(i,j+1,k)-ay(i,j,k)*rv(i,j,k)
     +              + rw(i,j,k+1)-rw(i,j,k)-wtavg(k+1)+wtavg(k))*ds
     +              - 0.67*((um(i+1,j,k)-um(i,j,k))
     +              + (vm(i,j+1,k)-vm(i,j,k))
     +              + (wm(i,j,k+1)-wm(i,j,k)))*ds/delt
#endif
            enddo
         enddo
      enddo
c
c compute normalization factor
c
c        call MPI_ALLREDUCE(aforcl,aforc,1,
c     &        MPI_REAL8,MPI_SUM,comm, ierr)
c
c        aforc = aforc/dble(itotal*jtotal*kdim)
c        if(iope) write(*,*)"aforc,aforcl = ",aforc,aforcl
c        if(iope) write(*,*)"u,v = ",u(5,5,5),v(5,5,5)
c        if(iope) write(*,*)"forc = ",forc(5,5,5),forc(5,8,5)
c
c #if defined(HBOT)
c      do i=1,idim
c        do j=1,jdim
c          do k=1,kdim
c          enddo
c        enddo
c      enddo
c#endif
        
#if defined(RADIATE)
c
c update lower pressure boundary before calling solver
c
c 
      call radbc(n)
c
#endif
#if defined(RADTOP)
c
c do the top boundary radiation condition
c
      call radbct(n)
#endif
c
c set up the pressure boundaries
c
      do j=0,jdim1
        do i=0,idim1
#if !defined(RADTOP) 
c           p(i,j,kdim1) =  p(i,j,kdim) - (rw(i,j,kdim1)-wtavg(kdim1))*dz
          p(i,j,kdim1) =  p(i,j,kdim)
c          p(i,j,kdim1) =  0.0
#endif

#if defined(HBOT)
c          olev = kdist(i,j)
c          do k=0,kdist(i,j)-1
c            p(i,j,k) = p(i,j,olev)
c          enddo
#endif
#if !defined(RADIATE)
c          p(i,j,0) = 0.0
c          p(i,j,0) = p(i,j,1) + (rw(i,j,1)-wtavg(1))*dz
          p(i,j,0) = p(i,j,1)
#endif

        enddo
      enddo

c      call topbotb(MDPRES,p,n)
      call imag(p)
#if defined(HTOP)

      do i=1,idim
      do j=0,jdim1
        do k=0,kdim1
          if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
            p(i-1,j,k) = p(i,j,k)
          endif
        enddo
      enddo
      enddo
#endif
      call imag(p)
c
#if defined(CONRES)
c
c reset background density
c
c         rho0=0.0
c         do k=1,kdim
c            drhol=0.0
c            do j=1,jdim
c               do i=1,idim
c                  s=sal(i,j,k)
c                  drhol=drhol+iestbl(s,th(i,j,k),5)
c               enddo
c            enddo
c            call MPI_ALLREDUCE(drhol,drho,1,
c     &        MPI_REAL8,MPI_SUM,comm, ierr)

c            rhobm1(k)=drho/(jtotal*itotal)
c            rhob(k)=1000.0+rhobm1(k)
c            if(iope) then
c              write(*,*)"k,rhob,t,s ",k,rhob(k),th(1,1,k),sal(1,1,k)
c            endif
c            rho0=rho0+rhob(k)
c         enddo
c         rhob(0)=rhob(1)
c         rhob(kdim1) = rhob(kdim)
c         rho0=rho0/float(kdim)
c         do j=1,jdim
c            do k=1,kdim
c               do i=1,idim
c                  s=sal(i,j,k)
cc                  rhop(i,j,k)=iestbl(s,th(i,j,k),2*k)-rhobm1(k)
c                  rhop(i,j,k)=iestbl(s,th(i,j,k),5)-rhobm1(k)
c               enddo
c            enddo
c            do i=1,idim
c               rhop(i,j,0)=rhop(i,j,1)
c               rhop(i,j,kdim1)=rhop(i,j,kdim)
c            enddo
c         enddo
c         call imag(rhop)

c
c compute pressure using conjugate residual
c
      call conres(n,forc,nerr)
c      if(iope) write(*,*)"sor "

c      call sorvdx(n,forc,nerr)
c      call imag(p)
#else

c
      do k=2,kdim-1
c
c coefficents for non-zero d/dy
c
         a(k) = 1.0
         b(k) =-2.0
         c(k) = 1.0
      enddo
      do k=2,kdim-1
         do j=1,jdim
            do i=1,idim
              ft(i,j,k)=-forc(i,j,k)
            enddo
         enddo
      enddo
c     
#ifdef NEWSOLVE
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            y1(i,j,k) = -cmplx(forc(i,j,k),0.0)
          enddo
        enddo
      enddo
      do j=1,jdim
        do i=1,idim
#if defined(RADIATE)
          gzb(i,j) = p(i,j,0)
#else
          gzb(i,j) = (rw(i,j,1)-wtavg(1))*dz
#endif
#if defined(RADTOP)
          gzt(i,j) = p(i,j,kdim1)
#else
          gzt(i,j) = -(rw(i,j,kdim1)-wtavg(kdim1))*dz
c          gzt(i,j) = 0.5*((abs(us(kdim1)+um(i,j,kdim))+
c     *               abs(vs(kdim1)+vm(i,j,kdim)))**2 +
c     *               vm(i,j,kdim)**2+ um(i,j,kdim)**2)
#endif
        enddo
      enddo
c
c
      ddx = dx
      ddz = dz
      iimax = idim
      jjmax = jdim
      kkmax = kdim
      kk2max = 2*kdim
      call poislk(y1,gzb,gzt,iimax,jjmax,kkmax,kk2max,ddx,ddz)
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            p(i,j,k) = real(y1(i,j,k))
          enddo
        enddo
      enddo
#else

      do j=1,jdim
         do i=1,idim
#if defined(RADIATE)
            ft(i,j,1)  = -forc(i,j,1)-p(i,j,0)
#else
            ft(i,j,1)  = -forc(i,j,1)-(rw(i,j,1)-wtavg(1))*dz
#endif
#if defined(RADTOP)
            ft(i,j,kdim)  = -forc(i,j,kdim)-p(i,j,kdim1)
#else
            ft(i,j,kdim) = -forc(i,j,kdim) +
     *           (rw(i,j,kdim1)-wtavg(kdim1))*dz
#endif
         enddo
      enddo
c
#if defined(RADIATE)
      b(1) = -2.0
#else
      b(1) = -1.0
#endif
      a(1) =  0.0
      c(1) = 1.0
c
#if defined(RADTOP)
       b(kdim) = -2.0
#else
       b(kdim) = -1.0
#endif

      a(kdim) = 1.0
      c(kdim) = 0.0
c
c call the 3-dimensional poission solver
c
      call pois3d(0,idim,1.0,0,jdim,1.0,1,kdim,a,b,c,
     *     idim,jdim,ft,ierr,tmparr)
c
      if(ierr.ne.0)then
         write(lunit,*)'pois3d error code = ',ierr
         nerr=1
         return
      endif
c
      do j=1,jdim
#if defined(PRINTDEBUG)
         pmax=0.0
#endif
         do k=1,kdim
            do i=1,idim
               p(i,j,k) = ft(i,j,k)
#if defined(PRINTDEBUG)
               pmax=amax1(pmax,abs(p(i,j,k)))
#endif
            enddo
         enddo
c
#if defined(PRINTDEBUG)
         write(*,*) 'pmax is',pmax,', for n=',n,', j=',j
#endif
c
      enddo
c
c endif for newsolve
c
#endif
c
c conres
c
#endif
c      call topbotb(p,MDPRES,n)
c      call imag(p)
c
      call deltime (cputime)
      prestim=prestim+cputime
c
c step u and w forward in time by 2*delt
c
#if defined(HTOP)

      do i=0,idim
      do j=0,jdim1
        do k=0,kdim1
          if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
            p(i-1,j,k) = p(i,j,k)
          endif
        enddo
      enddo
      enddo
#endif

#if defined(RADTOP)
      do k=1,kdim1
#else
      do k=1,kdim
#endif
#if defined(ANALY)
         wp_ave(k) = 0.0
c         p_ave(k) = 0.0
#endif
         do j=1,jdim
c
c set sponge coefficient to 1/8 of domain size with a 5 minute
c relaxation time period
c
#if defined(SPONGE)
             ii = iblkstart(my_pe+1)+i
             irelax = jtotal/8
             if(ii.lt.irelax/2) then
               relax = delt/(20.*(irelax/2-ii+1))
             else if(ii.gt.irelax/2.and.ii.lt.irelax) then
               relax = delt/(20.*(ii-irelax/2+1))
             else
               relax = 0.0
             endif
#endif

            do i=1,idim
                xpgr=2.*(p(i,j,k)-p(i-1,j,k))/(dx(i)+dx(i-1))
                u(i,j,k) = u(i,j,k)-delt*xpgr*bax(i,j,k)

c     
                ypgr=(p(i,j,k)-p(i,j-1,k))*ody
                v(i,j,k) = v(i,j,k)-delt*ypgr*vay(i,j,k)
c
                zpgr=(p(i,j,k)-p(i,j,k-1))*odz*vay(i,j,k)
                w(i,j,k)=w(i,j,k)-delt*zpgr
c            if(my_pe.eq.4) then
c              if(i.eq.4.and.j.eq.8.and.k.eq.236) then
c                write(*,*)"ru, xpgr, ypgr, zpgr " ,
c     *                ru(i,j,k), xpgr,ypgr,zpgr
c              endif
c            endif


#if defined(ANALY)
c               if(mod(n,ndiag).eq.0) then
c                  p_ave(k) = p_ave(k) - (um(i,j,k)+uave)*xpgr
c     *                       - (vm(i,j,k)+vave)*ypgr
c                  p_ave(k) = p_ave(k) - (um(i,j,k)-u_b2d(i,k))*xpgr
c     *                       - (vm(i,j,k)-v_b2d(i,k))*ypgr
                  p_ave(k) = p_ave(k) - (um(i,j,k)-u_bar(k))*xpgr
     *                       - (vm(i,j,k)-v_bar(k))*ypgr
                  wp_ave(k) =  wp_ave(k)-wm(i,j,k)
     *                         *(zpgr-wtavg(k))
c
c                  if(k.ge.kdist(i,j)) then
c                  wf_ave(k) = wf_ave(k)-ypgr
c                  endif
c               
c               endif

#endif
c #if defined(ICE) || defined(HBOT)
c               else
c                u(i,j,k) = 0.0
c                v(i,j,k) = 0.0
c                w(i,j,k) = 0.0
c               endif
c #endif
#if defined(RADIATE)
              if(k.le.8) then
              fac = 0.1/float(k)
              u(i,j,k) = u(i,j,k) + fac*(umm(i+1,j,k)+umm(i-1,j,k)-
     *           4.*umm(i,j,k)+umm(i,j+1,k)+umm(i,j-1,k)) 
              v(i,j,k) = v(i,j,k) + fac*(vmm(i+1,j,k)+vmm(i-1,j,k)-
     *           4.*vmm(i,j,k)+vmm(i,j+1,k)+vmm(i,j-1,k))
              w(i,j,k) = w(i,j,k) + fac*(wmm(i+1,j,k)+wmm(i-1,j,k)-
     *           4.*wmm(i,j,k)+wmm(i,j+1,k)+wmm(i,j-1,k))
              endif
#endif
#if defined(RADTOP)
              if(k.ge.kdim-8) then
              fac = 0.12/float(kdim1+1-k)
              u(i,j,k) = u(i,j,k) + fac*(umm(i+1,j,k)+umm(i-1,j,k)-
     *           4.*umm(i,j,k)+umm(i,j+1,k)+umm(i,j-1,k))
              v(i,j,k) = v(i,j,k) + fac*(vmm(i+1,j,k)+vmm(i-1,j,k)-
     *           4.*vmm(i,j,k)+vmm(i,j+1,k)+vmm(i,j-1,k))
              w(i,j,k) = w(i,j,k) + fac*(wmm(i+1,j,k)+wmm(i-1,j,k)-
     *           4.*wmm(i,j,k)+wmm(i,j+1,k)+wmm(i,j-1,k))
              endif
#endif
            enddo
         enddo
c
      enddo
#if defined(ANALY)
c      if(mod(n,ndiag).eq.0) then
        wp_ave(kdim1) = 0.0
        do k=1,kdim
          p_ave(k) = p_ave(k) + 0.5*(wp_ave(k)+wp_ave(k+1))
          dpesg(k) = dpesg(k) + 0.5*(wp_ave(k)+wp_ave(k+1))
        enddo
c      endif
#endif
c
      call topbotb(w,MDWVEL,n)
      call topbotb(v,MDVVEL,n)
      call topbotb(u,MDUVEL,n)
      call imag(u)
      call imag(v)
      call imag(w)
#if defined(OPENBC)
      call obcmom(rv,n)
#endif
c
c flux output routines
c
#if defined(AVEFLUX)
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
c
           uaa(i,j,k) = uaa(i,j,k)+u(i,j,k)
           vaa(i,j,k) = vaa(i,j,k)+v(i,j,k)
           waa(i,j,k) = waa(i,j,k)+w(i,j,k)
           taa(i,j,k) = taa(i,j,k)+th(i,j,k)
           saa(i,j,k) = saa(i,j,k)+sal(i,j,k)
c
           uc = u(i,j,k)+u(i+1,j,k)
           vc = v(i,j,k)+v(i,j+1,k)
           wc = w(i,j,k)+w(i,j,k+1)
c
           uwa(i,j,k) = uwa(i,j,k)+0.25*uc*wc
           vwa(i,j,k) = vwa(i,j,k)+0.25*vc*wc
           uva(i,j,k) = uva(i,j,k)+0.25*uc*vc
c
           uta(i,j,k) = uta(i,j,k)+uc*th(i,j,k)*0.5
           vta(i,j,k) = vta(i,j,k)+vc*th(i,j,k)*0.5
           wta(i,j,k) = wta(i,j,k)+wc*th(i,j,k)*0.5
c
           usa(i,j,k) = usa(i,j,k)+uc*sal(i,j,k)*0.5
           vsa(i,j,k) = vsa(i,j,k)+vc*sal(i,j,k)*0.5
           wsa(i,j,k) = wsa(i,j,k)+wc*sal(i,j,k)*0.5

          enddo
        enddo
      enddo
#endif
#if defined(HTOP)
      do k=0,kdim1
         do j=0,jdim1
            do i=0,idim1
               u(i,j,k) = u(i,j,k)*vay(i,j,k)
               v(i,j,k) = v(i,j,k)*vay(i,j,k)
               w(i,j,k) = w(i,j,k)*vay(i,j,k)
            enddo
         enddo
      enddo
#endif

#if defined(CHANNEL)
        if(jblkstart(my_pe+1).eq.1) then
          do i=0,idim1
            do k=0,kdim1
              v(i,0,k) = 0.0
              v(i,1,k) = 0.0
            enddo
          enddo
        endif
        if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
          do i=0,idim1
            do k=0,kdim1
c              v(i,jdim,k) = 0.0
              v(i,jdim1,k) = 0.0
            enddo
          enddo
        endif
#endif

#if defined(DISCHARGE)
           jj = jblkstart(my_pe+1)+j
           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
               u(i,j,k) = disch/(dphip*2.*plumw*dxx*dxx)
               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
             endif
           endif
           endif
#endif

c
c     umm    um    u
c      +     +     +
c      t    t    t+1
c 
c copy u to um
c and compute time average velocity field in u,v, and w 
c
      do k=0,kdim
         do j=0,jdim1
            do i=0,idim1
c
#if defined(STOKES)
               um(i,j,k)=u(i,j,k)
               vm(i,j,k)=v(i,j,k)
               u(i,j,k) = (um(i,j,k)+umm(i,j,k))*0.5
#if defined(ICE)
     *             +us(k)*str_tmask(i,j)
#else
c     *             +us(k)
#endif
               v(i,j,k) = (vm(i,j,k)+vmm(i,j,k))*0.5
c     *             +vs(k)
#else
               um(i,j,k)=u(i,j,k)
               vm(i,j,k)=v(i,j,k)
               u(i,j,k) = (um(i,j,k)+umm(i,j,k))*0.5
               v(i,j,k) = (vm(i,j,k)+vmm(i,j,k))*0.5
#endif
               wm(i,j,k)=w(i,j,k)
               w(i,j,k) = (wm(i,j,k)+wmm(i,j,k))*0.5
c
               u(i,j,k) = u(i,j,k)*ax(i,j,k)
               v(i,j,k) = v(i,j,k)*ay(i,j,k)
               w(i,j,k) = w(i,j,k)*way(i,j,k)
             enddo
           enddo
        enddo
#if defined(CHANNEL)
        if(jblkstart(my_pe+1).eq.1) then
          do i=0,idim1
            do k=0,kdim1
              v(i,0,k) = 0.0
              v(i,1,k) = 0.0
            enddo
          enddo
        endif
        if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
          do i=0,idim1
            do k=0,kdim1
c              v(i,jdim,k) = 0.0
              v(i,jdim1,k) = 0.0
            enddo
          enddo
        endif
#endif

        do k=0,kdim1
          do j=1,jdim
            do i=1,idim
#if defined(HTOPO)
               if(hdx(i,j,k).gt.10.0)then
                 w(i,j,k+1) = 0.0
                 wm(i,j,k+1) = 0.0
                 wmm(i,j,k+1) = 0.0
                 w(i,j,k) = 0.0
                 wm(i,j,k) = 0.0
                 wmm(i,j,k) = 0.0
                 u(i,j,k) = 0.0
                 um(i,j,k) = 0.0
                 umm(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 vm(i,j,k) = 0.0
                 vmm(i,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 um(i+1,j,k) = 0.0
                 umm(i+1,j,k) = 0.0
                 v(i,j+1,k) = 0.0
                 vm(i,j+1,k) = 0.0
                 vmm(i,j+1,k) = 0.0
               endif
#endif
#if defined(HTOP)
               if(ax(i,j,k).eq.0.0.and.ax(i+1,j,k).eq.0.0)then
                 w(i,j,k) = 0.0
                 wm(i,j,k) = 0.0
                 wmm(i,j,k) = 0.0
                 w(i,j,k+1) = 0.0
                 wm(i,j,k+1) = 0.0
                 wmm(i,j,k+1) = 0.0
                 u(i,j,k) = 0.0
                 um(i,j,k) = 0.0
                 umm(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 vm(i,j,k) = 0.0
                 vmm(i,j,k) = 0.0
c                 u(i+1,j,k) = 0.0
c                 um(i+1,j,k) = 0.0
c                 umm(i+1,j,k) = 0.0
c                 v(i,j+1,k) = 0.0
c                 vm(i,j+1,k) = 0.0
c                 vmm(i,j+1,k) = 0.0
               endif
#if defined(DISCHARGE)
           jj = jblkstart(my_pe+1)+j
           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
               u(i,j,k) = disch/(dphip*2.*plumw*dxx*dxx)
               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)

               th(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
               tke(i-1,j,k) = sedinit
             endif
           endif
           endif
#endif

#endif

#if defined(HBOT)
               if(k.gt.kdist(i,j).and.k.ne.kdim1)then
                 w(i,j,k+1) = 0.0
                 wm(i,j,k+1) = 0.0
                 wmm(i,j,k+1) = 0.0
                 w(i,j,k) = 0.0
                 wm(i,j,k) = 0.0
                 wmm(i,j,k) = 0.0
                 u(i,j,k) = 0.0
                 um(i,j,k) = 0.0
                 umm(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 vm(i,j,k) = 0.0
                 vmm(i,j,k) = 0.0
c                 um(i+1,j,k) = 0.0
c                 umm(i+1,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 um(i+1,j,k) = 0.0
                 umm(i+1,j,k) = 0.0
                 v(i,j+1,k) = 0.0
                 vm(i,j+1,k) = 0.0
                 vmm(i,j+1,k) = 0.0
               endif
               w(i,j,kdist(i,j)+1) = 0.0
#endif
#if defined(ICE)
               if(k.gt.kdist(i,j))then
                 w(i,j,k) = 0.0
                 wm(i,j,k) = 0.0
                 wmm(i,j,k) = 0.0
                 u(i,j,k) = 0.0
                 um(i,j,k) = 0.0
                 umm(i,j,k) = 0.0
                 v(i,j,k) = 0.0
                 vm(i,j,k) = 0.0
                 vmm(i,j,k) = 0.0
                 um(i+1,j,k) = 0.0
                 umm(i+1,j,k) = 0.0
                 u(i+1,j,k) = 0.0
                 v(i,j+1,k) = 0.0
                 vm(i,j+1,k) = 0.0
                 vmm(i,j+1,k) = 0.0
               endif
#endif

            enddo
         enddo
      enddo
      do j=0,jdim1
        do i=0,idim1
          way(i,j,kdim1) = 0.0
          w(i,j,kdim1) = 0.0
          wm(i,j,kdim1) = 0.0
          wmm(i,j,kdim1) = 0.0
        enddo
      enddo
c     
      do k=0,kdim1
         do j=0,jdim1
            do i=0,idim1
               umm(i,j,k) = um(i,j,k)
               vmm(i,j,k) = vm(i,j,k)
               wmm(i,j,k) = wm(i,j,k)
#if defined(HTOP)
c               wm(i,j,k) = way(i,j,k)*wm(i,j,k)
#endif
          enddo
        enddo
      enddo
#if defined(CHANNEL)
        if(jblkstart(my_pe+1).eq.1) then
          do i=0,idim1
            do k=0,kdim1
              v(i,0,k) = 0.0
              vm(i,0,k) = 0.0
              vmm(i,0,k) = 0.0
              v(i,1,k) = 0.0
              vm(i,1,k) = 0.0
              vmm(i,1,k) = 0.0
            enddo
          enddo
        endif
        if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
          do i=0,idim1
            do k=0,kdim1
c              v(i,jdim,k) = 0.0
              v(i,jdim1,k) = 0.0
              vm(i,jdim1,k) = 0.0
              vmm(i,jdim1,k) = 0.0
            enddo
          enddo
        endif
#endif

c
c
c      do j=0,jdim1
c         do i=0,idim1
c            vm(i,j,kdim1) = -vm(i,j,kdim)
c            vmm(i,j,kdim1) = -vmm(i,j,kdim)
c            umm(i,j,kdim1) = -umm(i,j,kdim)
c            um(i,j,kdim1) = -um(i,j,kdim)
c         enddo
c      enddo

#if defined(RADIATE) || defined(RADTOP)
      do j=0,jdim1
         do i=0,idim1
            vm(i,j,kdim1) = vm(i,j,kdim)
            vmm(i,j,kdim1) = vmm(i,j,kdim)
            umm(i,j,kdim1) = umm(i,j,kdim)
            um(i,j,kdim1) = um(i,j,kdim)
            umm(i,j,0) = umm(i,j,1)
            um(i,j,0) = um(i,j,1)
            vmm(i,j,0) = vmm(i,j,1)
            vm(i,j,0) = vm(i,j,1)
c
c constant gradient for stress condition
c
c            um(i,j,kdim1) = 2.*um(i,j,kdim)-um(i,j,kdim-1)
c
         enddo
      enddo
#endif
c
c set up slope weight
c
      sangle = atan(mslope)
      angwgt = 1./sin(sangle)
c
c vertical
c      mslope = 1.0
c      angwgt = 1.0
c
c discharge set
c
#if defined(DISCHARGE)
      if(n*delt.lt.120.) then
         disch = dischinit*n*delt/120.
      else
         disch = dischinit
      endif
#else
      disch = 0.
#endif
c      disch = 150.
c compute spd
c
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
            spd1(i,j,k) = sqrt(umm(i+1,j,k)**2+vmm(i+1,j,k)**2+
     *          wmm(i+1,j,k)**2)
          enddo
        enddo
      enddo
      call imag(spd1)
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
c            spd(i,j,k) = (spd1(i+1,j,k)+spd1(i,j,k))*0.5
            spd(i,j,k) = spd1(i+1,j,k)
          enddo
        enddo
      enddo
      call imag(spd)

c
      do i=1,idim
c#if defined(UNDERCUT)
c       if(iblkstart(my_pe+1)+i.lt.30) then
c#else
c       if(iblkstart(my_pe+1).le.1) then
c#endif
        do j=0,jdim
          jj = jblkstart(my_pe+1)+j
c          if(jj.le.745.or.jj.ge.846) then
          do k=1,kdim
c          dxcube = dxx*dxx*dxx
c          if(n*delt.gt.40.) disch = 0.0
c          do k = 1,kdim
c           if(n*delt.lt.6) then
c             dis = rand(i,j)*disch
c           else
c             dis = disch
c           endif
c           dis = disch*delt/400.
           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c           if(bway(i,j,k).eq.1) then
#if defined(DISCHARGE)
           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then

             if(k.lt.dphip+1) then
               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               ax(i,j,k) = 1.0
               ay(i,j,k) = 1.0
c               vay(i,j,k) = 1.0
c               way(i,j,k) = 1.0
c               way(i,j,k+1) = 1.0
c               if(vay(i,j,k+1).eq.0.0) way(i,j,k+1)=0.0

               th(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
               tke(i-1,j,k) = sedinit
             endif
           endif
#endif
c
c input discharge plume
c
c           dcor(i,j,k) = wm(i+1,j,k)
c           dpres(i,j,k) = vm(i+1,j,k)
c           sgf(i,j,k) = th(i,j,k)
           sgf(i,j,k) = tke(i,j,k)
c           tke(i,j,k) = rhop(i+1,j,k)+rhob(k)
c           tke(i,j,k) = w(i+1,j,k)
#if defined(MELT)              
c
c use 4 open boxes adjacent to ice
c
c            if(i.lt.idim) then
            if(mslope.ge.0) then
             if(ax(i,j,k+1).ne.0.0) then
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k+1)+v(ii+1,j,k+1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k+1)+u(ii+1,j,k+1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k+1)+w(ii+1,j,k+1)
               vedge = spd(i,j,k)+spd(i,j,k+1)+
     *                 spd(i+1,j,k)+spd(i+1,j,k+1)
             else
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k-1)+v(ii+1,j,k-1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k-1)+u(ii+1,j,k-1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k-1)+w(ii+1,j,k-1)
               vedge = spd(i,j,k)+spd(i,j,k-1)+
     *                 spd(i+1,j,k)+spd(i+1,j,k-1)

             endif
            else
             if(ax(ii,j,k-1).ne.0.0) then
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k-1)+v(ii+1,j,k-1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k-1)+u(ii+1,j,k-1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k-1)+w(ii+1,j,k-1)
               vedge = spd(i,j,k)+spd(i,j,k-1)+
     *                 spd(i+1,j,k)+spd(i+1,j,k-1)

             else
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k+1)+v(ii+1,j,k+1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k+1)+u(ii+1,j,k+1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k+1)+w(ii+1,j,k+1)
               vedge = spd(i,j,k)+spd(i,j,k+1)+
     *                 spd(i+1,j,k)+spd(i+1,j,k+1)

             endif
            endif
            vbc = 0.25*vedge
c            sgf(i,j,k) = vbc
c            sgf(i,j,k) = th(i+1,j,k)

c            else 
c            ii = i
c            if(mslope.ge.0) then
c             if(ax(i,j,k).ne.0.0) then
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k-1)+v(ii+1,j,k-1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k-1)+u(ii+1,j,k-1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k-1)+w(ii+1,j,k-1)
c             else
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k+1)+v(ii+1,j,k+1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k+1)+u(ii+1,j,k+1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k+1)+w(ii+1,j,k+1)
c             endif
c            else
c             if(ax(ii,j,k-1).ne.0.0) then
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k-1)+v(ii+1,j,k-1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k-1)+u(ii+1,j,k-1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k-1)+w(ii+1,j,k-1)
c             else
c              vedge = v(ii,j,k)+v(ii+1,j,k)+v(ii+1,j,k+1)+v(ii+1,j,k+1)
c              uedge = u(ii,j,k)+u(ii+1,j,k)+u(ii+1,j,k+1)+u(ii+1,j,k+1)
c              wedge = w(ii,j,k)+w(ii+1,j,k)+w(ii+1,j,k+1)+w(ii+1,j,k+1)
c             endif
c            endif
c            vbc = 0.5*vedge
c            ubc = 0.5*uedge
c            wbc = 0.5*wedge

c            endif
c
             
c
#if defined(MCPHEE)
             ustarice = sqrt(cd*v(i,j,k)*v(i,j,k)
     *          +cd*w(i,j,k)*w(i,j,k)
     *          +cd*0.005*0.005)
c             fice = 0.0
             coef = 1.57*sqrt(ustarice*znot/kmb)
             phiturb = 2.5*log(0.5/znot)
c             phiheat = phiturb+2.5*log(znot*0.4*ustarice/khb)
c             phisal = phiturb+2.5*log(znot*0.4*ustarice/kms)

             phiheat = phiturb + coef*(kmb/khb)**0.6667
             phisal = phiturb + coef*(kmb/kms)**0.6667
             c2 = lh*phiheat/phisal
c
c Jenkins
c
c             gammat = 0.011
c             gammas = 3.1e-4
c
c             c2 = gammas*lh/gammat

c
c             b = t+phiheat*fice/ustar+c2
             bb = th(i,j,k)+c2
             aa = m
             cc = -c2*sal(i,j,k)

             salwall = (-bb+sqrt(bb*bb-4.*aa*cc))/(2.*aa)


             wice = ustarice*gams*(sal(i,j,k)-salwall)/salwall

c
c bulk aerodynamic flux of heat
c
c            wice =(th(i,j,k)+0.0054*sal(i,j,k))*cd*
c     *            sqrt(w(i,j,k)**2+v(i,j,k)**2+0.002)/lh
c            wice = wice-wice*ran2(icrseed)/50.
c
c add 2 m/d background melt
c
c            wice = wice+2.0/86400.
c            wice = wice-wice*ran2(icrseed)/50.
#endif
#if defined(THREEEQ)
             zdep = (kdim-k)*dz+0.5*dz
c             aa = lamb1*(gamt*cw-gams*ci)
c             bb = gams*ci*(lamb1*sal(i,j,k)-lamb2-lamb3*zdep-(latf/ci))-
c     *            gamt*cw*(th(i,j,k)-lamb2-lamb3*zdep)
c             cc = gams*sal(i,j,k)*(ci*(lamb2+lamb3*zdep)+latf)
c
c revised
             aa = lamb1
             bb = lamb2+lamb3*zdep-gams*lh/gamt-th(i,j,k)
             cc = sal(i,j,k)*gams*lh/gamt
c
             salwall = (-bb-sqrt(bb*bb-4.*aa*cc))/(2.*aa)
c
c             ustarice = sqrt(cd*v(i,j,k)*v(i,j,k)
c     *          +cd*w(i,j,k)*w(i,j,k)
c     *          +cd*0.001*0.001)
             ustarice = sqrt(cd*vbc*vbc
     *          +cd*0.001*0.001)

             wice = rhob(1)*ustarice*gams*(sal(i,j,k)-salwall)/
     *              (salwall*920.)
             if(iope.and.j.eq.1.and.k.eq.50) then
                write(*,*)"wice,salwall,ustarice ",wice,salwall,ustarice
                write(*,*)"i,s,th ",i,w(i,j,k),sal(i,j,k),th(i,j,k)
             endif

#endif
#if defined(SCHULZ)
             zdep = (kdim-k)*dz+0.5*dz
c             aa = lamb1*(gamt*cw-gams*ci)
c             bb = gams*ci*(lamb1*sal(i,j,k)-lamb2-lamb3*zdep-(latf/ci))-
c     *            gamt*cw*(th(i,j,k)-lamb2-lamb3*zdep)
c             cc = gams*sal(i,j,k)*(ci*(lamb2+lamb3*zdep)+latf)
c
c revised
             aa = lamb1
             bb = lamb2+lamb3*zdep-0.07*lh-th(i,j,k)
             cc = sal(i,j,k)*0.07*lh
c
             salwall = (-bb-sqrt(bb*bb-4.*aa*cc))/(2.*aa)
c
c             ustarice = sqrt(cd*v(i,j,k)*v(i,j,k)
c     *          +cd*w(i,j,k)*w(i,j,k)
c     *          +cd*0.001*0.001)
             ustarice = sqrt(cd*vbc*vbc
     *          +cd*0.00001*0.00001)
             if(vbc.lt.0.05) then
               tranv = 0.0002
             else
               vbc = vbc-0.05
               ustarice = sqrt(cd*vbc*vbc
     *          +cd*0.00001*0.00001)

               tranv = 0.0002+2.2e-2*ustarice
             endif

             wice = rhob(1)*0.07*tranv*(sal(i,j,k)-salwall)/
     *              (salwall*920.)
             if(iope.and.j.eq.1.and.k.eq.50) then
                write(*,*)"wice,salwall,ustarice ",wice,salwall,ustarice
                write(*,*)"i,s,th ",i,w(i,j,k),sal(i,j,k),th(i,j,k)
             endif

#endif

c
c             if(n*delt.lt.900.) then
c               wice = wice-wice*ran2(icrseed)/50.
c             endif
#if defined(DISCHARGE)
c
c set wice to zero in discharge area
c
c           if(jj.gt.jtotal/2-30.and.jj.lt.jtotal/2+30) then
c           if(jj.gt.jtotal/2+200.and.jj.lt.jtotal/2+300) then
c           if(jj.gt.jtotal/2-50.and.jj.lt.jtotal/2+50) then
           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
               wice = 0.0
             endif
           endif
#endif

#endif
#if defined(NOMELT)
           wice = 0.0
#endif

c
c prescribe melt at 2e-6
c
c             wice = 2.0e-6

c
c 0.9 cm/hr
c            if(k.gt.5) then
c              if(n*delt.lt.1.) then
c                wice = n*delt*0.009/(3600.)
c              else
c                wice = 0.009/3600.
c              endif
cc              wice = wice+wice*ran2(icrseed)/20.
c            else
c              wice = 0.0
c            endif
c             wicen = 1.12*khb*(th(1,j,k)-th(0,j,k))*odx/lh
c             if(iope) then
c             if(k.eq.100.and.j.eq.1)write(*,*)"wice,s,t w",
c     *         wice,wicen,sal(i,j,k),sal(i-1,j,k),ax(i,j,k),ax(i+1,j,k)
c             endif
c            wice = wicen
c
c
c set wice to  5 m day
c
c endif for ax option


c           if(ax(i,j,k).eq.0.0.and.ax(i,j,k).ne.0.0) then
c            wice = 6./86400.
c
c add 4 m d-1 background melt
c
c            wice = wice+4./86400.
c            um(i,j,k) = wice
c
            heatflux = -1.*wice*lh/dxx
            saltflux = -1.*wice*sal(i,j,k)/dxx
c
c            dpres(i,j,k) = wice
c            sgf(i,j,k) = wice
            dpres(i,j,k) = vm(i+1,j,k)
           dcor(i,j,k) = th(i+1,j,k)
c            dcor(i,j,k) = um(i+1,j,k)
            dvadvv(i,j,k) = sal(i+1,j,k)
c            dvadvv(i,j,k) = vm(i,j,k)
c            
c
            if(th(i,j,k).lt.-1.5) then
                  heatflux = 0.0
                  saltflux = 0.0
            endif
            if(sal(i,j,k).lt.0.5) then
                  saltflux = 0.0
                  heatflux = 0.0
            endif
c            tke(i,j,k) = wice
c            sal(i,j,k) = sal(i,j,k)+saltflux*delt
            if(sal(i,j,k).lt.0.0) sal(i,j,k) = 0.0
c            th(i,j,k) = th(i,j,k)+heatflux*delt
            if(th(i,j,k).lt.-0.054*sal(i,j,k)) then
                th(i,j,k) = -0.054*sal(i,j,k)
            endif
           endif
         enddo
        enddo
c        endif
        enddo

        call imag(th)
        call imag(sal)
c
c diffuse heat and salt from ice edge
c
#if defined(DIFFMELT)
        do i=1,idim
         do j=1,jdim
          do k=1,kdim
                 th(i,j,k) = th(i,j,k) + delt*
     *           khb*((th(i+1,j,k)+th(i-1,j,k)-2.*th(i,j,k))/
     *             (dx(i)*dx(i))+(th(i,j-1,k)+
     *             th(i,j+1,k)+th(i,j,k+1)+th(i,j,k-1)-4.*th(i,j,k))
     *             *tdxs)
c
                sal(i,j,k) = sal(i,j,k) + delt*
     *           kms*((sal(i+1,j,k)+sal(i-1,j,k)-2.*sal(i,j,k))/
     *             (dx(i)*dx(i))+(sal(i,j-1,k)+
     *             sal(i,j+1,k)+sal(i,j,k+1)+sal(i,j,k-1)-4.*sal(i,j,k))
     *             *tdxs)

          enddo
         enddo
        enddo
c        do i=1,idim
c         do j=1,jdim
c          do k=1,kdim
c             th(i,j,k) = thp(i,j,k)
c             sal(i,j,k) = salp(i,j,k)
c          enddo
c         enddo
c        enddo
#endif
        if(iblkstart(my_pe+1).eq.1) then
             i=1
             j=1
             k=3
             if(iope) then
             if(k.eq.3.and.j.eq.1)write(*,*)"out1,s,t w",
     *         wice,sal(i,j,k),sal(i-1,j,k),
     *         th(i,j,k),th(i-1,j,k),wm(i,j,k)
             endif
        endif


        call imag(th)
        call imag(sal)
        call obcmom(rv,n)


c 
      call deltime(cputime)
      veltim=veltim+cputime
c 
c               umm   u    um
c      +    .    +    .    +
c     t-1        t  t+1/2 t+1
c 
c advance theta, salt, and tke from t to t+1
c 
#if defined(DISCHARGE)
        do k=1,kdim
         do j=1,jdim
            do i=1,idim

           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
           jj = jblkstart(my_pe+1)+j

           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               ax(i,j,k) = 1.0
               th(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
               tke(i-1,j,k) = sedinit
             endif
           endif
          endif
         enddo
         enddo
      enddo
      call imag(th)
      call imag(sal)
      call imag(tke)
      call imag(u)
      call imag(um)
      call imag(umm)
#endif

c
#if defined(USEMPI)
      call MPI_BARRIER(comm,mpierr)
#endif
#if defined(LEVEQ)
      call leveq(n,th,thp,MDTEMP)
#endif
c
#if defined(LEER)
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=0,idim1
c            kheat(i,j,k) = kheat(i,j,k)+khb
c          enddo
c        enddo
c      enddo        
#if defined(HTOP) || defined(ICE)
      call leer(n,th,thp,MDTEMP)
#else
      call leer(n,th,thp,MDTEMP)
#endif
#endif
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=0,idim1
c            kheat(i,j,k) = kheat(i,j,k)-khb
c          enddo
c        enddo
c      enddo

       do k=1,kdim
         do j=1,jdim
           tb = 0.0
           do i=1,idim
             tb = tb+th(i,j,k)
           enddo
           t_bar(j,k) = tb/dble(idim)
         enddo
       enddo
c       do k=2,kdim
c         do j=1,jdim
c           do i=1,idim
c             wt_ave(k) = wt_ave(k)+wm(i,j,k)*(th(i,j,k)-t_bar(j,k)
c     *                   +th(i,j,k-1)-t_bar(j,k-1))*0.5
c           enddo
c         enddo
c       enddo
c      do k=0,kdim1
c         do j=0,jdim1
c            do i=0,idim1
c               th(i,j,k)=thp(i,j,k)
c            enddo
c         enddo
c      enddo
c 
c      call imag(th)
c 
c test
c
c         if(my_pe.eq.0.and.n.eq.1) then
c           do i=0,idim1
c             do j=0,jdim1
c           do k=0,kdim1
c             sal(i,j,k) = 30.0
c               enddo
c             enddo
c           enddo
c         endif
c         call imag(sal)
c
c
c
      if (advsal) then
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=0,idim1
c            kheat(i,j,k) = kheat(i,j,k)
c          enddo
c        enddo
c      enddo

#if defined(LEER) && defined(SNGLSALT)
#if defined(HTOP) || defined(ICE)
         call leer(n,sal,salp,MDSAL)
#else
         call leer(n,sal,salp,MDSAL)
#endif
#endif
#if defined(LEVEQ) && defined(SNGLSALT)
         call leveq(n,sal,salp,MDSAL)
#endif
#if defined(LEER) && !defined(SNGLSALT)
         call leerdp(n,sal,salp,MDSAL)
#endif
#if defined(LEVEQ) && !defined(SNGLSALT)
         call leveqdp(n,sal,salp,MDSAL)
#endif
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
             if(salp(i,j,k).lt.1.0) salp(i,j,k) = 1.0
c            kheat(i,j,k) = kheat(i,j,k)
             if(salp(i,j,k).lt.0.8) then
                 write(*,*)"sal too small"
     *       ,i+iblkstart(my_pe+1)-1,j+jblkstart(my_pe+1)-1,k,my_pe
              write(*,*) "salp ",salp(i,j,k),salp(i-1,j,k),salp(i+1,j,k)
               write(*,*) "sal ",sal(i,j,k),sal(i-1,j,k),sal(i+1,j,k)
               write(*,*) "u ",u(i,j,k),u(i-1,j,k),u(i+1,j,k)
               write(*,*) "v ",v(i,j,k),v(i,j-1,k),v(i,j+1,k)
               write(*,*) "w ",wm(i,j,k),wm(i,j,k+1)
               write(*,*) "w ",w(i,j,k),w(i,j,k+1)
               write(*,*) "ax ",ax(i,j,k),ax(i-1,j,k),ax(i+1,j,k)
               write(*,*) "ay ",ay(i,j,k),ay(i,j-1,k),ay(i,j+1,k)
               write(*,*) "vay ",vay(i,j,k),vay(i,j-1,k),vay(i,j+1,k)
               write(*,*) "way ",way(i,j,k),way(i,j,k-1),way(i,j,k+1)
             endif
          enddo
        enddo
      enddo

#if defined(OPENBC)
c      call obcscal(sal,salp)
#endif
         do k=0,kdim1
            do j=0,jdim1

#if defined(OPENSPONGE)
c
           ii = iblkstart(my_pe+1)+i
           fmax = itotal - itotal/5
           irelax = ii - fmax
           if(ii.gt.fmax) then
             relax = delt*(sin(1.5708*irelax/fmax)**2)/10.
c             relax = delt/(20.*(irelax/2-ii+1))
           else
             relax = 0.0
           endif

#else
          relax = 0.0
#endif
               do i=0,idim1
c                   sal(i,j,k)=salp(i,j,k)-(us(k)+u_bar(k)+uave)*dsadx
#if defined(OPENSPONGE)
#if defined(ICE)
                  if(k.le.kdist(i,j)) then
                    sal(i,j,k) = salp(i,j,k)
     *               + relax*(salbar(k)-sal(i,j,k))

                  endif
#else
                  sal(i,j,k) = salp(i,j,k)
     *               + relax*(salbar(k)-sal(i,j,k))
#endif

#if defined(HBOT)
                  if(k.gt.kdist(i,j)) then
                    sal(i,j,k) = salp(i,j,k)
     *               + relax*(salbar(k)-sal(i,j,k))
                  endif
#else
                  sal(i,j,k) = salp(i,j,k)
     *               + relax*(salbar(k)-sal(i,j,k))
#endif

#else

#if defined(ICE)
                  if(k.le.kdist(i,j)) then
                    sal(i,j,k) = salp(i,j,k)
                  endif
#else
                  sal(i,j,k) = salp(i,j,k)
#endif

#if defined(HBOT)
                  if(k.gt.kdist(i,j)) then
                    sal(i,j,k) = salp(i,j,k)
c                    if(k.eq.kdim) then
c                      sal(i,j,kdim) = sal(i,j,kdim)-0.01*delt*odz*
c     *                        4.*6.*sal(i,j,kdim)/86400.
c                    endif
                  endif
#else
                  sal(i,j,k) = salp(i,j,k)
#endif
c
c sponge endif
c
#endif
                  if(iope) then
                    if(i.eq.5.and.j.eq.2.and.k.eq.4) then
                      write(*,*)"ubar,rhop ",u_bar(k),rhop(i,j,80)
                    endif
                  endif
               enddo
            enddo
         enddo
#if defined(DISCHARGE)
c           do k=1,kdim
c           do j=1,jdim
c           do i=1,idim
c
c           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c           jj = jblkstart(my_pe+1)+j
c
c           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
c             if(k.lt.dphip+1) then
cc               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
cc               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
cc               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               fracmelt  = delt*disch/(dphip*2*plumw*dxx*dxx)
c               th(i,j,k) = th(i,j,k)*(dxx-fracmelt)+0.5*fracmelt
c               sal(i,j,k) = sal(i,j,k)*(dxx-fracmelt)+1.0*fracmelt
c             endif
c           endif
c          endif
c          enddo
c          enddo
c          enddo
#endif

c 
#if defined(SNGLSALT)
c
         call imag(sal)

#else
         call imagdp(sal)
#endif

      end if
c
#if defined(OPENBC)
c      call obcscal(th,thp)
#endif
      do k=0,kdim1
        do j=0,jdim1
#if defined(OPENSPONGE)
           ii = iblkstart(my_pe+1)+i
           fmax = itotal - itotal/5
           irelax = ii - fmax
           if(ii.gt.fmax) then
             relax = delt*(sin(1.5708*irelax/fmax)**2)/10.
c             relax = delt/(20.*(irelax/2-ii+1))
           else
             relax = 0.0
           endif
#else
          relax = 0.0
#endif

c            th(i,j,k)=thp(i,j,k)-(us(k)+u_bar(k)+uave)*dtadx
c
          do i=1,idim
#if defined(OPENSPONGE)
#if defined(ICE)
            if(k.le.kdist(i,j))then
              th(i,j,k) = thp(i,j,k)
     *          + relax*(thbar(k)-th(i,j,k))
            endif
#else
            th(i,j,k) = thp(i,j,k)
     *          + relax*(thbar(k)-th(i,j,k))
#endif

#if defined(HBOT)
                  if(k.gt.kdist(i,j)) then
                    th(i,j,k) = thp(i,j,k)
     *              + relax*(thbar(k)-th(i,j,k))
                  endif
#else
            th(i,j,k) = thp(i,j,k)
     *          + relax*(thbar(k)-th(i,j,k))
#endif
#else
#if defined(FRONT2)
c
c add frontal temperature advection
c
c            if(k.gt.kdim-6) then
c              thp(i,j,k) = thp(i,j,k)+0.001*(u(i,j,k)+u(i+1,j,k))
c     *                   *0.5*delt
c            else
c              thp(i,j,k) = thp(i,j,k)+0.0004*(u(i,j,k)+u(i+1,j,k))
c     *                   *0.5*delt
c              thp(i,j,k) = thp(i,j,k)+0.0001*(u(i,j,k)+u(i+1,j,k))
c     *                   *0.5*delt
c              thp(i,j,k) = thp(i,j,k)+0.00002*(u(i,j,k)+u(i+1,j,k))
c              thp(i,j,k) = thp(i,j,k)+0.000216*(u(i,j,k)+u(i+1,j,k))
c     *                   *0.5*delt
c              thp(i,j,k) = thp(i,j,k)+0.00005*(u(i,j,k)+u(i+1,j,k))
c     *                   *0.5*delt
c              thp(i,j,k) = thp(i,j,k)+0.00002*(u(i,j,k)+u(i+1,j,k))
c            endif
#endif

#if defined(ICE)
            if(k.le.kdist(i,j))then
              th(i,j,k) = thp(i,j,k)
            endif
#else
            th(i,j,k) = thp(i,j,k)
#endif

#if defined(HBOT)
                  if(k.ge.kdist(i,j)) then
                    th(i,j,k) = thp(i,j,k)
                  endif
#else
            th(i,j,k) = thp(i,j,k)
#endif

#endif

	  enddo
        enddo
      enddo

      call imag(th)
c
c reset background density
c
c         rho0=0.0
c         drho=0.0
c         do k=1,kdim
c            drhol=0.0
c            do j=1,jdim
c               do i=1,idim
c                  s=sal(i,j,k)
cc                  drhol=drhol+iestbl(s,th(i,j,k),5)
c                  drhol = drhol+sigthta(s,th(i,j,k),5)
c               enddo
c            enddo
c            call MPI_ALLREDUCE(drhol,drho,1,
c     &        MPI_REAL8,MPI_SUM,comm, ierr)
c
c            rhobm1(k)=drho/(jtotal*itotal)
c            rhob(k)=1000.0+rhobm1(k)
cc            if(iope) then
cc              write(*,*)"k,rhob,t,s ",k,rhob(k),th(1,1,k),sal(1,1,k)
cc            endif
c            rho0=rho0+rhob(k)
c         enddo
c         rho0=rho0/float(kdim)
c         rhob(0) = rhob(1)
c         rhob(kdim1) = rhob(kdim)
c

c        if(iblkstart(mype+1).eq.1) then
c             i=1
c             j=1
c             k=150
c             if(iope) then
c             if(k.eq.150.and.j.eq.1)write(*,*)"out2,s,t w",
c     *         wice,sal(i,j,k),sal(i-1,j,k),
c     *         th(i,j,k),th(i-1,j,k),wm(i,j,k),rhobm1(1)
c             endif
c        endif

#if defined(FRONT)
c      do k=0,kdim1
c        do j=0,jdim1
c          th(0,j,k) = th(1,j,k)
c          th(idim1,j,k) = th(idim,j,k)
c        enddo
c      enddo
#endif
   
#if defined(MPOND)
c      call mpond(n)
#endif
#if defined(ICE)
      if(my_pe.eq.10) then
        write(*,*)"th(10,10,kdim) before",th(10,10,kdim)
      endif
      call icemod(n,nerr)
c      if(mod(int(n*delt),60).eq.0) then
c      do i=1,idim
c        do j=1,jdim
c          if(hice(i,j).eq.0.0) then
c          do kk=1,kdim
c                   th(i,j,kk)=th(i,j,kk)+
c     *                  60.*swheatf*rad(kk)/(cp*
c     *                  rho0*dz)
c
c          enddo
c          th(i,j,kdim)=th(i,j,kdim)+
c     *         60.*heatf*randn(i,j)/(cp*rho0*dz)
c          th(i,j,kdim1) = th(i,j,kdim)
c          endif
c        enddo
c      enddo
c      endif
      call imag(sal)
      call imag(th)
#if defined(ICEWALL)
      call wallice(n)
#endif
#if !defined(WARMICE)
      call imag(cfr)
      call leer(n,cfr,cfrp,MDCFR)
      do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
            cfr(i,j,k) = cfrp(i,j,k)
          enddo
        enddo
      enddo
#endif
#endif
      if (usesgs.and.advtke) then
        do k=0,kdim1
        do j=0,jdim1
          do i=0,idim1
c
c add settling velocity for sediment
c
            w(i,j,k) = w(i,j,k)-0.02
          enddo
        enddo
        enddo

#if defined(LEVEQ)
         call leveq(n,tke,tkep,MDTKE)
#endif
#if defined(LEER)
         call leer(n,tke,tkep,MDTKE)
#endif
c
c dont need to update tke, it will be performed by deardorff
c 
         call imag(tkep)
         do k=0,kdim1
         do j=0,jdim1
          do i=0,idim1
c
c remove settling velocity for sediment
c
            w(i,j,k) = w(i,j,k)+0.02
          enddo
        enddo
      enddo

      end if
#if defined(TRACER)
c 
c advance tracers from t to t+1
c 
      do itrc=1,ntrace
#if defined(LEVEQ)
         call leveq(n,trc(0,0,0,itrc),trcp,MDTRCR+itrc)
#endif
#if defined(LEER)
         call leer(n,trc(0,0,0,itrc),trcp,MDTRCR+itrc)
#endif
         do k=0,kdim1
            do j=1,jdim
               do i=0,idim1
                  trc(i,j,k,itrc)=trcp(i,j,k)
               enddo
            enddo
         enddo
c     
         call imag(trc(0,0,0,itrc))
        
      enddo
#endif
c 
      call deltime(cputime)
      advtim=advtim+cputime
c
c      if(havefloats)then
c         call flt_update(n)
c         call deltime(cputime)
c         floattim=floattim+cputime
c      endif
c 
c compute perturbation density
c 
      call denupdte
c
c adjust for sediment load
c
c      do j=1,jdim
c        do k=1,kdim
c          do i=1,idim
c            rhot = rhobm1(k) + rhop(i,j,k)
c            rhot = rhot+(1.0-rhot/2650.)*tke(i,j,k)
c            rhop(i,j,k) = rhot-rhobm1(k)
c          enddo
c       enddo
c      enddo
c         testt = iestbl(sal(4,4,4),th(4,4,4),5)
c         siggt = sigthta(sal(4,4,4),th(4,4,4),5)
c         write(*,*)"my_pe iestbl ",my_pe,testt,sal(4,4,4),
c     *     th(4,4,4),siggt
c
c         do j=1,jdim
c            do k=1,kdim
c               do i=1,idim
c                  kki = 5
c                  s=sal(i,j,k)
cc                  rhop(i,j,k)=iestbl(s,th(i,j,k),2*k)-rhobm1(k)
c                  rhop(i,j,k)=sigthta(s,th(i,j,k),kki)-rhobm1(k)
c               enddo
c            enddo
c            do i=1,idim
c               rhop(i,j,0)=rhop(i,j,1)
c               rhop(i,j,kdim1)=rhop(i,j,kdim)
c            enddo
c         enddo
c         call imag(rhop)
c
c make theta a tracer
c
c      do k=0,kdim1
c        do j=0,jdim1
c          do i=0,idim1
c            rhop(i,j,k) = 0.0
c          enddo
c        enddo
c      enddo
c
c compute new potential energy and effect of subgrid mixing.
c
#if defined(ANALY)
      do k=1,kdim
        pel(k) = 0.0
        zz = (k-1)*dz + 0.5*dz
        do j=1,jdim
          do i=1,idim
            pel(k) = pel(k) + (rhop(i,j,k)+rhobm1(k))*zz*9.81/rhob(k)
          enddo
        enddo
      enddo
c
c do processor summation
c
#if defined(USEMPI)
      call MPI_REDUCE(pel,pen,kdim,
     &        MPI_REAL8,MPI_SUM,mastertask,
     &        comm, ierr)
      call MPI_BARRIER(comm,nerr)
#else
      do k=1,kdim
        pen(k) = pel(k)
      enddo
#endif
c
c compute average rate of change due to subgrid term
c
c      if(iope) then
c        do k=1,kdim
c          dpesg(k) = dpesg(k)+(pen(k)-peo(k))/delt
c        enddo
c      endif
#endif
c 
c end of time step b.c. updates
c 
      call bc_update(n)
c 
c call subgrid scale model update
c
c if advecting tke, then new field is in tkep and calling order
c is reversed for deardorff so that tkep goes into tke
c 
      if (usesgs) then
c 
         if (advtke.and.usedeardorff) then
            call deardorff(n,tkep,tke)
            do k=0,kdim1
               do j=1,jdim
                  do i=0,idim1
                     kmom(i,j,k)=kmomp(i,j,k)
                  enddo
               enddo
            enddo
         else if(usedeardorff) then
            call deardorff(n,tke,tkep)
            do j=1,jdim
               do k=0,kdim1
                  do i=0,idim1
                     tke(i,j,k)=tkep(i,j,k)
                     kmom(i,j,k)=kmomp(i,j,k)
                  enddo
               enddo
            enddo
         else if(.not.usedeardorff)then
c
c if not using deardorff then use Metais
c
            call metais(n)
c
c if not using tke then use smagorinsky
c
c     call smag(n)           
c
         end if
c 
         if(advtke) then
         do j=0,jdim1
           do i=0,idim1
             do k=0,kdim1
cc
cc include 0.01 m/s settling velocity
cc
c               tke(i,j,k) = tkep(i,j,k)-0.01*tkep(i,j,k)*delt
               tke(i,j,k) = tkep(i,j,k)
             enddo
           enddo
         enddo
         call imag(tke)
         endif
         call imag(kmom)
         call imag(kheat)
c     
         kmm = 0.1*dz*dz/delt
         do k=1,kdim
            do j=1,jdim
               do i=1,idim
                  if(kmom(i,j,k).gt.kmm.or.kmom(i,j,k).lt.0.0)then
                     write (lunit,*) 'Momentum diff. incorrect at '
     *                    ,i,' ',k,' with value  ',kmom(i,j,k),kmm
                     nerr=1
                     return 
                  end if
               enddo
            enddo
         enddo
      else
         do k=0,kdim1
            do j=0,jdim1
               do i=0,idim1
                  kmom(i,j,k)=kmb
                  kheat(i,j,k)=khb
               enddo
            enddo
         enddo
      end if
c 
c
      if(usetopref) then
        do k=0,kdim1
          do j=0,jdim1
            do i=0,idim1
              umm(i,j,k) = umm(i,j,k)+uave
              um(i,j,k) = um(i,j,k)+uave
              vmm(i,j,k) = vmm(i,j,k)+vave
              vm(i,j,k) = vm(i,j,k)+vave
            enddo
          enddo
        enddo
      endif
#if defined(CHANNEL)
        if(jblkstart(my_pe+1).eq.1) then
          do i=0,idim1
            do k=0,kdim1
              v(i,0,k) = 0.0
              v(i,1,k) = 0.0
            enddo
          enddo
        endif
        if(jblkstart(my_pe+1)+jdim.ge.jtotal) then
          do i=0,idim1
            do k=0,kdim1
c              v(i,jdim,k) = 0.0
              v(i,jdim1,k) = 0.0
              vm(i,jdim1,k) = 0.0
              vmm(i,jdim1,k) = 0.0
            enddo
          enddo
        endif
#endif
      call deltime (cputime)
      updtetim=updtetim+cputime
c     
      return
      end
