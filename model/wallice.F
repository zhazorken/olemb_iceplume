      subroutine wallice(n)
c
c This routine calculates which grid volumes are adjacent to ice in the
c horizontal direction, for example at a lead edge.  The routine figures
c out the lateral heat and salinity flux between the ice and the water.
c Currently, it is set for melting ice and assumes the heat flux from 
c the ice is negligible.
c
#if defined(USEMPI)
#include "mpif.h"
#endif

#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "ice.inc"
#include "finitv.inc"
c
      integer i,j,k,n,ii,jj
      real cd,odepth,vstar,spd,drag
c
c first loop over all points and find those volumes that have an ice
c edge
c
      cd = 0.16/(log(odepth/0.005))**2
      do i=1,idim
        odepth = 0.5*dx(i)
        odx = 1./dx(i)
        do j=1,jdim
          do k=1,kdist(i,j)
c
c if k is greater than neighbor kdist, then there is
c ice between this volume and the adjacent point.
c Compute ustar based on the vertical velocity and call
c iceflux to figure out heat/salt from from melting.
c
            if (kdist(i-1,j).lt.k) then
                ii = i-1
                jj = j
                spd = wm(i,j,k)**2+vm(i,j,k)**2
                vstar = sqrt(spd*cd)
                call iceflux(i,j,k,ii,jj,vstar)
c
c include drag term on velocity
c
                drag = -cd*(vm(i,j,k))*abs(vm(i,j,k))
                vm(i,j,k) = vm(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt
c
            endif
            if (kdist(i+1,j).lt.k) then        
                ii = i+1
                jj = j
                spd = wm(i,j,k)**2+vm(i,j,k)**2
                vstar = sqrt(spd*cd)
                call iceflux(i,j,k,ii,jj,vstar)
c
                drag = -cd*(vm(i,j,k))*abs(vm(i,j,k))
                vm(i,j,k) = vm(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt

            endif
            if (kdist(i,j-1).lt.k) then        
                ii = i
                jj = j-1
                spd = wm(i,j,k)**2+um(i,j,k)**2
                vstar = sqrt(spd*cd)
c                if(i.eq.idim/2) then
c                  write(*,*)"wallice vstar,spd,ii,jj",vstar,spd,ii,jj
c                endif
                call iceflux(i,j,k,ii,jj,vstar)
                drag = -cd*(um(i,j,k))*abs(um(i,j,k))
                um(i,j,k) = um(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt

            endif
            if (kdist(i,j+1).lt.k) then        
                ii = i
                jj = j+1
                spd = wm(i,j,k)**2+um(i,j,k)**2
                vstar = sqrt(spd*cd)
                call iceflux(i,j,k,ii,jj,vstar)
                drag = -cd*(um(i,j,k))*abs(um(i,j,k))
                um(i,j,k) = um(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt
            endif
          enddo
        enddo
      enddo
      call imag(sal)
      call imag(th)
      return
      end

      subroutine iceflux(i,j,k,ii,jj,vstar)

#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "ice.inc"
#include "communicate.inc"
c
      real deldepth,thy,saly,fice,znot,vstar,lfusion,lh
      real m,dss,c2,phiheat,phiturb,phisal,a,b,c,sw
      real wice,coef,kond,cd,odepth,swallflux,hwallflux
      integer i,j,k,ii,jj

c
c
        odepth = 0.5*dz
        deldepth = dz
        znot = 0.005
        lfusion = 3.34e5
        lh = lfusion/cp
        kms = 7.4e-10
        m = 0.054
c
        thy = th(i,j,k)
        saly = sal(i,j,k)
        if(vstar.gt.0.000001) then

c
          fice = 0.0
c
          coef = 1.57*sqrt(vstar*znot/kmb)
          phiturb = 2.5*log(odepth/znot)
          phiheat = phiturb + coef*(kmb/khb)**0.66667
          phisal = phiturb + coef*(kmb/kms)**0.66667
          dss = 0.0
          c2 = (lh+dss)*phiheat/phisal
c
c
c solve quadratic equation
c
          b = thy+phiheat*fice/vstar - 
c          b = thy - 
     *            m*sice(ii,jj) + c2
c
          a = m
c
c          c = -thy*sice(ii,jj) - phiheat*fice*
c     *            sice(ii,jj)/vstar - c2*saly
          c = -thy*sice(ii,jj) -  c2*saly
c
          sw = (-b + sqrt(b*b - 4.*a*c))/(2.*a)
c
c compute ice increase
c
          if(abs(sw-sice(ii,jj)).gt.0.5) then
            wice = vstar*(saly - sw)/(phisal*
     *             (sw - sice(ii,jj)))
          else
            wice = vstar*(thy+m*sw)/(lh*phiheat)
          endif
c
        else
          wice = 0.0
        endif
        
c
        swallflux = -1.*wice*(saly-sice(ii,jj))/deldepth
        hwallflux = -1.*wice*lh/deldepth 
c     *                     delt*(fice+wice*dss)/deldepth
#if defined(ICEWALL)
        wallsf(i,j,k) = swallflux*deldepth
        wsf = wsf + wallsf(i,j,k)
        wallhf(i,j,k) = hwallflux*deldepth 
        whf = whf + wallhf(i,j,k)
        sal(i,j,k) = sal(i,j,k)+delt*swallflux
        th(i,j,k) = th(i,j,k)+delt*hwallflux
        wallarea = wallarea+1.
#endif

c        if(i.eq.idim/2.and.k.gt.35) then
c          write(*,*)"hwallflux, th,sice,j ",hwallflux,th(i,j,k),
c     *       sice(ii,jj),j
c          write(*,*)"wice,vstar is ",wice,vstar
c        endif
        return
        end
