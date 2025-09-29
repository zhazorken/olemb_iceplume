      subroutine mpond(n)
c
c mpond simulates fluxes at the edge of a melt pond.
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
#include "cyclic.inc"
#include "ice.inc"
#include "finitv.inc"
#include "communicate.inc"
c
      integer i,j,k,n,ii,jj,kk
      real cd,odepth,vstar,spd,drag
      real swave
c
c first loop over all points and find those volumes that have an ice
c edge
c
      odepth = 0.5*dx
      cd = 0.16/(log(odepth/0.005))**2
      do i=1,idim
        do j=1,jdim
c
c first calculate the heat flux in the pond from surface fluxes and
c solar radiation
c
c if the bottom intersects the surface, then we are out of the pond and
c on to the ice surface
c
          if(hbot(i,j).lt.kdim*dz) then
            if(mod(n,int(3./delt)).eq.0) then
              if(swheatf.gt.0.0) then
                swave = 0.0
                do kk=kdist(i,j),kdim
                   th(i,j,kk)=th(i,j,kk)+
     *                  3.*swheatf*rad(kk)/(cp*
     *                  rho0*dz)
                   swave = swheatf*rad(kk)+swave
                enddo
                if(my_pe.eq.8) then
                  if(i.eq.200.and.j.eq.8) then
                     write(*,*)"solar heating ", n*delt,swave, heatf
                  endif
                endif
              endif
              th(i,j,kdim)=th(i,j,kdim)+
     *             3.*heatf/(cp*rho0*dz)
              th(i,j,kdim1) = th(i,j,kdim)
c
c add in surface melt water of 1 cm per day with lead coverage
c of 5% (19 parts for 1 part lead)
c take half of that to account for melt ponds
c
c              sal(i,j,kdim) = sal(i,j,kdim)-0.01*delt*odz*
c     *                        19.*6.*sal(i,j,kdim)/86400.
c     *                        0.5*19.*6.*sal(i,j,kdim)/86400.
              sal(i,j,kdim1) = sal(i,j,kdim)
c
            endif
c
c call ice flux to get fluxes of s and t on the bottom of the pond
c
            call mpondflux(i,j,kdist(i,j),i,j,uistar(i,j))

          endif
c

          if(kdist(i,j).ne.kdim) then
c
            do k=kdist(i,j),kdim
c
c if k is greater than neighbor kdist, then there is
c ice between this volume and the adjacent point.
c Compute ustar based on the vertical velocity and call
c mpondflux to figure out heat/salt from from melting.
c
              if (kdist(i-1,j).gt.k) then
                ii = i-1
                jj = j
                spd = wm(i,j,k)**2+vm(i,j,k)**2
                vstar = sqrt(spd*cd)
                call mpondflux(i,j,k,ii,jj,vstar)
c
c include drag term on velocity
c
                drag = -cd*(vm(i,j,k))*abs(vm(i,j,k))
                vm(i,j,k) = vm(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt
c
              endif
              if (kdist(i+1,j).gt.k) then        
                ii = i+1
                jj = j
                spd = wm(i,j,k)**2+vm(i,j,k)**2
                vstar = sqrt(spd*cd)
                call mpondflux(i,j,k,ii,jj,vstar)
c
                drag = -cd*(vm(i,j,k))*abs(vm(i,j,k))
                vm(i,j,k) = vm(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt

              endif
              if (kdist(i,j-1).gt.k) then        
                ii = i
                jj = j-1
                spd = wm(i,j,k)**2+um(i,j,k)**2
                vstar = sqrt(spd*cd)
c                if(i.eq.idim/2) then
c                  write(*,*)"wallice vstar,spd,ii,jj",vstar,spd,ii,jj
c                endif
                call mpondflux(i,j,k,ii,jj,vstar)
                drag = -cd*(um(i,j,k))*abs(um(i,j,k))
                um(i,j,k) = um(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt

              endif
              if (kdist(i,j+1).gt.k) then        
                ii = i
                jj = j+1
                spd = wm(i,j,k)**2+um(i,j,k)**2
                vstar = sqrt(spd*cd)
                call mpondflux(i,j,k,ii,jj,vstar)
                drag = -cd*(um(i,j,k))*abs(um(i,j,k))
                um(i,j,k) = um(i,j,k)+drag*odx*delt
                drag = -cd*(wm(i,j,k))*abs(wm(i,j,k))
                wm(i,j,k) = wm(i,j,k)+drag*odx*delt
              endif
            enddo
          endif
c
        enddo
      enddo
      call imag(sal)
      call imag(th)
      return
      end

      subroutine mpondflux(i,j,k,ii,jj,vstar)

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
      integer i,j,k,ii,jj,ix,jy

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
        do ix=1,idim
          do jy=1,jdim
            sice(ix,jy) = 4.0
          enddo
        enddo
        if(vstar.lt.0.000001) vstar = 0.0000011
        if(vstar.ge.0.000001) then

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
        if(i.eq.ii.and.j.eq.jj) then
c
c this is the bottom
c
          iceh(i,j) = hwallflux*deldepth
          icesal(i,j) = swallflux*deldepth
          sicebc = sicebc + swallflux*deldepth
          hicebc = hicebc + hwallflux*deldepth
          wallsf(i,j,k) = swallflux*deldepth
          wallhf(i,j,k) = hwallflux*deldepth 
          bvst = bvst + vstar
          bottoma = bottoma + 1.
        else
c
c otherwise side walls
c
          wallsf(i,j,k) = swallflux*deldepth
          wsf = wsf + wallsf(i,j,k)
          wallhf(i,j,k) = hwallflux*deldepth 
          whf = whf + wallhf(i,j,k)
          wvst = wvst + vstar
          wallarea = wallarea+1.
        endif
        sal(i,j,k) = sal(i,j,k)+delt*swallflux
        th(i,j,k) = th(i,j,k)+delt*hwallflux

        if(i.eq.idim/2.and.k.gt.25) then
          write(*,*)"hwallflux, th,sice,j ",hwallflux,th(i,j,k),
     *       sice(ii,jj),j
          write(*,*)"wice,vstar is ",wice,vstar
        endif
        return
        end
