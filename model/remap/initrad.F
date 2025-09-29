      subroutine initrad
c
c   $Id: initrad.F,v 1.4 1998/11/13 20:52:22 mikeh Exp $
c
c This routine initializes a solar irradiance value as a function of
c depth using fractional irradiance data.  These values are then
c interpolated to a vertical column with grid spacing matching the current
c model run.  The grid column is at the temperature levels.
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "cyclic.inc"
#include "communicate.inc"
c
      real depth,irrad(0:kdim1)
      integer kz,kk
      real tdepth
      real p1,p2,p3,p4
      real k1,k2,k3,k4

c
c values taken from Pegau and Paulson, "The summertime evolution of an
c Arctic lead"
      p1 = 0.481
      p2 = 0.194
      p3 = 0.123
      p4 = 0.202
c
      k1 = 0.18
      k2 = 3.25
      k3 = 27.5
      k4 = 300.

      rad(kdim) = 0.0
      do kz=0,kdim-1
         depth = (kdim-kz)*dz
c
c tropical profile from Bill Smyth
c
c         rad(kz) = .329*exp(-depth/27.2)+.437*
c     *                exp(-depth/2.93)
c
         irrad(kz) = p1*(1.0-exp(-k1*depth)) +
     *             p2*(1.0-exp(-k2*depth)) +
     *             p3*(1.0-exp(-k3*depth)) +
     *             p4*(1.0-exp(-k4*depth)) 
      enddo
      do kk=kdim,1,-1
        rad(kk) = irrad(kk-1)-irrad(kk)
        if(iope) then
          write(*,*)"Rad ",kk, rad(kk)
        endif
      enddo
#if defined(MPOND)
      do kz=0,kdim
        depth = (kdim+kz)*dz 
        irrad(kz) = p1*(1.0-exp(-k1*depth)) +
     *             p2*(1.0-exp(-k2*depth)) +
     *             p3*(1.0-exp(-k3*depth)) +
     *             p4*(1.0-exp(-k4*depth))
      enddo
      do kk=1,kdim
        rad(kk) = rad(kk) + 0.7*(irrad(kk)-irrad(kk-1))
        if(iope) then
          write(*,*)"Radu ",kk, rad(kk)
        endif
      enddo
#endif
 
c
      return
      end
