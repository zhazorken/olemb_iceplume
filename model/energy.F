ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Title:       energy
c  Author:      Eric D. Skyllingstad
c  Date:        5/19/95, 4/24/97, 1/6/2010
c  Version:     1.0 f77
c  Description:
c
c Routine for calculating energy of small scale and baroclinic scales
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine energy(n)
c
c   $Id: energy.F,v 1.4 1999/08/16 20:22:31 eds Exp $
c
c  include files:

#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "moddef.inc"
#include "cyclic.inc"
#include "finitv.inc"
#include "ice.inc"
#include "communicate.inc"


c  declare variables:
c
      integer i,j,k,n,ncoun,olev
c
      real up(0:idim1,0:jdim1,0:kdim1),vp(0:idim1,0:jdim1,0:kdim1)
      real wp(0:idim1,0:jdim1,0:kdim1)
      real pp(0:idim1,0:jdim1,0:kdim1)
      real rhopp(0:idim1,0:jdim1,0:kdim1)
c
c  Setup constants
c  ds=dx=dy=dz  (not valid if dx, dy, and dz are not equal)
c
c first filter the velocities three times using the laplacian
c
      ncoun = 15 
      call fivept(um,uf,ncoun)
      call fivept(vm,vf,ncoun)
      call fivept(wm,wf,ncoun)
      call fivept(th,rhopf,ncoun)
      call fivept(p,pf,ncoun)
c
c first compute terms for fine scale
c
      do k=0,kdim1
        do j=0,jdim1
          do i=0,jdim1
            up(i,j,k) = um(i,j,k)-uf(i,j,k)
            vp(i,j,k) = vm(i,j,k)-vf(i,j,k)
            wp(i,j,k) = wm(i,j,k)-wf(i,j,k)
            pp(i,j,k) = p(i,j,k)-pf(i,j,k)
c            rhopp(i,j,k) = rhop(i,j,k)-rhopf(i,j,k)
            u2_ave(k) = u2_ave(k)+up(i,j,k)**2
            v2_ave(k) = v2_ave(k)+vp(i,j,k)**2
            w2_ave(k) = w2_ave(k)+wp(i,j,k)**2
          enddo
        enddo
      enddo
c
c compute four shear production terms
c
      do k=1,kdim
      do j=1,jdim
      do i=1,idim
c
c horizontal shear production terms
c
        vpupum(k) = (up(i,j,k)+up(i,j-1,k))*(vp(i,j,k)+vp(i-1,j,k))*
     *           (uf(i,j,k)-uf(i,j-1,k))*ody*0.25 + vpupum(k)
c
        vpupvm(k) = (up(i,j,k)+up(i,j-1,k))*(vp(i,j,k)+vp(i-1,j,k))*
     *           (vf(i,j,k)-vf(i-1,j,k))*ody*0.25 + vpupvm(k)
c
c vertical shear production terms
c
        upwpum(k) = (up(i,j,k)+up(i,j,k-1))*(wp(i,j,k)+wp(i-1,j,k))*
     *           (uf(i,j,k)-uf(i,j,k-1))*odz*0.25 + upwpum(k)
c
        vpwpvm(k) = (vp(i,j,k)+vp(i,j,k-1))*(wp(i,j,k)+wp(i,j-1,k))*
     *           (vf(i,j,k)-vf(i,j,k-1))*odz*0.25 + vpwpvm(k)
c
c buoyancy
c
        wrhop(k) = wp(i,j,k)*(rhopp(i,j,k)+rhopp(i,j,k-1))*grav/
     *              rhob(k)+ wrhop(k)
c
c        tke(i,j,k) = vpupum(k)+vpupvm(k)+upwpum(k)+vpwpvm(k)
      enddo
      enddo
      enddo
      
c
      return
      end
c
c filter subroutine
c
      subroutine fivept(ui,un,n)
#include "param.inc"
c
c This routine does a five point average  filter n times
c
      real ui(0:idim1,0:jdim1,0:kdim1),un(0:idim1,0:jdim1,0:kdim1)
      real uo(0:idim1,0:jdim1,0:kdim1)
      integer n,nn,i,j,k

c first copy over unfiltered values
c
      do k=0, kdim1
        do j=0,jdim1
          do i=0,idim1
            uo(i,j,k) = ui(i,j,k)
          enddo
        enddo
      enddo

c
c compute filtered fields

      do nn=1,n
        do k=0,kdim1
          do j=1,jdim
            do i=1,idim
              un(i,j,k)  = (uo(i+1,j,k)+uo(i,j,k)+uo(i-1,j,k) +
     *                     uo(i+1,j+1,k)+uo(i-1,j-1,k)+
     *                     uo(i+1,j-1,k)+uo(i-1,j+1,k)+
     *                     uo(i,j+1,k)+uo(i,j-1,k))/9.0
            enddo
          enddo
        enddo
c
c update boundaries 
c
        call imag(un)
c
c copy to old fields for next pass
c
        if(nn.ne.n) then
          do k=0, kdim1 
            do j=0,jdim1 
              do i=0,idim1 
                uo(i,j,k) = un(i,j,k)
              enddo
            enddo 
          enddo
        endif
      enddo
c
c all done
c
      return
      end
