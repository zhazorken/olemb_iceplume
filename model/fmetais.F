ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Title:       metais
c  Author:      Eric D. Skyllingstad, Don W. Denbo
c  Date:        5/19/95, 4/24/97
c  Version:     1.0 f77
c  Description:
c
c  Subgrid parameterization based on the method of Metais and Lesieur (1992)
c  modified by Ducros et al. (1996) (JFM, 326, pp.1-36)
c  Uses structure function with  low freq. flow filtered out.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine metais(n)
c
c   $Id: fmetais.F,v 1.4 1999/08/16 20:22:31 eds Exp $
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
      real duw(jdim),dvw(jdim),dww(jdim),due(jdim),dve(jdim)
      real dwe(jdim)
      real uso,uno,vso,vno,wso,wno
      real ws1,kmm
      real vt(0:idim1,0:jdim1,0:kdim1),wt(0:idim1,0:jdim1,0:kdim1)
      real ut(0:idim1,0:jdim1,0:kdim1)
      real uc(0:idim1,0:jdim1,0:kdim1),vc(0:idim1,0:jdim1,0:kdim1)
      real wc(0:idim1,0:jdim1,0:kdim1)
c
c  Setup constants
c  ds=dx=dy=dz  (not valid if dx, dy, and dz are not equal)
c
c first filter the velocities three times using the laplacian
c
      ncoun = 3
      call laplace(um,vm,wm,ut,vt,wt,ncoun)
c
c next compute average velocities at the diffusion grid point
c
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
            uc(i,j,k) = (ut(i+1,j,k)+ut(i,j,k))*0.5
            vc(i,j,k) = (vt(i,j+1,k)+vt(i,j,k))*0.5
            wc(i,j,k) = (wt(i,j,k+1)+wt(i,j,k))*0.5
          enddo
        enddo
      enddo
c
      call imag(uc)
      call imag(vc)
      call imag(wc)

c
c now march across grid and compute modulus
c
      do k=1,kdim
        do j=1,jdim
c
c initialize western differences
c
          duw(j) = (uc(1,j,k)-uc(0,j,k))**2
          dvw(j) = (vc(1,j,k)-vc(0,j,k))**2
          dww(j) = (wc(1,j,k)-wc(0,j,k))**2
        enddo
c
        kmm =0.2*dz*dz/delt
        do i=1,idim
c
c initialize north south values
c
          uno =  (uc(i,1,k)-uc(i,0,k))**2
          vno =  (vc(i,1,k)-vc(i,0,k))**2
          wno =  (wc(i,1,k)-wc(i,0,k))**2
c
          do j=1,jdim
c
c update east values
c
c
            due(j) = duw(j)
            dve(j) = dvw(j)
            dwe(j) = dww(j)
c
c update west values
c
            duw(j) = (uc(i+1,j,k)-uc(i,j,k))**2
            dvw(j) = (vc(i+1,j,k)-vc(i,j,k))**2
            dww(j) = (wc(i+1,j,k)-wc(i,j,k))**2
c
c update south values
c
            uso = uno
            vso = vno
            wso = wno
c
c do new north differences
c
            uno = (uc(i,j+1,k)-uc(i,j,k))**2
            vno = (vc(i,j+1,k)-vc(i,j,k))**2
            wno = (wc(i,j+1,k)-wc(i,j,k))**2
c
c now compute modulus
c
            ws1 = 0.25*(uso+uno+vso+vno+wso+wno+
     &            duw(j)+dvw(j)+dww(j)+due(j)+dve(j)+dwe(j))
c
c compute eddy viscosities and dissipation rate (stored in tke)
c Assume a Prandtl number of 0.6
c
            kmom(i,j,k) = 0.000845*dy*sqrt(ws1)
            kmom(i,j,k) = min(kmom(i,j,k),kmm)
            kheat(i,j,k) = 1.6667*kmom(i,j,k)
c
c
c
          enddo
        enddo
      enddo
c
c increase surface km to similarity value, assuming constant flux
c between first point and surface
c
      do i=1,idim
        do j=1,jdim
#if defined(HBOT)
c          olev = kdist(i,j)
          olev = kdist(i,j)
          if(olev.lt.kdim) then
            kmom(i,j,olev)=kmom(i,j,olev)+0.2*uistarh(i,j)*0.2*(0.5*dz)
            kmom(i,j,olev)=min(kmom(i,j,olev),kmm)
c            kheat(i,j,olev) = 1.6667*kmom(i,j,olev)
            if(olev.gt.1) then
              kmom(i,j,olev-1) = kmom(i,j,olev-1) + 
     *          0.2*uistarh(i,j)*0.2*(0.5*dz+dz)
                kmom(i,j,olev-1) = min(kmom(i,j,olev-1),kmm)
c              kheat(i,j,olev-1) = 1.6667*kmom(i,j,olev-1)
           endif
          endif

#endif
#if defined(BOTTOM) 
          kmom(i,j,1) = kmom(i,j,1) + 0.2*uistar(i,j)*0.2*dz
          kmom(i,j,1) = min(kmom(i,j,1),kmm)
          kmom(i,j,2) = kmom(i,j,2) + 0.2*uistar(i,j)*0.2*dz*2.
          kmom(i,j,2) = min(kmom(i,j,2),kmm)
          kheat(i,j,1) = 1.6667*kmom(i,j,1)
c          kheat(i,j,2) = 1.6667*kmom(i,j,2)
#endif
#if !defined(RADTOP) && !defined(ICE)
          kmom(i,j,kdim) = kmom(i,j,kdim) + 0.2*uistart(i,j)*0.2*dz
          kheat(i,j,kdim) = 1.6667*kmom(i,j,kdim)
#endif
#if defined(ICE)
          olev = kdist(i,j)
          kmom(i,j,olev) = kmom(i,j,olev) + 0.2*uistar(i,j)*0.2*dz
          kheat(i,j,olev) = 1.6667*kmom(i,j,olev)
#endif
        enddo
      enddo
#if defined(OPENBC)
#if !defined(UNDERCUT)
c      if(iblkstart(my_pe+1).le.1) then

        do k=1,kdim
c          if(k.gt.18) then
c          do j=1,jdim
c           kmom(19,j,k) = kmom(19,j,k)+0.2*sqrt(abs(dragw(j,k)))*0.2*dxx
c          kmom(19,j,k) = min(kmom(19,j,k),kmm)
c           kmom(20,j,k) = kmom(20,j,k)+0.2*sqrt(abs(dragw(j,k)))*0.2*dxx
c          kmom(20,j,k) = min(kmom(20,j,k),kmm)
c          enddo
cc          else
          do j=1,jdim
          do i=1,idim
            if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
             kmom(i,j,k)=kmom(i,j,k)+0.2*sqrt(abs(dragw(j,k)))*0.2*dxx
             kmom(i,j,k)=kmom(i,j,k)+0.2*sqrt(abs(dragvw(j,k)))*0.2*dxx
             kmom(i,j,k) = min(kmom(i,j,k),kmm)
             kheat(i,j,k) = 1.6667*kmom(i,j,k)
            endif
          enddo
          enddo
        enddo
#endif
#endif

c
c do top/bottom boundaries
c
      call topbotb(kmom,MDKM,n)
      call topbotb(kheat,MDKH,n)
c
c  return to calling routine:
c
      return
      end
c
c laplacian subroutine
c
      subroutine laplace(ui,vi,wi,un,vn,wn,n)
#include "param.inc"
c
c This routine does a laplacian filter n times
c
      real ui(0:idim1,0:jdim1,0:kdim1),vi(0:idim1,0:jdim1,0:kdim1)
      real wi(0:idim1,0:jdim1,0:kdim1),uo(0:idim1,0:jdim1,0:kdim1)
      real vo(0:idim1,0:jdim1,0:kdim1),wo(0:idim1,0:jdim1,0:kdim1)
      real vn(0:idim1,0:jdim1,0:kdim1),wn(0:idim1,0:jdim1,0:kdim1)
      real un(0:idim1,0:jdim1,0:kdim1)
      integer n,nn,i,j,k

c first copy over unfiltered values
c
      do k=0, kdim1
        do j=0,jdim1
          do i=0,idim1
            uo(i,j,k) = ui(i,j,k)
            vo(i,j,k) = vi(i,j,k)
            wo(i,j,k) = wi(i,j,k)
          enddo
        enddo
      enddo

c
c compute filtered fields

      do nn=1,n
        do k=0,kdim1
          do j=1,jdim
            do i=1,idim
              un(i,j,k)  = uo(i+1,j,k)-4.*uo(i,j,k)+uo(i-1,j,k) +
     *                     uo(i,j+1,k)+uo(i,j-1,k)
              vn(i,j,k)  = vo(i+1,j,k)-4.*vo(i,j,k)+vo(i-1,j,k) +
     *                     vo(i,j+1,k)+vo(i,j-1,k)
              wn(i,j,k)  = wo(i+1,j,k)-4.*wo(i,j,k)+wo(i-1,j,k) +
     *                     wo(i,j+1,k)+wo(i,j-1,k)
            enddo
          enddo
        enddo
c
c update boundaries 
c
        call imag(un)
        call imag(vn)
        call imag(wn)
c
c copy to old fields for next pass
c
        if(nn.ne.n) then
          do k=0, kdim1 
            do j=0,jdim1 
              do i=0,idim1 
                uo(i,j,k) = un(i,j,k)
                vo(i,j,k) = vn(i,j,k)
                wo(i,j,k) = wn(i,j,k)
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
