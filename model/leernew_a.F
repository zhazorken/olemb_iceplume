      subroutine leerftop(n,phi,phip,ivar)
      integer n,ivar
c
c   $Id: leer.F,v 1.6 1998/10/21 21:56:39 dwd Exp $
c
c  **  This subroutine calculates a positive, definite
c  **  horizontal and vertical advection based on the
c  **  a corrected upstream scheme.  Input arrays are:
c  **
c  **     n   - time step
c  **     phi - input old scalar values, must be positive
c  **     phip - new output scalar values
c  **     ivar - variable type
c
c  ** common blocks and variable types may need redefinition
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "olemtke.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif
#include "communicate.inc"
#include "finitv.inc"
c
      integer i,j,k,ip2,jp2,im2,jm2,kp2
      real phip(0:idim1,0:jdim1,0:kdim1),phi(0:idim1,0:jdim1,0:kdim1)
#if defined(USEMPI)
      integer ier
      real pnorth(0:idim1,0:kdim1)
      real psouth(0:idim1,0:kdim1)
      real aynorth(0:idim1,0:kdim1)
      real aysouth(0:idim1,0:kdim1)
      real peast(0:jdim1,0:kdim1)
      real pwest(0:jdim1,0:kdim1)
      real axeast(0:jdim1,0:kdim1)
      real axwest(0:jdim1,0:kdim1)

#endif
      real phix(0:idim1,0:jdim1,0:kdim1)
      real phiy(0:idim1,0:jdim1,0:kdim1)
      real phiz(0:idim1,0:jdim1,0:kdim1)
      real wsb(0:idim1,0:jdim1,0:kdim1)
      real phim,delc,fluxu1,fluxu0,fluxw1,fluxw0,fluxvb
      real fluxwd,fluxud,fluxub,fluxwb,ub,vb,wb,ucross,wcross,vcross
      real sgs,eps,dphi
      real kright,kleft,kback,kfront,ktop,kbot,b1,b2,b3
      integer olev
c
#include "moddef.inc"
c
c first do the v term
c
c compute bootstrap flux on edge if j equals 1
c
c get a copy of the north boundary second level boundary points
c

c small number
      eps = 1.0e-7
#if defined(USEMPI) 
      call northshift(phi,pnorth) 
      call MPI_BARRIER(comm,ier)
      call southshift(phi,psouth)
      call MPI_BARRIER(comm,ier)
      call northshift(ay,aynorth) 
      call MPI_BARRIER(comm,ier)
      call southshift(ay,aysouth) 
      call MPI_BARRIER(comm,ier)
      call westshift(phi,pwest)
      call MPI_BARRIER(comm,ier)
      call eastshift(phi,peast)
      call MPI_BARRIER(comm,ier)
      call westshift(ax,axwest)
      call MPI_BARRIER(comm,ier)
      call eastshift(ax,axeast)
      call MPI_BARRIER(comm,ier)


#endif
c
c define difference arrays
c
      do k=1,kdim
        do j=1,jdim
          do i=1,idim
            phix(i,j,k) = (phi(i+1,j,k)-phi(i,j,k))
     *                    *ax(i+1,j,k)*ax(i,j,k)
            phiy(i,j,k) = (phi(i,j+1,k)-phi(i,j,k))
     *                    *ay(i,j+1,k)*ay(i,j,k)
           
#if defined(ICE)
            if(k.ge.kdist(i,j)) then
              phiz(i,j,k) = 0.0
            else
              phiz(i,j,k) = (phi(i,j,k+1)-phi(i,j,k))
            endif
#endif
#if defined(HBOT)
            if(k.ge.kdist(i,j)) then
              phiz(i,j,k) = 0.0
            else
              phiz(i,j,k) = (phi(i,j,k+1)-phi(i,j,k))
            endif
#endif
#if !defined(HBOT) && !defined(ICE)
             phiz(i,j,k) = (phi(i,j,k+1)-phi(i,j,k))
#endif
          enddo
        enddo
      enddo
      call imag(phix)
      call imag(phiy)
      call imag(phiz)
      do j=0,jdim1
        do i=0,idim1
          phix(i,j,kdim1) = 0.0
          phiy(i,j,kdim1) = 0.0
          phiz(i,j,kdim1) = 0.0
          phiz(i,j,kdim) = 0.0
          phix(i,j,0) = 0.0
          phiy(i,j,0) = 0.0
          phiz(i,j,0) = 0.0
        enddo
      enddo
      j=1
      do k=1,kdim
         do i=1,idim
c
c
            delc = 0.0
            wb = (w(i,j-1,k)+w(i,j,k)+w(i,j,k+1)+w(i,j-1,k+1))
     *           *0.25
            if(wb.gt.0) then
               wcross = 0.5*dtdz*wb*phiz(i,j,k-1)
            else
               if(k.eq.kdim) then
                  wcross = 0.0
               else
                  wcross = 0.5*dtdz*wb*phiz(i,j,k)
               endif
            endif
c            wcross = 0.0
            ub = (u(i,j-1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))
     *           *0.25
            if(ub.gt.0) then
               ucross = 0.5*dtdx*ub*phix(i-1,j,k)
            else
               ucross = 0.5*dtdx*ub*phix(i,j,k)
            endif
c            ucross = 0.0
            if(v(i,j,k).gt.0) then
               jm2 = jdim-1
#if defined(USEMPI)
               phim=(phi(i,j-1,k)-psouth(i,k))*phiy(i,j-1,k)
     *              *ay(i,j-1,k)*aysouth(i,k)
#else
               phim=(phi(i,j-1,k)-phi(i,jm2,k))*phiy(i,j-1,k)
     *               *ay(i,jm2,k)*ay(i,j-1,k)
#endif
#if defined(USEMPI)
               dphi = (phi(i,j,k)-psouth(i,k))
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
#else
               dphi = phi(i,j,k)-phi(i,jm2,k)
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
#endif
               endif
               fluxvb = v(i,j,k)*(phi(i,j-1,k)+0.5*(1.-v(i,j,k)*dtdy)*
     *              delc - wcross - ucross)*ay(i,j,k)
            else
               phim=phiy(i,j-1,k)*phiy(i,j,k)
               dphi= (phi(i,j+1,k)-phi(i,j-1,k))
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
               endif
               fluxvb = v(i,j,k)*(phi(i,j,k)-0.5*(1.+v(i,j,k)*dtdy)*
     *              delc - wcross - ucross)*ay(i,j,k)
            endif
            
            fluxvd(i,k) = fluxvb
         enddo
      enddo
c
      do j=1,jdim
         do k=1,kdim
            do i=1,idim
c
c do horizontal terms first
c
c
c first, the upstream flux
c
               delc = 0.0
               wb = (w(i,j+1,k)+w(i,j,k)+w(i,j,k+1)+w(i,j+1,k+1))
     *              *0.25
               if(wb.gt.0) then
                  wcross = 0.5*dtdz*wb*phiz(i,j+1,k-1)
               else
                  if(k.eq.kdim) then
                     wcross = 0.0
                  else
                     wcross = 0.5*dtdz*wb*phiz(i,j+1,k)
                  endif
               endif
c               wcross = 0.0
               ub = (u(i,j+1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j+1,k))
     *              *0.25
               if(ub.gt.0) then
                  ucross = 0.5*dtdx*ub*phix(i-1,j+1,k)
               else
                  ucross = 0.5*dtdx*ub*phix(i,j+1,k)
               endif
c               ucross = 0.0
               if(v(i,j+1,k).gt.0) then
                  phim=phiy(i,j-1,k)*phiy(i,j,k)
                  dphi = (phi(i,j+1,k)-phi(i,j-1,k))
                  if(phim.gt.0.and.dphi.ne.0.0)then
                     delc = 2.*phim/dphi
                  endif
                  fluxv1(i,k) = ay(i,j+1,k)*v(i,j+1,k)*(phi(i,j,k)+
     *                 0.5*(1.-v(i,j+1,k)*dtdy)*delc - wcross - ucross)
               else
                  jp2 = j+2
                  
#if defined(USEMPI)
                  if(j.eq.jdim)then
                    phim=phiy(i,j,k)*(pnorth(i,k)- 
     *                    phi(i,j+1,k))*ay(i,j+1,k)*aynorth(i,k)
                    dphi = (pnorth(i,k)-phi(i,j,k)) 
                    if(phim.gt.0.and.dphi.ne.0.0)then 
                       delc = 2.*phim/dphi
                    endif
                 else
                    phim=phiy(i,j,k)*(phi(i,jp2,k)-
     *                   phi(i,j+1,k))*ay(i,jp2,k)*ay(i,j+1,k)
                    dphi = (phi(i,jp2,k)-phi(i,j,k))
                    if(phim.gt.0.and.dphi.ne.0.0)then
                       delc = 2.*phim/dphi
                    endif
                 endif
#else 
                 if(j.eq.jdim)jp2 = 2
                 phim=phiy(i,j,k)*(phi(i,jp2,k)- phi(i,j+1,k))
     *                   *ay(i,jp2,k)*ay(i,j+1,k)
                 dphi = (phi(i,jp2,k)-phi(i,j,k))
                 if(phim.gt.0.and.dphi.ne.0.0)then
                    delc = 2.*phim/dphi
                 endif
#endif
                 fluxv1(i,k) = ay(i,j+1,k)*v(i,j+1,k)*(phi(i,j+1,k)-
     *                0.5*(1.+v(i,j+1,k)*dtdy)*delc - wcross - ucross)
              endif
c              if(k.gt.kdist(i,j)) then
                 fluxv1(i,k) = 0.0
c              endif
c
c next the downstream flux
c
              fluxv0(i,k) = fluxvd(i,k)
              fluxvd(i,k) = fluxv1(i,k)
c
c compute v flux
c
                phip(i,j,k)=phi(i,j,k)-dtdy*(fluxv1(i,k)-fluxv0(i,k))
c     
           enddo
        enddo
      enddo


c     
      do j=1,jdim
         do k=1,kdim 
c     
c compute bootstrap flux on edge
c
            i = 1
            delc = 0.0
            wb = (w(i-1,j,k)+w(i,j,k)+w(i,j,k+1)+w(i-1,j,k+1))
     *           *0.25
            if(wb.gt.0) then
               wcross = 0.5*dtdz*wb*phiz(i,j,k-1)
            else
               if(k.eq.kdim) then
                  wcross = 0.0
               else
                  wcross = 0.5*dtdz*wb*phiz(i,j,k)
               endif
            endif
c            wcross = 0.0
            vb = (v(i-1,j,k)+v(i,j,k)+v(i,j+1,k)+v(i-1,j+1,k))
     *           *0.25
            if(vb.gt.0) then
               vcross = 0.5*dtdy*vb*phiy(i,j-1,k)
            else
               vcross = 0.5*dtdy*vb*phiy(i,j,k)
            endif
c            vcross = 0.0
c
            if(u(i,j,k).gt.0) then
               im2 = idim-1
               phim=(phi(i-1,j,k)-pwest(j,k))*phix(i-1,j,k)
     *              *axwest(j,k)*ax(i-1,j,k)
               dphi = (phi(i,j,k)-pwest(j,k))
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
               endif
               fluxub = u(i,j,k)*(phi(i-1,j,k)
     *              +0.5*(1.-u(i,j,k)*dtdx)*
     *              delc - wcross - vcross)*ax(i,j,k)
            else
               phim=phix(i-1,j,k)*phix(i,j,k)
               dphi = (phi(i+1,j,k)-phi(i-1,j,k))
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
               endif
               fluxub = u(i,j,k)*(phi(i,j,k)-0.5*(1.+u(i,j,k)*dtdx)*
     *              delc - wcross - vcross)*ax(i,j,k)
            endif

            fluxud = fluxub

            do i=1,idim 
c
c do horizontal terms first
c
c
c first, the upstream flux
c
               delc = 0.0
               wb = (w(i+1,j,k)+w(i,j,k)+w(i,j,k+1)+w(i+1,j,k+1))
     *              *0.25
               if(wb.gt.0) then
                  wcross = 0.5*dtdz*wb*phiz(i+1,j,k-1)
               else
                  if(k.eq.kdim) then
                     wcross = 0.0
                  else
                     wcross = 0.5*dtdz*wb*phiz(i+1,j,k)
                  endif
               endif
c               wcross = 0.0
               vb = (v(i+1,j,k)+v(i,j,k)+v(i,j+1,k)+v(i+1,j+1,k))
     *              *0.25   
               if(vb.gt.0) then
                  vcross = 0.5*dtdy*vb*phiy(i+1,j-1,k)
               else
                  vcross = 0.5*dtdy*vb*phiy(i+1,j,k)
               endif
c               vcross = 0.0
               if(u(i+1,j,k).gt.0) then
                  phim=phix(i-1,j,k)*phix(i,j,k)
                  dphi = (phi(i+1,j,k)-phi(i-1,j,k))
                  if(phim.gt.0.and.dphi.ne.0.0)then
                    if(phi(i+1,j,k).eq.phi(i-1,j,k)) then
                      write(*,*) "phi diff zero"
                      write(*,*) phim,idim
                      write(*,*) phix(i-1,j,k),phix(i,j,k)
                      write(*,*) phi(i+1,j,k),phi(i,j,k),phi(i-1,j,k)
                      write(*,*) ax(i+1,j,k),ax(i,j,k),ax(i-1,j,k)
                      write(*,*) iblkstart(my_pe+1),i,j,k
                    endif
                     delc = 2.*phim/(phi(i+1,j,k)-phi(i-1,j,k))
                  endif
                  fluxu1 = ax(i+1,j,k)*u(i+1,j,k)*(phi(i,j,k)+
     *                 0.5*(1.-u(i+1,j,k)*dtdx)*delc - wcross - vcross)
               else
                  ip2 = i+2
                  if(i.eq.idim)then
                   phim=phix(i,j,k)*(peast(j,k)-
     *                 phi(i+1,j,k))*ax(i+1,j,k)*axeast(j,k)
                    dphi = (peast(j,k)-phi(i,j,k))
                    if(phim.gt.0.and.dphi.ne.0.0)then
                       delc = 2.*phim/dphi
                    endif
                    fluxu1 = ax(i+1,j,k)*u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
                  else
                    phim=phix(i,j,k)*(phi(ip2,j,k)-
     *                 phi(i+1,j,k))*ax(ip2,j,k)*ax(i+1,j,k)
                    dphi = (phi(ip2,j,k)-phi(i,j,k))
                    if(phim.gt.0.and.dphi.ne.0.0)then
                       delc = 2.*phim/dphi
                    endif
                    fluxu1 = ax(i+1,j,k)*u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
                  endif

c                  if(i.eq.idim)ip2 = 2
c                  phim=phix(i,j,k)*(phi(ip2,j,k)- phi(i+1,j,k))
c     *                 *ax(ip2,j,k)
c                  if(phim.gt.0)then 
c                     delc = 2.*phim/(phi(ip2,j,k)-phi(i,j,k)) 
c                  endif 
c                  fluxu1 = ax(i+1,j,k)*u(i+1,j,k)*(phi(i+1,j,k)-
c     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
               endif
c               if(k.gt.kdist(i,j)) then
                 fluxu1 = 0.0
c               endif
c     
c     next the downstream flux
c     
               fluxu0 = fluxud
               fluxud = fluxu1
c
c compute u flux
c
                phip(i,j,k)=phip(i,j,k)-dtdx*(fluxu1-fluxu0)
c


            enddo
         enddo

c
c Now  do the vertical terms
c
         do i=1,idim
c
c compute the bottom boot strap
c
            k=1
            delc = 0.0
            ub = (u(i+1,j,k-1)+u(i,j,k-1)+u(i,j,k)+u(i+1,j,k))*0.25
            if(ub.gt.0) then
               ucross = 0.5*dtdx*ub*phix(i-1,j,k)
            else   
               ucross = 0.5*dtdx*ub*phix(i,j,k)
            endif  
c            ucross = 0.0
            vb = (v(i,j+1,k-1)+v(i,j,k-1)+v(i,j,k)+v(i,j+1,k))*0.25
            if(vb.gt.0) then
               vcross = 0.5*dtdy*vb*phiy(i,j-1,k)
            else  
               vcross = 0.5*dtdy*vb*phiy(i,j,k)
            endif
c            vcross = 0.0
c
c
            if(w(i,j,k).gt.0) then
               phim=phiz(i,j,k-1)*phiz(i,j,k)
               dphi = (phi(i,j,k+1)-phi(i,j,k-1))
               if(phim.gt.0.and.dphi.ne.0.0) then
                  delc = 2.*phim/dphi
               endif
               fluxwb = w(i,j,k)*(phi(i,j,k-1) + 0.5*(1.-w(i,j,k)*dtdz)*
     *              delc - ucross - vcross)
            else
               phim=phiz(i,j,k-1)*phiz(i,j,k)
               dphi = (phi(i,j,k+1)-phi(i,j,k-1))
               if(phim.gt.0.and.dphi.ne.0.0)then
                  delc = 2.*phim/dphi
               endif
c
c do upstream
c
               fluxwb = w(i,j,k)*(phi(i,j,k)-0.5*(1.+w(i,j,k)*dtdz)*
     *              delc - ucross - vcross)
            endif
#if defined(HBOT)
            if(k.gt.kdist(i,j)) then
              fluxwb = 0.0
            endif
#endif
            fluxwd = fluxwb
c
c for k=kdim the following is degenerate
c where phip(i,j,k) = phip(i,j,k) + dtdz*fluxw0
c
            do k=1,kdim
c
c first the upstream flux
c
               delc = 0.0
               ub = (u(i+1,j,k)+u(i,j,k)+u(i,j,k+1)+u(i+1,j,k+1))
     *              *0.25
               if(ub.gt.0) then 
                  ucross = 0.5*dtdx*ub*phix(i-1,j,k+1)
               else 
                  ucross = 0.5*dtdx*ub*phix(i,j,k+1)
               endif
c               ucross = 0.0
               vb = (v(i,j+1,k)+v(i,j,k)+v(i,j,k+1)+v(i,j+1,k+1))
     *              *0.25
               if(vb.gt.0) then
                  vcross = 0.5*dtdy*vb*phiy(i,j-1,k+1)
               else
                  vcross = 0.5*dtdy*vb*phiy(i,j,k+1)
               endif
c               vcross = 0.0
               if(w(i,j,k+1).gt.0) then
                  phim=(phiz(i,j,k-1))*phiz(i,j,k)
                  dphi = (phi(i,j,k+1)-phi(i,j,k-1))
                  if(phim.gt.0.and.dphi.ne.0.0) then
                     delc = 2.*phim/dphi
                  endif
                  fluxw1 = w(i,j,k+1)*(phi(i,j,k)+0.5*
     *                 (1.-w(i,j,k+1)*dtdz)*delc - ucross - vcross)
               else
                  kp2 = k+2
#if defined(HBOT)
                  if(kp2.gt.kdist(i,j)) kp2=kdist(i,j)
#else
                  if(kp2.gt.kdim) kp2=kdim
#endif
                  phim=phiz(i,j,k)*(phi(i,j,kp2)-phi(i,j,k+1))
     *                 *way(i,j,kp2)
                  dphi = (phi(i,j,kp2)-phi(i,j,k))
                  if(phim.gt.0.and.dphi.ne.0.0)then
                     delc = 2.*phim/dphi
                  endif
                  fluxw1 = w(i,j,k+1)*(phi(i,j,k+1)-0.5*
     *                 (1.+w(i,j,k+1)*dtdz)*delc - ucross - vcross)
               endif
#if defined(HBOT)
               if(k.ge.kdist(i,j)) then
                 fluxw1 = 0.0
               endif
#endif

c     
c next the downstream flux
c
               fluxw0 = fluxwd
               fluxwd = fluxw1
c
               if(k.le.kdist(i,j)) then
                 phip(i,j,k)=phip(i,j,k)-dtdz*(fluxw1-fluxw0)
               endif
c     
            enddo
         enddo
      enddo

c
c  **  Call boundary update routine using the new phip
c      no sgs call.
c
c      call image(phip)
#if defined(USEMPI)
      call MPI_BARRIER(comm,ier)
#endif
      call topbotb(phip,MDNULL,n)
c
c  **  This call to bcupdte updates phip at 0 and idim1 and o and kdim1
c  **  The top and bottom image points are probably not used at this
c  **  stage.
c
c
c  **  do compressiblity correction
c
      do k=1,kdim
         do j=1,jdim
            do i=1,idim
c               if(k.le.kdist(i,j)) then
               phip(i,j,k)=phip(i,j,k)+((ax(i+1,j,k)*u(i+1,j,k)
     *              -ax(i,j,k)*u(i,j,k))*dtdx+
     *              (w(i,j,k+1)-w(i,j,k))*dtdz
     *              + (ay(i,j+1,k)*v(i,j+1,k)-ay(i,j,k)*v(i,j,k))*dtdy)
     *              *phi(i,j,k)
c               endif
            enddo
         enddo
      enddo

c
c finish up the processing by mixing through turbulence parameterization.
c Don't do this from the TKE fields.  TKE is diffused in the sgs update
c routine.
c
      if(ivar.eq.MDTEMP)then
          do j=1,jdim
            do i=1,idim
c simple box
#if defined(HBOT)
               do k=1,kdist(i,j)
#else
               do k=1,kdim
#endif
                  kright = kheat(i+1,j,k)+kheat(i,j,k)+2.*khb
                  kleft = kheat(i,j,k)+kheat(i-1,j,k)+2.*khb
                  b1=(kright*phix(i,j,k)-
     *                 kleft*phix(i-1,j,k))*0.5*odx*odx
c
                  kback = kheat(i,j+1,k)+kheat(i,j,k)+2.*khb
                  kfront = kheat(i,j,k)+kheat(i,j-1,k)+2.*khb
                  b2=(kback*phiy(i,j,k)-
     *                 kfront*phiy(i,j-1,k))*0.5*ody*ody
c
                  ktop = kheat(i,j,k+1)+kheat(i,j,k)+2.*khb
#if defined(HBOT)
                  if(k.ge.kdist(i,j)) then
                    ktop = 0.0
                  endif
#endif
                  sgs=b1+b2
                  phip(i,j,k)=phip(i,j,k)+delt*sgs
c
                  if(k.eq.kdim) then
                    ktop = 0.0
                  endif
                  wsb(i,j,k) = ktop*(phi(i,j,k+1)-phi(i,j,k))
     *                *0.5*odz*odz
                  sgs=b1+b2
                  phip(i,j,k)=phip(i,j,k)+delt*sgs
                  tke(i,j,k) = tke(i,j,k)+b1+b2
               enddo
               wsb(i,j,kdist(i,j)) = 0.0
               wsb(i,j,kdim) = 0.0
            enddo
         enddo
      else if(ivar.eq.MDSAL) then
          do j=1,jdim
            do i=1,idim
c simple box
#if defined(HBOT)
               do k=1,kdist(i,j)
#else
               do k=1,kdim
#endif
                  kright = kheat(i+1,j,k)+kheat(i,j,k) + 2.*kms
                  kleft = kheat(i,j,k)+kheat(i-1,j,k) + 2.*kms
                  b1=(kright*phix(i,j,k)-
     *                 kleft*phix(i-1,j,k))*0.5*odx*odx
c
                  kback = kheat(i,j+1,k)+kheat(i,j,k) + 2.*kms
                  kfront = kheat(i,j,k)+kheat(i,j-1,k) + 2.*kms
                  b2=(kback*phiy(i,j,k)-
     *                 kfront*phiy(i,j-1,k))*0.5*ody*ody
c
                  ktop = kheat(i,j,k+1)+kheat(i,j,k) + 2.*kms
#if defined(HBOT)
                  if(k.ge.kdist(i,j)) then
                    ktop = 0.0
                  endif
#endif
                  if(k.eq.kdim) then
                    ktop = 0.0
                  endif
                  wsb(i,j,k) = ktop*(phi(i,j,k+1)-phi(i,j,k))
     *                *0.5*odz*odz
                  sgs=b1+b2
                  phip(i,j,k)=phip(i,j,k)+delt*sgs
               enddo
               wsb(i,j,kdist(i,j)) = 0.0
               wsb(i,j,1) = 0.0
               wsb(i,j,kdim) = 0.0
            enddo
         enddo
      endif
      do i=1,idim
        do j=1,jdim
#if defined(HBOTOLD) 
          do k=kdist(i,j)+1,kdim
#else
          do k=1,kdist(i,j)
#endif
            phip(i,j,k) = 
     *         phip(i,j,k)+delt*(wsb(i,j,k)-wsb(i,j,k-1))
          enddo
        enddo
      enddo

c
c update boundary values
c
      call topbotb(phip,ivar,n)
      call imag(phip)
c
c
c  **  get out of here normally
c
      return
      end
