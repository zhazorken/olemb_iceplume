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
#include "finitv.inc"
#include "communicate.inc"
c
      integer i,j,k,ip2,jp2,im2,jm2,kp2
      real phip(0:idim1,0:jdim1,0:kdim1),phi(0:idim1,0:jdim1,0:kdim1)
#if defined(USEMPI)
      integer ier
      real pnorth(0:idim1,0:kdim1)
      real psouth(0:idim1,0:kdim1)
#endif
      real phim,delc,fluxu1,fluxu0,fluxw1,fluxw0,fluxvb
      real fluxwd,fluxud,fluxub,fluxwb,ub,vb,wb,ucross,wcross,vcross
      real sgs
      real kright,kleft,kback,kfront,ktop,kbot,b1,b2,b3
c
#include "moddef.inc"
c
c first do the v term
c
c compute bootstrap flux on edge if j equals 1
c
c get a copy of the north boundary second level boundary points
c
#if defined(USEMPI) 
      call northshift(phi,pnorth) 
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
      call southshift(phi,psouth)
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif
      j=1
      do k=1,kdim
         do i=1,idim
c
c
            delc = 0.0
            wb = (w(i,j-1,k)+w(i,j,k)+w(i,j,k+1)+w(i,j-1,k+1))
     *           *0.25
            if(wb.gt.0) then
               wcross = 0.5*dtdz*wb*(phi(i,j,k)-phi(i,j,k-1))
            else
               if(k.eq.kdim) then
                  wcross = 0.0
               else
                  wcross = 0.5*dtdz*wb*(phi(i,j,k+1)-phi(i,j,k))
               endif
            endif
            ub = (u(i,j-1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))
     *           *0.25
            if(ub.gt.0) then
               ucross = 0.5*dtdx*ub*(phi(i,j,k)-phi(i-1,j,k))
            else
               ucross = 0.5*dtdx*ub*(phi(i+1,j,k)-phi(i,j,k))
            endif
            if(v(i,j,k).gt.0) then
               jm2 = jdim-1
#if defined(USEMPI)
               phim=(phi(i,j-1,k)-psouth(i,k))*(phi(i,j,k)
#else
               phim=(phi(i,j-1,k)-phi(i,jm2,k))*(phi(i,j,k)
#endif
     *              -phi(i,j-1,k))
               if(phim.gt.0)then
#if defined(USEMPI)
                  delc = 2.*phim/(phi(i,j,k)-psouth(i,k))
#else
                  delc = 2.*phim/(phi(i,j,k)-phi(i,jm2,k))
#endif
               endif
               fluxvb = v(i,j,k)*(phi(i,j-1,k)+0.5*(1.-v(i,j,k)*dtdy)*
     *              delc - wcross - ucross)
            else
               phim=(phi(i,j,k)-phi(i,j-1,k))*(phi(i,j+1,k)-
     *              phi(i,j,k))
               if(phim.gt.0)then
                  delc = 2.*phim/(phi(i,j+1,k)-phi(i,j-1,k))
               endif
               fluxvb = v(i,j,k)*(phi(i,j,k)-0.5*(1.+v(i,j,k)*dtdy)*
     *              delc - wcross - ucross)
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
                  wcross = 0.5*dtdz*wb*(phi(i,j+1,k)-phi(i,j+1,k-1))
               else
                  if(k.eq.kdim) then
                     wcross = 0.0
                  else
                     wcross = 0.5*dtdz*wb*(phi(i,j+1,k+1)-phi(i,j+1,k))
                  endif
               endif
               ub = (u(i,j+1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j+1,k))
     *              *0.25
               if(ub.gt.0) then
                  ucross = 0.5*dtdx*ub*(phi(i,j+1,k)-phi(i-1,j+1,k))
               else
                  ucross = 0.5*dtdx*ub*(phi(i+1,j+1,k)-phi(i,j+1,k))
               endif
               if(v(i,j+1,k).gt.0) then
                  phim=(phi(i,j,k)-phi(i,j-1,k))*(phi(i,j+1,k)-
     *                 phi(i,j,k))
                  if(phim.gt.0)then
                     delc = 2.*phim/(phi(i,j+1,k)-phi(i,j-1,k))
                  endif
                  fluxv1(i,k) = v(i,j+1,k)*(phi(i,j,k)+
     *                 0.5*(1.-v(i,j+1,k)*dtdy)*delc - wcross - ucross)
               else
                  jp2 = j+2
                  
#if defined(USEMPI)
                  if(j.eq.jdim)then
                    phim=(phi(i,j+1,k)-phi(i,j,k))*(pnorth(i,k)- 
     *                    phi(i,j+1,k))
                    if(phim.gt.0)then 
                       delc = 2.*phim/(pnorth(i,k)-phi(i,j,k)) 
                    endif
                 else
                    phim=(phi(i,j+1,k)-phi(i,j,k))*(phi(i,jp2,k)-
     *                   phi(i,j+1,k))
                    if(phim.gt.0)then
                       delc = 2.*phim/(phi(i,jp2,k)-phi(i,j,k))
                    endif
                 endif
#else 
                 if(j.eq.jdim)jp2 = 2
                 phim=(phi(i,j+1,k)-phi(i,j,k))*(phi(i,jp2,k)-
     *                phi(i,j+1,k))
                 if(phim.gt.0)then
                    delc = 2.*phim/(phi(i,jp2,k)-phi(i,j,k))
                 endif
#endif
                 fluxv1(i,k) = v(i,j+1,k)*(phi(i,j+1,k)-
     *                0.5*(1.+v(i,j+1,k)*dtdy)*delc - wcross - ucross)
              endif
c
c next the downstream flux
c
              fluxv0(i,k) = fluxvd(i,k)
              fluxvd(i,k) = fluxv1(i,k)
c
c compute v flux
c
              phip(i,j,k)=phi(i,j,k)-dtdy*vol(i,j,k)*
     *            (ay(i,j+1,k)*fluxv1(i,k)-ay(i,j,k)*fluxv0(i,k))
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
               wcross = 0.5*dtdz*wb*(phi(i,j,k)-phi(i,j,k-1))
            else
               if(k.eq.kdim) then
                  wcross = 0.0
               else
                  wcross = 0.5*dtdz*wb*(phi(i,j,k+1)-phi(i,j,k))
               endif
            endif
            vb = (v(i-1,j,k)+v(i,j,k)+v(i,j+1,k)+v(i-1,j+1,k))
     *           *0.25
            if(vb.gt.0) then
               vcross = 0.5*dtdy*vb*(phi(i,j,k)-phi(i,j-1,k))
            else
               vcross = 0.5*dtdy*vb*(phi(i,j+1,k)-phi(i,j,k))
            endif
c
            if(u(i,j,k).gt.0) then
               im2 = idim-1
               phim=(phi(i-1,j,k)-phi(im2,j,k))*(phi(i,j,k)-
     *              phi(i-1,j,k))
               if(phim.gt.0)then
                  delc = 2.*phim/(phi(i,j,k)-phi(im2,j,k))
               endif
               fluxub = u(i,j,k)*(phi(i-1,j,k)+0.5*(1.-u(i,j,k)*dtdx)*
     *              delc - wcross - vcross)
            else
               phim=(phi(i,j,k)-phi(i-1,j,k))*(phi(i+1,j,k)-
     *              phi(i,j,k))
               if(phim.gt.0)then
                  delc = 2.*phim/(phi(i+1,j,k)-phi(i-1,j,k))
               endif
               fluxub = u(i,j,k)*(phi(i,j,k)-0.5*(1.+u(i,j,k)*dtdx)*
     *              delc - wcross - vcross)
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
                  wcross = 0.5*dtdz*wb*(phi(i+1,j,k)-phi(i+1,j,k-1))
               else
                  if(k.eq.kdim) then
                     wcross = 0.0
                  else
                     wcross = 0.5*dtdz*wb*(phi(i+1,j,k+1)-phi(i+1,j,k))
                  endif
               endif
               vb = (v(i+1,j,k)+v(i,j,k)+v(i,j+1,k)+v(i+1,j+1,k))
     *              *0.25   
               if(vb.gt.0) then
                  vcross = 0.5*dtdy*vb*(phi(i+1,j,k)-phi(i+1,j-1,k))
               else
                  vcross = 0.5*dtdy*vb*(phi(i+1,j+1,k)-phi(i+1,j,k))
               endif
               if(u(i+1,j,k).gt.0) then
                  phim=(phi(i,j,k)-phi(i-1,j,k))*(phi(i+1,j,k)-
     *                 phi(i,j,k))
                  if(phim.gt.0)then
                     delc = 2.*phim/(phi(i+1,j,k)-phi(i-1,j,k))
                  endif
                  fluxu1 = u(i+1,j,k)*(phi(i,j,k)+
     *                 0.5*(1.-u(i+1,j,k)*dtdx)*delc - wcross - vcross)
               else
                  ip2 = i+2
                  if(i.eq.idim)ip2 = 2
                  phim=(phi(i+1,j,k)-phi(i,j,k))*(phi(ip2,j,k)-
     *                 phi(i+1,j,k))
                  if(phim.gt.0)then 
                     delc = 2.*phim/(phi(ip2,j,k)-phi(i,j,k)) 
                  endif 
                  fluxu1 = u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
               endif
c     
c     next the downstream flux
c     
               fluxu0 = fluxud
               fluxud = fluxu1
c
c compute u flux
c
               phip(i,j,k)=phip(i,j,k)-dtdx*vol(i,j,k)*
     *              (ax(i+1,j,k)*fluxu1-ax(i,j,k)*fluxu0)
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
               ucross = 0.5*dtdx*ub*(phi(i,j,k)-phi(i-1,j,k))
            else   
               ucross = 0.5*dtdx*ub*(phi(i+1,j,k)-phi(i,j,k))
            endif  
            vb = (v(i,j+1,k-1)+v(i,j,k-1)+v(i,j,k)+v(i,j+1,k))*0.25
            if(vb.gt.0) then
               vcross = 0.5*dtdy*vb*(phi(i,j,k)-phi(i,j-1,k))
            else  
               vcross = 0.5*dtdy*vb*(phi(i,j+1,k)-phi(i,j,k))
            endif
c
c
            if(w(i,j,k).gt.0) then
               phim=(phi(i,j,k)-phi(i,j,k-1))*(phi(i,j,k+1)-
     *              phi(i,j,k))
               if(phim.gt.0) then
                  delc = 2.*phim/(phi(i,j,k+1)-phi(i,j,k-1))
               endif
               fluxwb = w(i,j,k)*(phi(i,j,k-1) + 0.5*(1.-w(i,j,k)*dtdz)*
     *              delc - ucross - vcross)
            else
               phim=(phi(i,j,k)-phi(i,j,k-1))*(phi(i,j,k+1)-
     *              phi(i,j,k))
               if(phim.gt.0)then
                  delc = 2.*phim/(phi(i,j,k+1)-phi(i,j,k-1))
               endif
c
c do upstream
c
c     delc = 0.0
               fluxwb = w(i,j,k)*(phi(i,j,k)-0.5*(1.+w(i,j,k)*dtdz)*
     *              delc - ucross - vcross)
            endif
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
                  ucross = 0.5*dtdx*ub*(phi(i,j,k+1)-phi(i-1,j,k+1))
               else 
                  ucross = 0.5*dtdx*ub*(phi(i+1,j,k+1)-phi(i,j,k+1))
               endif
               vb = (v(i,j+1,k)+v(i,j,k)+v(i,j,k+1)+v(i,j+1,k+1))
     *              *0.25
               if(vb.gt.0) then
                  vcross = 0.5*dtdy*vb*(phi(i,j,k+1)-phi(i,j-1,k+1))
               else
                  vcross = 0.5*dtdy*vb*(phi(i,j+1,k+1)-phi(i,j,k+1))
               endif
               if(w(i,j,k+1).gt.0) then
                  phim=(phi(i,j,k)-phi(i,j,k-1))*(phi(i,j,k+1)-
     *                 phi(i,j,k))
                  if(phim.gt.0) then
                     delc = 2.*phim/(phi(i,j,k+1)-phi(i,j,k-1))
                  endif
                  fluxw1 = w(i,j,k+1)*(phi(i,j,k)+0.5*
     *                 (1.-w(i,j,k+1)*dtdz)*delc - ucross - vcross)
               else
                  kp2 = k+2
                  if(kp2.gt.kdim) kp2=kdim
                  phim=(phi(i,j,k+1)-phi(i,j,k))*(phi(i,j,kp2)-
     *                 phi(i,j,k+1))
                  if(phim.gt.0)then
                     delc = 2.*phim/(phi(i,j,kp2)-phi(i,j,k))
                  endif
                  fluxw1 = w(i,j,k+1)*(phi(i,j,k+1)-0.5*
     *                 (1.+w(i,j,k+1)*dtdz)*delc - ucross - vcross)
               endif
c     
c next the downstream flux
c
               fluxw0 = fluxwd
               fluxwd = fluxw1
c
               phip(i,j,k)=phip(i,j,k)-dtdz*vol(i,j,k)*(fluxw1-fluxw0)
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
      call MPI_BARRIER(MPI_COMM_WORLD,ier)
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
c      do k=1,kdim
c         do j=1,jdim
c            do i=1,idim
c               phip(i,j,k)=phip(i,j,k)+vol(i,j,k)*
c     *              ((ax(i+1,j,k)*u(i+1,j,k)-ax(i,j,k)*u(i,j,k))*dtdx+
c     *              (w(i,j,k+1)-w(i,j,k))*dtdz
c     *              + (ay(i,j+1,k)*v(i,j+1,k)-ay(i,j,k)*v(i,j,k))*dtdy)
c     *              *phi(i,j,k)
c            enddo
c         enddo
c      enddo
c
c finish up the processing by mixing through turbulence parameterization.
c Don't do this from the TKE fields.  TKE is diffused in the sgs update
c routine.
c
      if(ivar.ne.MDTKE)then
         do k=1,kdim
            do j=1,jdim
               do i=1,idim
                  kright = kheat(i+1,j,k)+kheat(i,j,k)
                  kleft = kheat(i,j,k)+kheat(i-1,j,k)
                  b1=(ax(i+1,j,k)*kright*(phi(i+1,j,k)-phi(i,j,k))-
     *                 ax(i,j,k)*kleft*(phi(i,j,k)-phi(i-1,j,k)))*
     *                 0.5*odx*odx
c
                  kback = kheat(i,j+1,k)+kheat(i,j,k)
                  kfront = kheat(i,j,k)+kheat(i,j-1,k)
                  b2=(ay(i,j+1,k)*kback*(phi(i,j+1,k)-phi(i,j,k))-
     *                 ay(i,j,k)*kfront*(phi(i,j,k)-phi(i,j-1,k)))*
     *                 0.5*ody*ody
c
                  ktop = kheat(i,j,k+1)+kheat(i,j,k)
                  if(k.ge.kdist(i,j)) ktop = 0.0
                  kbot = kheat(i,j,k)+kheat(i,j,k-1)
#if defined(RADIATE)
                  if(k.eq.1)then
                     kbot = 0.0
                  endif
#endif
                  b3=(ktop*(phi(i,j,k+1)-phi(i,j,k))-
     *                 kbot*(phi(i,j,k)-phi(i,j,k-1)))*0.5*odz*odz
                  sgs=vol(i,j,k)*(b1+b2)+b3
                  phip(i,j,k)=phip(i,j,k)+delt*sgs
               enddo
            enddo
         enddo
      endif
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
