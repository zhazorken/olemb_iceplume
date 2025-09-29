      subroutine leer(n,phi,phip,ivar)
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
#include "communicate.inc"
#include "finitv.inc"
c
      integer i,j,k,ip2,jp2,im2,jm2,kp2
      real phip(0:idim1,0:jdim1,0:kdim1),phi(0:idim1,0:jdim1,0:kdim1)
#if defined(USEMPI)
      integer ier
      real pnorth(0:idim1,0:kdim1)
      real psouth(0:idim1,0:kdim1)
      real peast(0:jdim1,0:kdim1)
      real pwest(0:jdim1,0:kdim1)
#endif
      real phim,delc,fluxu1,fluxu0,fluxw1,fluxw0,fluxvb
      real fluxwd,fluxud,fluxub,fluxwb,ub,vb,wb,ucross,wcross,vcross
      real sgs,dtsub
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
#if defined(DISCHARGE)
      do k=1,kdim
        do j=0,jdim1
          do i=1,idim

           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
           jj = jblkstart(my_pe+1)+j

           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
c               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               th(i-1,j,k) = 0.5
               thp(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
               salp(i-1,j,k) = 1.0
               tke(i-1,j,k) = sedinit
               tkep(i-1,j,k) = sedinit
             endif
           endif
          endif
          enddo
        enddo
      enddo
#endif
      call imag(sal)
      call imag(th)

#if defined(USEMPI) 
      call northshift(phi,pnorth) 
      call MPI_BARRIER(comm,ier)
      call southshift(phi,psouth)
      call MPI_BARRIER(comm,ier)
      call westshift(phi,pwest)
      call MPI_BARRIER(comm,ier)
      call eastshift(phi,peast)
      call MPI_BARRIER(comm,ier)
#endif

      j=1
      do k=1,kdim
         do i=1,idim
c
c
c
            dtdx = delt/dx(i)
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
            wcross = 0.0
            ub = (u(i,j-1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j-1,k))
c     *           ax(i,j-1,k)*ax(i,j,k)*ax(i+1,j,k)*ax(i+1,j-1,k)*0.25
            if(ub.gt.0) then
               ucross = 0.5*dtdx*ub*(phi(i,j,k)-phi(i-1,j,k))
            else
               ucross = 0.5*dtdx*ub*(phi(i+1,j,k)-phi(i,j,k))
            endif
            ucross = 0.0
c
            if(iblkstart(my_pe+1)+i-1.eq.1) ucross = 0.0
c
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
#if defined(CHANNEL)
            if(jblkstart(my_pe+1).eq.1) fluxvb = 0.0
#endif
            
            fluxvd(i,k) = fluxvb
         enddo
      enddo
c
      do j=1,jdim
         do k=1,kdim
            do i=1,idim
c
c make sure discharge is in force
c
c            if(ivar.eq.MDTKE.and.k.ne.1) then
c               w(i,j,k) = w(i,j,k)-0.04*way(i,j,k)
c            endif
c               disch = 600.
c               jj = jblkstart(my_pe+1)+j
c               if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c               if(jj.gt.jtotal/2-30.and.jj.lt.jtotal/2+30) then
c                 if(k.lt.5) then
c                   u(i,j,k) = disch/(4*60.*dxx*dxx)
c                   um(i,j,k) = disch/(4*60.*dxx*dxx)
c                   umm(i,j,k) = disch/(4*60.*dxx*dxx)
c                   ax(i,j,k) = 1.0
c                   th(i-1,j,k) = 0.5
c                   sal(i-1,j,k) = 1.0
c    
c                 endif
c               endif
c               endif


c
c do horizontal terms first
c
c
c first, the upstream flux
c
c 
               dtdx = delt/dx(i)
c

               delc = 0.0
               wb = (w(i,j+1,k)+w(i,j,k)+w(i,j,k+1)+w(i,j+1,k+1))
     *              *0.25
               if(wb.gt.0) then
                  wcross = 0.5*dtdz*wb*(phi(i,j+1,k)-phi(i,j+1,k-1))
                  wcross = wcross*vay(i,j+1,k)*vay(i,j+1,k-1)
               else
                  if(k.eq.kdim) then
                     wcross = 0.0
                  else
                     wcross = 0.5*dtdz*wb*(phi(i,j+1,k+1)-phi(i,j+1,k))
                  endif
                  wcross = wcross*vay(i,j+1,k)*vay(i,j+1,k-1)
               endif
               ub = (u(i,j+1,k)+u(i,j,k)+u(i+1,j,k)+u(i+1,j+1,k))
c     *           *0.25*ax(i,j+1,k)*ax(i,j,k)*ax(i+1,j,k)*ax(i+1,j+1,k)
               if(ub.gt.0) then
                  ucross = 0.5*dtdx*ub*(phi(i,j+1,k)-phi(i-1,j+1,k))
               else
                  ucross = 0.5*dtdx*ub*(phi(i+1,j+1,k)-phi(i,j+1,k))
               endif
               ucross = ucross*bax(i,j+1,k)*bax(i-1,j+1,k)*
     *                  bax(i+1,j+1,k)

c
               if(iblkstart(my_pe+1)+i-1.eq.1) then
                   ucross = 0.0
                   wcross = 0.0
               endif
c
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
#if defined(CHANNEL)
c            if(jblkstart(my_pe+1)+jdim.ge.jtotal) fluxv1(i,k) = 0.0
#endif
c
c next the downstream flux
c
              fluxv0(i,k) = fluxvd(i,k)
              fluxvd(i,k) = fluxv1(i,k)
c
#if defined(CHANNEL)
c                if(jblkstart(my_pe+1)+j.ge.jtotal) then
c                   fluxv1(i,k) = 0.0
c                else if(jblkstart(my_pe+1)+j.eq.2) then
c                   fluxv0(i,k) = 0.0
c                endif
#endif
c
c compute v flux
c
              phip(i,j,k)=phi(i,j,k)-dtdy*(fluxv1(i,k)-fluxv0(i,k))
c     *           *vay(i,j,k)
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
c 
            dtdx = delt/dx(i)
c
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
#if defined(CHANNEL)
               if(jblkstart(my_pe+1)+j.ge.jtotal) then
                  vcross = 0.0
               endif
#endif
            endif
c
            if(u(i,j,k).gt.0) then
               im2 = idim-1
#if defined(USEMPI)
               phim=(phi(i-1,j,k)-pwest(j,k))*(phi(i,j,k)-
#else
               phim=(phi(i-1,j,k)-phi(im2,j,k))*(phi(i,j,k)-
#endif
     *              phi(i-1,j,k))
               if(phim.gt.0)then
#if defined(USEMPI)
                  delc = 2.*phim/(phi(i,j,k)-pwest(j,k))
#else
                  delc = 2.*phim/(phi(i,j,k)-phi(im2,j,k))
#endif
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
c#if defined(DISCHARGE)

c           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
c           jj = jblkstart(my_pe+1)+j
c
c           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
c             if(k.lt.dphip+1) then
c               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               fluxub = u(i,j,k)*phi(i-1,j,k)
c             endif
c           endif
c          endif
c#endif


c            if(iblkstart(my_pe+1).eq.1.and.k.gt.2) then
c              fluxub = 0.0
c            endif

c            fluxud = ax(i,j,k)*fluxub
            fluxud = fluxub

            do i=1,idim 
c
c do horizontal terms first
c
c
c first, the upstream flux
c
c 
               dtdx = delt/dx(i)
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
               wcross = wcross*vay(i+1,j,k)*vay(i+1,j,k-1)

               vb = (v(i+1,j,k)+v(i,j,k)+v(i,j+1,k)+v(i+1,j+1,k))
     *              *0.25   
               if(vb.gt.0) then
                  vcross = 0.5*dtdy*vb*(phi(i+1,j,k)-phi(i+1,j-1,k))
               else
                  vcross = 0.5*dtdy*vb*(phi(i+1,j+1,k)-phi(i+1,j,k))
#if defined(CHANNEL)
                  if(jblkstart(my_pe+1)+j.ge.jtotal) then
                    vcross = 0.0
                  endif
#endif
               endif
               if(u(i+1,j,k).gt.0) then
                  phim=(phi(i,j,k)-phi(i-1,j,k))*(phi(i+1,j,k)-
     *                 phi(i,j,k))
                  if(phim.gt.0)then
                     delc = 2.*phim/(phi(i+1,j,k)-phi(i-1,j,k))
                  endif
                  if(bax(i+1,j,k).eq.0.0) then
                    delc = 0.0
                    wcross = 0.0
                    vcross = 0.0
                  endif
                  fluxu1 = u(i+1,j,k)*(phi(i,j,k)+
     *                 0.5*(1.-u(i+1,j,k)*dtdx)*delc - wcross - vcross)
               else
                  ip2 = i+2
#if defined(USEMPI)
                  if(i.eq.idim)then
                   phim=(phi(i+1,j,k)-phi(i,j,k))*(peast(j,k)-
     *                 phi(i+1,j,k))
                    if(phim.gt.0)then
                       delc = 2.*phim/(peast(j,k)-phi(i,j,k))
                    endif
                    fluxu1 = u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
                  else
                    phim=(phi(i+1,j,k)-phi(i,j,k))*(phi(ip2,j,k)-
     *                 phi(i+1,j,k))
                    if(phim.gt.0)then
                       delc = 2.*phim/(phi(ip2,j,k)-phi(i,j,k))
                    endif
                    fluxu1 = u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
                  endif

#else
                  if(i.eq.idim)ip2 = 2
                  phim=(phi(i+1,j,k)-phi(i,j,k))*(phi(ip2,j,k)-
     *                 phi(i+1,j,k))
                  if(phim.gt.0)then 
                     delc = 2.*phim/(phi(ip2,j,k)-phi(i,j,k)) 
                  endif 
                  fluxu1 = u(i+1,j,k)*(phi(i+1,j,k)-
     *                 0.5*(1.+u(i+1,j,k)*dtdx)*delc - wcross - vcross)
#endif
               endif
c #if defined(DISCHARGE)

c           if(bax(i+1,j,k).eq.0.0.and.bax(i+2,j,k).ne.0.0) then
c           jj = jblkstart(my_pe+1)+j
c
c           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
c             if(k.lt.dphip+1) then
c               u(i+1,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               fluxu1 = u(i+1,j,k)*phi(i,j,k)
c             endif
c           endif
c          endif
c #endif

c     
c     next the downstream flux
c     
               fluxu0 = fluxud
               fluxud = fluxu1
c
c compute u flux
c
               phip(i,j,k)=phip(i,j,k)-dtdx*(ax(i+1,j,k)*fluxu1-
     *             ax(i,j,k)*fluxu0)
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
c 
            dtdx = delt/dx(i)
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
            ucross = 0.0
            vcross = 0.0
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
            if(ivar.eq.MDTEMP) then
              wt_ave(k) = wt_ave(k)+fluxwd
            else if(ivar.eq.MDSAL) then
              ws_ave(k) = ws_ave(k)+fluxwd
            endif

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
                  ucross = ucross*bax(i,j,k+1)*bax(i-1,j,k+1)
               else 
                  ucross = 0.5*dtdx*ub*(phi(i+1,j,k+1)-phi(i,j,k+1))
                  ucross = ucross*bax(i+1,j,k+1)*bax(i,j,k+1)
               endif
c
               if(iblkstart(my_pe+1)+i-1.eq.1) ucross = 0.0
c
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
                  if(kp2.gt.kdim) kp2=kdim1
                  phim=(phi(i,j,k+1)-phi(i,j,k))*(phi(i,j,kp2)-
     *                 phi(i,j,k+1))
                  if(phim.gt.0)then
                     delc = 2.*phim/(phi(i,j,kp2)-phi(i,j,k))
                  endif
                  fluxw1 = w(i,j,k+1)*(phi(i,j,k+1)-0.5*
     *                 (1.+w(i,j,k+1)*dtdz)*delc - ucross - vcross)
               endif
               if(k.eq.kdim) fluxw1 = 0.0
c     
c next the downstream flux
c
               fluxw0 = fluxwd
               fluxwd = fluxw1
               if(ivar.eq.MDTEMP) then
                 wt_ave(k+1) = wt_ave(k+1)+fluxw1
               else if(ivar.eq.MDSAL) then
                 ws_ave(k+1) = ws_ave(k+1)+fluxw1
c               else if(ivar.eq.MDQC) then
c                 wqc_ave(k+1) = wqc_ave(k+1)+fluxw1
               endif
c
               phip(i,j,k)=phip(i,j,k)-dtdz*(fluxw1-
     *              fluxw0)
c     *           *vay(i,j,k)
c     
            enddo
         enddo
      enddo
      do k=1,kdim
        do j=0,jdim1
          do i=1,idim
#if defined(DISCHARGE)

           if(bax(i,j,k).eq.0.0.and.bax(i+1,j,k).ne.0.0) then
           jj = jblkstart(my_pe+1)+j

           if(jj.gt.jtotal/2-plumw.and.jj.lt.jtotal/2+plumw) then
             if(k.lt.dphip+1) then
c               u(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               umm(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
c               um(i,j,k) = disch/(dphip*2*plumw*dxx*dxx)
               th(i-1,j,k) = 0.5
               thp(i-1,j,k) = 0.5
               sal(i-1,j,k) = 1.0
               salp(i-1,j,k) = 1.0
               tke(i-1,j,k) = sedinit
               tkep(i-1,j,k) = sedinit
             endif
           endif
          endif
#endif
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
c 
               dtdx = delt/dx(i)
c
               phip(i,j,k)=phip(i,j,k)+((u(i+1,j,k)-
     *               u(i,j,k))*dtdx+
     *              (w(i,j,k+1)-w(i,j,k))*dtdz
     *              + (v(i,j+1,k)-v(i,j,k))*dtdy)
     *                *phi(i,j,k)
c     *           *vay(i,j,k)
cc     *              (rhob(k+1)*w(i,j,k+1)-rhob(k)*w(i,j,k))*dtdz*2.0/
cc     *              (rhob(k+1)+rhob(k))
cc     *              + (v(i,j+1,k)-v(i,j,k))*dtdy-
cc     *              (wm(i,j,k)+wm(i,j,k+1))*(rhob(k+1)-rhob(k))*dtdz/
cc     *              (rhob(k+1)+rhob(k)))*phi(i,j,k)
c
            enddo
         enddo
      enddo
c
c finish up the processing by mixing through turbulence parameterization.
c Don't do this from the TKE fields.  TKE is diffused in the sgs update
c routine.
c
c      do j=1,jdim
c        do i=1,idim
c          do k=2,kdim-30
c              dtsub = 0.05*sin(3.141*float(k)/float(kdim-30))
c              phip(i,j,k) = phip(i,j,k)-0.5*dtsub*
c     *          (phi(i,j,k+1)-phi(i,j,k-1))*dtdz
c          enddo
cc          phip(i,j,kdim) = phip(i,j,kdim)+(phip(i,j,kdim-1)-
cc     *        phi(i,j,kdim-1))
c        enddo
c      enddo

#if defined(SUBLIM)
      if(mod(n,120).eq.0) then
      do j=1,jdim
        do i=1,idim
          do k=2,kdim-1
c            dtsub = 1.75e-6*(k*dz-0.5*dz)
            dtsub = 5.e-6*(k*dz-0.5*dz)
            phip(i,j,k) = phip(i,j,k)+0.5*dtsub*120.*
     *          (phi(i,j,k+1)-phi(i,j,k-1))*dtdz
c     *           *vay(i,j,k)

          enddo
          phip(i,j,kdim) = phip(i,j,kdim)+(phip(i,j,kdim-1)-
     *        phi(i,j,kdim-1))
     *           *vay(i,j,k)
        enddo
      enddo
      endif
#endif

c      if(ivar.ne.MDTKE)then
         do k=1,kdim
            do j=1,jdim
               do i=1,idim
c                  if(ivar.eq.MDTKE.and.k.ne.1) then
c                        w(i,j,k)=w(i,j,k)+0.04
c                  endif
c
                  odx = 1./dx(i)
c
                  kright = kheat(i+1,j,k)+kheat(i,j,k)
                  kleft = kheat(i,j,k)+kheat(i-1,j,k)
#if defined(HTOP)
                  kright = kright*vay(i,j,k)*vay(i+1,j,k)
                  kleft = kleft*vay(i,j,k)*vay(i-1,j,k)
#endif
                  if(iblkstart(my_pe+1)+i-1.eq.1) then
                    kleft = 0.0
                  endif
                  b1=(kright*(phi(i+1,j,k)-phi(i,j,k))-
     *                 kleft*(phi(i,j,k)-phi(i-1,j,k)))*0.5*odx*odx
c
                  kback = kheat(i,j+1,k)+kheat(i,j,k)
                  kfront = kheat(i,j,k)+kheat(i,j-1,k)
#if defined(HTOP)
                  kback = kback*ay(i,j+1,k)
                  kfront = kfront*ay(i,j-1,k)
#endif
                  b2=(kback*(phi(i,j+1,k)-phi(i,j,k))-
     *                 kfront*(phi(i,j,k)-phi(i,j-1,k)))*0.5*ody*ody
c
                  ktop = kheat(i,j,k+1)+kheat(i,j,k)
                  kbot = kheat(i,j,k)+kheat(i,j,k-1)
#if defined(HTOP)
                  ktop = ktop*bway(i,j,k+1)
                  kbot = kbot*bway(i,j,k)
#endif
#if defined(CHANNEL)
                  if(jblkstart(my_pe+1)+j.ge.jtotal) then
                    b2 = 0.0
                  else if(jblkstart(my_pe+1)+j.eq.2) then
                    b2 = 0.0
                  endif
#endif
c
#if defined(RADIATE)
                  if(k.eq.1)then
                     kbot = 0.0
                  endif
#endif
                  b3=(ktop*(phi(i,j,k+1)-phi(i,j,k))-
     *                 kbot*(phi(i,j,k)-phi(i,j,k-1)))*0.5*odz*odz
                  sgs=b1+b2+b3
                  phip(i,j,k)=phip(i,j,k)+delt*sgs
     *           *vay(i,j,k)
                  if(ivar.eq.MDTEMP) then
                    wts_ave(k+1) = wts_ave(k+1) -
     *                        ktop*(phi(i,j,k+1)-phi(i,j,k))*0.5*odz
                  else if(ivar.eq.MDSAL) then
                    wss_ave(k+1) = wss_ave(k+1) -
     *                        ktop*(phi(i,j,k+1)-phi(i,j,k))*0.5*odz
                  endif
                  if(ivar.eq.MDTKE) then
                    if(bax(i,j,k).ne.0.0) then
                     phip(i,j,k) = phip(i,j,k)+delt*
     *               0.04*(phi(i,j,k+1)-
     *                     phi(i,j,k))*odz
c     *           *vay(i,j,k)
                     if(k.eq.1) then
                       phip(i,j,k) = phip(i,j,k)-0.04*delt*
     *                   phi(i,j,k)*odz
c     *           *vay(i,j,k)
                     endif
c                     if(phip(i,j,k).lt.0.0) then
c                       phip(i,j,k) = 0.0
c                     endif
                    endif
                  endif
               enddo
            enddo
         enddo
#if defined(APLUME)
      if(ivar.eq.MDTEMP) then
         do i=0,idim1
           do j=0,jdim1
             jj = j+jblkstart(my_pe+1)
             if(iblkstart(my_pe+1)+i.gt.507.and.
     *          iblkstart(my_pe+1)+i.lt.518) then
                if(jj.gt.251.and.jj.lt.262) then
                 phip(i,j,1)=phip(i,j,1)-0.1*delt*(phi(i,j,1)-0.0)*odz
                endif
             endif
           enddo
         enddo
      else if(ivar.eq.MDSAL) then
         do i=0,idim1
           do j=0,jdim1
             jj = j+jblkstart(my_pe+1)
             if(iblkstart(my_pe+1)+i.gt.507.and.
     *          iblkstart(my_pe+1)+i.lt.517) then
                if(jj.gt.251.and.jj.lt.261) then
                 phip(i,j,1)=phip(i,j,1)-0.1*delt*(phi(i,j,1)-2.0)*odz
                endif
             endif
           enddo
         enddo
      endif
#endif


c      endif
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
