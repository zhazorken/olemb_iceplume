      subroutine icinit(ier)
c
c   $Id: icinit.F,v 1.7 1998/08/04 18:12:48 dwd Exp $
c
c     This routine parses the I.C. input cards and adds perturbations
c     to the temperature, vertical velocity, and tke fields.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "olemi.inc" 
#include "olemtke.inc"
#include "cyclic.inc"
#include "communicate.inc"

      character value*80
      real s,xtemp

      real ran2,ffsave,sigthta,rhos,qs,ws,ustar,ts
      real drhodtheta,tau
      real zi1,zz,zt,zw,iestbl
#if defined(USEMPI)
c      real phi1(0:idim1,0:jtotal1,0:kdim1)
#endif
      real*8 drho,drhol,filtave
      integer i,j,k
      integer unts(3),fosave,ier,ierr

      logical flgs(3),mkstst,fvsave
c
c parse I.C. flags  icfstr (icuseran)
c
      ier=0
c
      if(isrestart)then
         if(iope) then
           write(lunit,4)
 4         format(/'Restart -- IC flags and values ignored')
         endif
      else
         if(icfstr.eq.' ')then
            icflgs(1)=0
         else
            call parse(icfstr,value)
            call mks(value,xtemp,flgs,unts)
            if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
            icflgs(1)=int(xtemp)
         endif
         icuseran=icflgs(1).eq.1
c     
c parse I.C. values  icvstr (iczi, icramp, icseed)
c
         if(icvstr.eq.' ')then
            icvals(1)=100.0
         else
            call parse(icvstr,value)
            call mks(value,icvals(1),flgs,unts)
            if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 2000
         endif
         iczi=icvals(1)
c
         if(icvstr.eq.' ')then
            icvals(2)=0.1
         else
            call parse(icvstr,value)
            call mks(value,icvals(2),flgs,unts)
            if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         endif
         icramp=icvals(2)
c     
         if(icvstr.eq.' ')then
            icvals(3)=-444765
         else
            call parse(icvstr,value)
            call mks(value,icvals(3),flgs,unts)
            if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         endif
         icrseed=int(icvals(3))
#if defined(USEMPI)
         icrseed=int(icvals(3))*float(my_pe+1)/float(nprocs)
#endif
         write(*,*)"Random seed ",icrseed,my_pe
c
c write out IC flags and values
c
         if(iope) then
           write(lunit,6)
 6         format(/'I.C. flags:')
           if(icuseran)then
              write(lunit,3)1,'Use random perturbations'
 3            format(4x,i2,'. ',a)
           else
              write(lunit,3)1,'Do not use random perturbations'
           endif
c
           write(lunit,1) iczi,icramp,icrseed
 1         format(/'I.C. values:'/
     *          5x,'1. mixed layer depth = ',f12.3,' (m)'/
     *          5x,'2.  random amplitude = ',1pg12.3/
     *          5x,'3.       random seed = ',i12)
        endif
c
         if(icuseran)then
c
c create random field, use tkep to save space
c
c#if defined(USEMPI)
c            if(iope)then
c              do j=0,jtotal1
c                 do k=0,kdim1
c                    do i=0,idim1
c                       phi1(i,j,k)=(0.5-ran2(icrseed))*icramp
c                    enddo
c                 enddo
c              enddo
c            endif
c            call distrib_3d(tkep,phi1)
c#else
            do j=0,jdim1
              do k=0,kdim1
                do i=0,idim1
                  tkep(i,j,k)=(0.5-ran2(icrseed))*icramp
                enddo
              enddo
            enddo
c#endif
c
c filter field, use thp to save space
c
            fosave=filtord
            ffsave=filtfac
            fvsave=filtvert
c
            filtord=2
            filtfac=0.5
            filtvert=.true.
            call filter(tkep,thp)
            do j=1,jdim
               do k=1,kdim
                  do i=1,idim
                     thp(i,j,k)=tkep(i,j,k)+thp(i,j,k)
                  enddo
               enddo
            enddo
c
            filtord=fosave
            filtfac=ffsave
            filtvert=fvsave
c
c compute filtered average perturbation and demean filterd field
c
            do k=1,kdim
               filtave=0.0
               do j=1,jdim
                  do i=1,idim
                     filtave=filtave+thp(i,j,k)
                  enddo
               enddo
               filtave=filtave/(jdim*idim)
               do j=1,jdim
                  do i=1,idim
                     thp(i,j,k)=thp(i,j,k)-filtave
                  enddo
               enddo
            enddo
c
c compute t_star, w_star, u_star, and q
c
             qs=-heatf/(cp*rho0)
c
c compute the volumetric expansion at the surface
c
             s=sal(1,1,kdim)
             drhodtheta=-(sigthta(s,th(1,1,kdim)+0.001,0.0) - 
     *            sigthta(s,th(1,1,kdim),0.0))/0.001
c
             rhos=drhodtheta*qs
c
             if(rhos.le.0)then
                ws=0.0
                ts=0.0
             else
                ws=(grav*iczi*rhos/rho0)**0.3333333
                ts=qs/ws
             endif
             tau = sqrt(ustr**2 + vstr**2)
             ustar = sqrt(tau/rho0)
c
             if(iope) then
               write(lunit,5)ws,ustar,ts,qs
 5          format(/5x,'              w_star = ',1pe12.3,' (m/s)'/
     *              5x,'              u_star = ',e12.3,' (m/s)'/
     *              5x,'              t_star = ',e12.3,' (degC)'/
     *              5x,'                  Qs = ',e12.3,' (degC m/s)'/)
            endif
c     
             zi1=0.844*iczi
             zz=1.2*iczi
c  
c---------Initialize variables, if data does not already exist (u,v,w,tke).
c
             do j=1,jdim
                do 145 i=1,idim
                   do 140 k=0,kdim1
                      zt=(kdim-k+0.5)*dz
                      zw=(kdim1-k)*dz
c
                      if(icuseran.and.(k.ne.kdim1))then
                         if(zt.le.zi1)then
                            th(i,j,k)=th(i,j,k)+thp(i,j,k)*(1-zt/zi1)*ts
                         endif
                         if(zw.le.zi1)then
                            w(i,j,k)=w(i,j,k)+thp(i,j,k)*(1-zw/zi1)*ws
                         endif
                         if(zt.le.zz)then
                            tke(i,j,k)=tke(i,j,k)+
     *                           sqrt(icramp*ws*ws*(1-zt/zz))
                         endif
                      endif
c     
                      um(i,j,k)=u(i,j,k)
                      umm(i,j,k)=u(i,j,k)
                      vm(i,j,k)=v(i,j,k)
                      vmm(i,j,k)=v(i,j,k)
                      if(k.ne.0) then
                         wm(i,j,k)=w(i,j,k)
                         wmm(i,j,k)=w(i,j,k)
                      endif
 140               continue
 145            continue
             enddo
             call imag(u)
             call imag(um)
             call imag(umm)
             call imag(v)
             call imag(vm)
             call imag(vmm)
             call imag(w)
             call imag(wm)
             call imag(wmm)
             call imag(tke)
c
c recompute rhobm1, rho0, rhop, and pressure
c
             rho0=0.0
             do 105 k=1,kdim
                drhol=0.0
                do 106 j=1,jdim
                   do 107 i=1,idim
                      s=sal(i,j,k)
                      drhol=drhol+iestbl(s,th(i,j,k),2*k)
 107               continue
 106            continue
                call MPI_ALLREDUCE(drhol,drho,1,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,comm, ierr)

                rhobm1(k)=drho/(jtotal*idim)
                rhob(k)=1000.0+rhobm1(k)
                rho0=rho0+rhob(k)
 105         continue
             rhob(0)=rhob(1)
             rho0=rho0/float(kdim)
             do j=1,jdim
                do 161 k=1,kdim
                   do 162 i=1,idim
                      s=sal(i,j,k)
                      rhop(i,j,k)=iestbl(s,th(i,j,k),2*k)-rhobm1(k)
 162               continue
 161            continue
                
                do 163 i=1,idim
                   rhop(i,j,0)=rhop(i,j,1)
                   rhop(i,j,kdim1)=rhop(i,j,kdim)
 163            continue
                
                do 180 i=1,idim
                   p(i,j,kdim)=0.0
                   do 170 k=kdim,2,-1
                      p(i,j,k-1)=p(i,j,k)+
     &                     0.5*grav*(rhop(i,j,k)+rhop(i,j,k-1))*dz
 170               continue
 180            continue
             enddo
#if defined(SNGLSALT)
             call imag(sal)
#else 
             call imagdp(sal)
#endif
             call imag(th)
             call imag(rhop)
             call imag(p)
          endif
          
       endif
       return
 2000  write(*,999)value(:30)
 999   format('Units error on I.C. card for entry: ',a)
       ier=1
       return
       end






