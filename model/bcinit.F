      subroutine bcinit (nstep,com1,com2,com3,ier)
c 
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "cyclic.inc"
#include "olemtke.inc"
#include "communicate.inc"
c 
      integer i,j,unts(3),ier,nstep,lnblk,len
      logical flgs(3),mkstst
      real ran2,xtemp
      real phi2(itotal,jtotal)
      real randt(idim,jdim)
      character*120 com1,com2,com3,value*80
#if defined(RADIATE) || defined(RADTOP)
      integer js
      real rki(idim),rkj(jtotal),rpi,drhodz
      data rpi/3.1415927/
#endif
c
c 
      ier=0
      com1='$Id: bcinit.F,v 1.6 1998/10/21 21:56:34 dwd Exp $'
      com2='Cyclic Condition B.C. at sides w=0 z=0,-H'
c
c parse B.C. flags   bcfstr
c          (uselatent,userandom,useswheat,usestress,userain)
c
      if(bcfstr.eq.' ')then
         bcflgs(1)=0
      else
         call parse(bcfstr,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         bcflgs(1)=int(xtemp)
      endif
      uselatent=bcflgs(1).eq.1
c
      if(bcfstr.eq.' ')then
         bcflgs(2)=1
      else
         call parse(bcfstr,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         bcflgs(2)=int(xtemp)
      endif
      bcuseran=bcflgs(2).eq.1
      if(bcuseran)then
         com3='Random heating at top  du/dz = 0, z=0'
      else
         com3='Uniform heating at top  du/dz = 0, z=0'
      endif
c
      if(bcfstr.eq.' ')then
         bcflgs(3)=1
      else
         call parse(bcfstr,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         bcflgs(3)=int(xtemp)
      endif
      useswheat=bcflgs(3).eq.1
c
      if(bcfstr.eq.' ')then
         bcflgs(4)=1
      else
         call parse(bcfstr,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         bcflgs(4)=int(xtemp)
      endif
      usestress=bcflgs(4).eq.1

      if(bcfstr.eq.' ')then
         bcflgs(5)=1
      else
         call parse(bcfstr,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
         bcflgs(5)=int(xtemp)
      endif
      userain=bcflgs(5).eq.1
c
c parse B.C. values  bcvstr (stop, ranseed, ranamp, ranoff, cool, latent)
c
      call parse(bcvstr,value)
      call mks(value,bcvals(1),flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 2000
      ncool=int(bcvals(1)/delt+0.5)
c
      call parse(bcvstr,value)
      call mks(value,bcvals(2),flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
      nf=int(bcvals(2))
c
      call parse(bcvstr,value)
      call mks(value,bcvals(3),flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
      ranamp=bcvals(3)
c
      call parse(bcvstr,value)
      call mks(value,bcvals(4),flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 2000
      nrndom=int(bcvals(4)/delt+0.5)
c
c Parse B.C. forcing bcfrcstr
c         (cool, latent, ustokes, vstokes, ustress, vstress)
c
c
c W m-2 == kg s-3   (W == kg m2 s-3)
c
      call parse(bcfrcstr,value)
      call mks(value,heatf,flgs,unts)
      if(mkstst(flgs,unts,.false.,.true.,.true.,0,1,-3))goto 2000
      bcfrce(1)=heatf
c
      call parse(bcfrcstr,value)
      call mks(value,latent,flgs,unts)
      if(mkstst(flgs,unts,.false.,.true.,.true.,0,1,-3))goto 2000
      if(.not.uselatent) then
         lheaton=.false.
         latent=0.0
      endif
      bcfrce(2)=latent
c
c N m-2 == kg m-1 s-2  (N == kg m s-2)
c
      call parse(bcfrcstr,value)
      call mks(value,ustk,flgs,unts)
      if(mkstst(flgs,unts,.true.,.false.,.true.,1,0,-1))goto 2000
      bcfrce(3)=ustk
c
      call parse(bcfrcstr,value)
      call mks(value,vstk,flgs,unts)
      if(mkstst(flgs,unts,.true.,.false.,.true.,1,0,-1))goto 2000
      bcfrce(4)=vstk
c
      call parse(bcfrcstr,value)
      call mks(value,ustr,flgs,unts)
      if(mkstst(flgs,unts,.true.,.true.,.true.,-1,1,-2))goto 2000
c
      call parse(bcfrcstr,value)
      call mks(value,vstr,flgs,unts)
      if(mkstst(flgs,unts,.true.,.true.,.true.,-1,1,-2))goto 2000
      if(.not.usestress) then
         uston=.false.
         vston=.false.
         ustr=0.0
         vstr=0.0
      endif
      bcfrce(5)=ustr
      bcfrce(6)=vstr
c
c interpret B.C. Flags, Values, and Forces:
c
      if(iope) then
c
        len=lnblk(com1,120)
        write(lunit,1) com1(:len)
 1      format(/'B.C. Comments:'/1x,a)
        len=lnblk(com2,120)
        write(lunit,6)com2(:len)
 6      format(1x,a)
        len=lnblk(com3,120)
        write(lunit,6)com3(:len)
c
        write(lunit,5)
 5      format(/'B.C. flags:')
        if(uselatent)then
           write(lunit,3)1,'Use latent heat flux'
 3         format(4x,i2,'. ',a)
        else
           write(lunit,3)1,'Do not use latent heat flux'
        endif
        if(bcuseran)then
           write(lunit,3)2,'Use random surface cooling'
        else
           write(lunit,3)2,'Use uniform surface cooling'
        endif
        if(useswheat) then
           write(lunit,3)3,'Use short wave heat flux'
           swheaton=.true.
        else
           write(lunit,3)3,'Do not use short wave heat flux'
        endif
        if(usestress) then
           write(lunit,3)4,'Use surface wind stress'
        else
           write(lunit,3)4,'Do not use surface wind stress'
        endif
        if(userain) then
           write(lunit,3)5,'Use rainfall forcing'
        else
           write(lunit,3)5,'Do not use rainfall forcing'
        endif
c     
c set defaults
c
        if(.not.useswheat) then
           swheaton=.false.
           swheatf=0.0
        endif
        if(.not.userain) then
           rainon=.false.
           rainf=0.0
        endif
c
        write(lunit,2) bcvals(1)/3600.,nf,ranamp,bcvals(4)/3600.
2       format(/'B.C. values:'/
     *     5x,'1.          turn cooling off at = ',f12.3,' (hours)'/
     *     5x,'2.                  random seed = ',i12/
     *     5x,'3. amplitude of random function = ',1pg12.3/
     *     5x,'4.   turn random forcing off at = ',0pf12.3,' (hours)')
c
        write(lunit,4) heatf, latent, ustk,vstk,ustr,vstr
 4      format(/'B.C. default forcing:'/
     *     5x,'1. maximum heat flux = ',1pg12.3,' (W m-2)'/
     *     5x,'2.  latent heat flux = ',g12.3,' (W m-2)'/
     *     5x,'3. u Stokes velocity = ',g12.3,' (m s-1)'/
     *     5x,'4. v Stokes velocity = ',g12.3,' (m s-1)'/
     *     5x,'5.     u wind stress = ',g12.3,' (N m-2)'/
     *     5x,'6.     v wind stress = ',g12.3,' (N m-2)')
c
      endif
c
c set up the random cooling
c 
      if(nstep.gt.nrndom.or.(.not.bcuseran))then
       if(iope)then
         do 16 j=1,jtotal
            do 17 i=1,itotal
               phi2(i,j)=1.0
 17         continue
 16      continue
       endif
#if defined(USEMPI)
       call distrib_2d(randt,phi2)
#else
       do j=0,jdim1
         do i=0,idim1
           randt(i,j) = phi2(i,j)
         enddo
       enddo
#endif
      else
       if(iope) then
         do 18 j=1,jtotal
            do 10 i=1,itotal
               phi2(i,j)=(1.0-ranamp)+ranamp*ran2(nf)
 10         continue
c
c filter in x-direction
c
            do 11 i=1,itotal-1
               phi2(i,j) = (2.*phi2(i,j)+phi2(i+1,j)+phi2(i-1,j))/4.
 11         continue
 18      continue
c
c filter in y-direction
c
         do 30 i=0,itotal
            do 31 j=2,jtotal-1
               phi2(i,j) = (2.*phi2(i,j)+phi2(i,j+1)+phi2(i,j-1))/4.
 31         continue
 30      continue
       endif
#if defined(USEMPI)
       call distrib_2d(randt,phi2)
#else
       do j=0,jdim1
         do i=0,idim1
           randt(i,j) = phi2(i,j)
         enddo
       enddo
#endif
      do j=1,jdim
         do i=1,idim
           randn(i,j) = randt(i,j)
         enddo
      enddo
      call imag2d(randn)

c       if(nprocs.gt.1)
c     *   write(*,*)"randn, mype ",my_pe,(randn(idim1,j),j=1,jdim)
      endif
c
#if !defined(RADIATE)
c 
c force the rigid bottom boundary condition
c 
      do 20 j=0,jdim1
         do 21 i=0,idim1
            w(i,j,1)=0.0
            wm(i,j,1)=0.0
            wmm(i,j,1)=0.0
 21      continue
 20   continue
#endif
#if !defined(RADTOP)
c
c force the rigid lid top boundary condition
c
      do j=0,jdim1
         do i=0,idim1
            w(i,j,kdim1)=0.0
            wm(i,j,kdim1)=w(i,j,kdim1)
            wmm(i,j,kdim1)=w(i,j,kdim1)
         enddo
      enddo
#endif

#if defined(RADIATE) || defined(RADTOP)
c
c setup power of 2 factor for fft routines
c
      nstage = nint(log(float(idim))/log(2.))
c
      if(iope) then
        write(lunit, *)'Using radiative lower boundary'
        write(lunit, *)'nstage is ',nstage
      endif
c
c initialize lower boundary wavenumber array
c
c do for x axis
c
      do 40 i=2,idim/2
          rki(i) = abs(2*sin(rpi*float(i-1)/float(idim))/dx)
          rki(idim-i+2) = abs(2.*sin(rpi*float(i-1)/float(idim))/dx)
40    continue
      rki(1) = rki(2)
      rki(idim/2+1) = abs(2/dx)
c
c do for y axis
c
      do 45 j=2,jtotal/2
          rkj(j) = abs(2*sin(rpi*float(j-1)/float(jtotal))/dy) 
          rkj(jtotal-j+2) = abs(2.*sin(rpi*float(j-1)/float(jtotal))/dy)
45    continue
      rkj(1) = rkj(2)
      rkj(jtotal/2+1)=abs(2/dy)
c
c compute the two-dimensional wavenumber
c
      do 50 j=1,jdim
         do 51 i=1,idim
            js = blockstart(my_pe+1)+j-1
            rk(i,j) = sqrt(rki(i)**2+rkj(js)**2)
51       continue
50    continue
#endif
#if defined(RADIATE)
c
c compute brunt vaisala value for lower boundary
c
      drhodz = (rhob(2)-rhob(1))/dz
      bvf = sqrt(-grav*drhodz*2./(rhob(1)+rhob(2)))
      if(iope) then
        write(lunit,*)'Bottom Brunt-Vaisala Frequency is ',bvf
      endif
c
c compute background dtdz and dsdz (low)
c
      dtdzlow = th(1,1,2)-th(1,1,1)
      dsdzlow = sal(1,1,2)-sal(1,1,1)
c
c also set th at bottom with constant gradient condition
c
      do 60 i=0,idim1
      do 60 j=0,jdim1
        th(i,j,0) = th(i,j,1)-(th(i,j,2)-th(i,j,1))
60    continue
c
#endif
#if defined(RADTOP)
c
c compute brunt vaisala value for upper boundary
c
      drhodz = (rhob(kdim)-rhob(kdim-1))/dz
      bvft = sqrt(-grav*drhodz*2./(rhob(kdim)+rhob(kdim-1)))
#if defined(USEMPI)
      if(iope) then
        write(*,*)'Top Brunt-Vaisala Frequency is ',bvft
      endif
#else
        write(*,*)'Top Brunt-Vaisala Frequency is ',bvft
#endif
c
c compute background dtdz and dsdz (low)
c
c      dtdzlow = th(1,1,2)-th(1,1,1)
      dtdztop = th(1,1,kdim)-th(1,1,kdim-1)
      dsdztop = sal(1,1,kdim)-sal(1,1,kdim-1)
c
c also set th at top with constant gradient condition
c
      do 60 i=0,idim1
      do 60 j=0,jdim1
        th(i,j,kdim1) = th(i,j,kdim)+(th(i,j,kdim)-th(i,j,kdim-1))
60    continue
c
#endif
      return
2000  write(*,999)value(:30)
999   format('Units error on B.C. card for entry: ',a)
      ier=1
      return
      end
