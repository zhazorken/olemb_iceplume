      subroutine ic_read(ifile)
c
c
c   $Id: ic_read.F,v 1.11 1998/10/21 21:56:37 dwd Exp $
c
c     This routine reads initial condition (IC) or restart (RS) file.
c     The only input argument is the name of the cdf file output.  
c     
c     It is assumed that at the point that this routine is called, that some 
c     parameters (idim,jdim,kdim,dx,dy, and dz) have already been set 
c     (probably from rdparm).  These values are checked with those found in the
c     IC file(i.e. values read in in this routine). FILE_TYPE is 'START' for a
c     IC cdf file, and FILE_TYPE is 'RESTART'for a restart file.  
c     If the compatibility test fails, program execuation is stopped with the 
c     following error codes (ier = 0 indicates file compatibility):
c
c     1 = field time not equal to zero
c     2 = dx(rdparm) not equal to dx(IC)
c     3 = dy(rdparm) not equal to dy(IC)
c     4 = dz(rdparm) not equal to dz(IC)
c     5 = imax(rdparm) not equal to imax(IC)
c     6 = jmax(rdparm) not equal to jmax(IC)
c     7 = kmax(rdparm) for equal to kmax(IC)
c
c     Flow control in this program is dependent upon 1) the results of the 
c     compatibility test, 2) the value of the cdf attribute FILE_TYPE 
c     (stored in character string ftype), 3) and the number of tracers (ntrace)
c
c     Failing the compatibility test results in program termination. Control
c     based on FILE_TYPE occurs periodically throughout program.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "moddef.inc"
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "olemi.inc" 
#include "olemtke.inc"
#include "netcdf.inc"
#include "cyclic.inc"
#include "olemcdf.inc"
#include "refdate.inc"
#include "fftw_f77_param.h"
#include "communicate.inc"
#include "ice.inc"
#include "finitv.inc"
c
      character ifile*80,name*15,hold*2,titlei*80
      integer cdf,itrc,ntracei,junk,ierr
      integer lnblk,j,ilen,i,k,ier,n
      integer rcode,varid(4),varidr(4),varidm(4)
      integer xwdim,zudim,yudim,tdim,xdim,ydim,zdim
      integer salid,thid,rhoid,rhopid,thbarid,salbarid
      integer uid,umid,ummid
      integer wid,wmid,wmmid
      integer vid,vmid,vmmid
#if defined(ICE)
      integer hiceid,ticeid,siceid,cfrid,vbrid
      integer stx,endx
#endif
#if defined(HBOT)
      integer hbotid
#endif
      integer tkid,prid,khid,kmid,trid(trcdim)
      integer lf_tid, hf_tid, str_tid, str_tbcid
      integer ictid,sindx
      integer imax,jmax,kmax,ii,jj
c
c temporary space for MPI covering total domain
c
      real phi(0:idim1,0:jdim1,0:kdim1)
      real phi2(itotal,jtotal)
      real htop,hwid,r,thinmin,vz
      real ang,rnew,xof
      real phi1d(1:kdim)
      real h_mtn,a_mtn
      real*8 drhol,drho
      character tunits*80,str*20


      integer coun(4),star(4),counr(4),counm(4)
      real starti,iestbl,s,zl(2*kdim)
#if !defined(SNGLSALT)
c      real ssal(0:idim1,0:jdim1,0:kdim1)
#endif

      character*80 bccom1t,bccom2t,bccom3t,sgscom1t,sgscom2t,sgscom3t
      character ifilet*80,bfilet*80
      integer mdflgst(20)
      real mdvalst(20)
      integer notimest,nofldt,nostartt,nostopt,nodtt,ndiagt
      integer diloct,djloct,dkloct,ns

      logical uread,vread,wread,tkeread
#if defined(ICE)
      logical hiceread,ticeread,siceread,cfrread,vbrread
#endif
#if defined(HBOT)
      logical hbotread
#endif
#if defined(ANALY)
      real*8 pel(kdim),zz
      real den(0:idim1,0:jdim1,0:kdim1)
#endif
      data uread,vread,wread,tkeread/.false.,.false.,.false.,.false./
      data hfton,lfton,strton/.false.,.false.,.false./
c
c set ice ridge values
c
#if defined(ICE)
c        htop = 3.0
c        htop = 18.0
        htop = 5.0
c        htop = 15.0
c          htop = 1.0
c        hwid = 35.0
        hwid = 25.0
c        hwid = 10.
        write(lunit,*)"Ridge height,width ",htop,hwid
#endif
c
c-------Open cdf file.
c
#if defined(PRINTDEBUG)
      write(*,*)"Opening initialization file"
#endif
      cdf=ncopn(ifile,ncnowrit,rcode)

c-------Get global attributes from cdf ic or rst file.
c-------Note: variables ntracei are
c-------  essentially discarded after the test for compatability
c-------  succeeds, and their "non-i" equivalents are used from
c-------  then on.
c
      call ncagtc(cdf,NCGLOBAL,'FILE_TYPE',ftype,7,rcode)
      call ncagtc(cdf,NCGLOBAL,'title',titlei,80,rcode) 
c
c define grid spacing and size here
c
      call ncagt(cdf,NCGLOBAL,'dx',dx,rcode)
      call ncagt(cdf,NCGLOBAL,'dy',dy,rcode)
      call ncagt(cdf,NCGLOBAL,'dz',dz,rcode)
      call ncagt(cdf,NCGLOBAL,'imax',imax,rcode)
      call ncagt(cdf,NCGLOBAL,'jmax',jmax,rcode)
      call ncagt(cdf,NCGLOBAL,'kmax',kmax,rcode)
c
      call ncagt(cdf,NCGLOBAL,'field_time',starti,rcode)
      call ncagt(cdf,NCGLOBAL,'ntrace',ntracei,rcode)
      call ncagtc(cdf,NCGLOBAL,'creation_date',itime,30,rcode)
      call ncagtc(cdf,NCGLOBAL,'iccom_1',com1,120,rcode)
      call ncagtc(cdf,NCGLOBAL,'iccom_2',com2,120,rcode)
      call ncagtc(cdf,NCGLOBAL,'iccom_3',com3,120,rcode)
      call ncagtc(cdf,NCGLOBAL,'version',iver,80,rcode)
      do itrc=1,ntracei
         write(hold,'(i1)')itrc
         call ncagtc(cdf,NCGLOBAL,'trccom_'//hold,trccom(itrc),120,
     &               rcode)
         call ncagtc(cdf,NCGLOBAL,'trcfile_'//hold,trfile(itrc),80,
     &               rcode)
      enddo
c Are there any other global attributes we want??
#if defined(PRINTDEBUG)
      write(*,*)"Global attributes are read ",my_pe
#endif

      if(ftype.eq.'START')ititle=titlei

c------Test for compatibility.
      ier=0
      if(starti.ne.0.0.and.ftype.eq.'START') then
         write(*,*)' initial time not equal to zero.'
         ier=1
      endif
c
c test for cubic grid, define generic spacing 'ds'
c
      ds=dx
      if(dx.ne.dy .or. dy.ne.dz)then
         write(*,*)' input grid not cubic',dx,dy,dz
         ier=2
      endif
c
      if(imax.ne.idim*nprocx.or.jmax.ne.jdim*nprocy.or.kmax.ne.kdim)
     *  then
         write(*,*)' imax, jmax, kmax not equal to dimensioned size'
         write(*,*)imax,nprocx*idim,jmax,nprocy*jdim,kmax,kdim
         ier=5
      endif
c
      if(iope) then
        write (lunit,30) dx,dy,dz,imax,jmax,kmax
#if defined(PRINTDEBUG)
        write (*,30) dx,dy,dz,imax,jmax,kmax
#endif
 30     format(/'Geometry from I.C. file:'/
     *     5x,'  dx = ',f12.1,' (meters)'/
     *     5x,'  dy = ',f12.1,' (meters)'/
     *     5x,'  dz = ',f12.1,' (meters)'/
     *     5x,'imax = ',i12/
     *     5x,'jmax = ',i12/
     *     5x,'kmax = ',i12/)
c
        if(ntrace.ne.ntracei) then
           write(*,*)' ntrace (from input file): ',ntrace
           write(*,*)' ntrace (from initial condition file) : ',ntracei
           ier=8
        endif
      endif

c------Stop program execution if incompatible .in and .ic inputs.
      if(ier.ne.0) then
         write(7,*)' initial condition file is incompatible with inputs'
         write(7,*)' ier : ',ier
         if(ier.eq.5) then
            write(7,40) idim,jdim,jdim
 40         format(/'Dimensioned geometry:'/
     *           5x,'  idim = ',i5/
     *           5x,'  jdim = ',i5/
     *           5x,'  kdim = ',i5/)
         endif
#if defined(USEMPI)         
         call MPI_ABORT(ier)
#endif
      endif
c
c setup the table driven equation of state
c
      do k=1,kdim
         zl(2*k)=(kdim-k+0.5)*dz
         zl(2*k-1)=(kdim-k+1.0)*dz
      enddo
c
      call mkies(thtamin,thtamax,salmin,salmax,zl,2*kdim)
c
c------Define dimensions.
c
      zudim=ncdid(cdf,'zu',rcode)
      xwdim=ncdid(cdf,'xw',rcode)
      yudim=ncdid(cdf,'yu',rcode)
      tdim=ncdid(cdf,'time',rcode)
      varid(1)=xwdim     
      varid(2)=yudim
      varid(3)=zudim
      varid(4)=tdim

      xdim=ncdid(cdf,'x',rcode)
      ydim=ncdid(cdf,'y',rcode)
      varidr(1)=xdim     
      varidr(2)=ydim
      varidr(3)=zudim
      varidr(4)=tdim

      zdim=ncdid(cdf,'z',rcode)
      varidm(1)=xwdim
      varidm(2)=yudim
      varidm(3)=zdim
      varidm(4)=tdim
#if defined(PRINTDEBUG)
      write(*,*) "Dimensions defined "
#endif
c
c read units and try to find reference time
c
      ictid = ncvid(cdf, 'time', rcode)
      icdate = .false.
c
      call ncagtc(cdf,ictid,'units', tunits, 80, rcode)
      sindx = index(tunits,'since')
      if(sindx .gt. 0) then
         do i=sindx+5,len(tunits)
            if(tunits(i:i).ne.' ') then
               call parsedate(tunits(i:), icref)
               icdate = .true.
               goto 100
            endif
         enddo
      endif
 100  if(icdate) then
         if(iope) then
           call frmtdate(icref, str)
           write(lunit,*) 'I.C. reference date = ', str
         endif
      endif
c
c------Get variable id's from cdf file (START or RESTART).
      call varnme(MDTEMP,name)
      ilen=lnblk(name,15)
      thid=ncvid(cdf,name(:ilen),rcode)

      call varnme(MDSAL,name)
      ilen=lnblk(name,15)
      salid=ncvid(cdf,name(:ilen),rcode)

      call varnme(MDRHOBM1,name)
      ilen=lnblk(name,15)
      rhoid=ncvid(cdf,name(:ilen),rcode)
#if defined(SPONGE)
      call varnme(MDSALBAR,name)
      ilen=lnblk(name,15)
      salbarid=ncvid(cdf,name(:ilen),rcode)

      call varnme(MDTHBAR,name)
      ilen=lnblk(name,15)
      thbarid=ncvid(cdf,name(:ilen),rcode)
#endif
c
#if defined(TRACER)
      do itrc=1,ntrace
         call varnme(MDTRCR+itrc,name)
         ilen=lnblk(name,15)
         trid(itrc)=ncvid(cdf,name(:ilen),rcode)
      enddo
#endif

c------Get/check for existence of initial u, v, w, or tke variables.
c------For example: if uread is true then the var u has already had
c------  an initial condition assigned to it in the ic cdf file, and
c------  so we should read those values in, not reinitialize it.

c------Turn off netCDF error checking.
      call ncpopt(0)

      call varnme(MDUVEL,name)
      ilen=lnblk(name,15)
      uid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         uread=.true.
      endif

      call varnme(MDWVEL,name)
      ilen=lnblk(name,15)
      wid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         wread=.true.
      endif

      call varnme(MDVVEL,name)
      ilen=lnblk(name,15)
      vid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         vread=.true.
      endif

      call varnme(MDTKE,name)
      ilen=lnblk(name,15)
      tkid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         tkeread=.true.
      endif

      call varnme(MDHFTMASK,name)
      ilen=lnblk(name,15)
      hf_tid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         hfton=.true.
      endif

      call varnme(MDLFTMASK,name)
      ilen=lnblk(name,15)
      lf_tid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         lfton=.true.
      endif

      call varnme(MDSTRTMASK,name)
      ilen=lnblk(name,15)
      str_tid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         strton=.true.
      endif

      call varnme(MDSTRTBC,name)
      ilen=lnblk(name,15)
      str_tbcid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         strton=.true.
      endif
#if defined(HBOT)
      call varnme(MDHBOT,name)
      ilen=lnblk(name,15)
      hbotid=ncvid(cdf,name(:ilen),rcode)
      if(rcode.ne.0) then
        hbotread=.true.
      endif
#endif

#if defined(ICE)

      call varnme(MDHICE,name)
      ilen=lnblk(name,15)
      hiceid=ncvid(cdf,name(:ilen),rcode)
      if(rcode.ne.0) then
        hiceread=.true.
      endif

      call varnme(MDSICE,name)
      ilen=lnblk(name,15)
      siceid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         siceread=.true.
      endif

      call varnme(MDTICE,name)
      ilen=lnblk(name,15)
      ticeid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         ticeread=.true.
      endif

      call varnme(MDVBR,name)
      ilen=lnblk(name,15)
      vbrid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         vbrread=.true.
      endif

      call varnme(MDCFR,name)
      ilen=lnblk(name,15)
      cfrid=ncvid(cdf,name(:ilen),rcode)
      if (rcode.eq.0) then
         cfrread=.true.
      endif
#endif

      
c------Turn netCDF error checking back on.
      call ncpopt(NCVERBOS+NCFATAL)

c------Write information to log file.
      if(ftype.eq.'START') then
         if(iope) then
           write(lunit,*)'*********START file*********'
#if defined(PRINTDEBUG)
           write(*,*)'*********START file*********'
#endif
         endif
         start=0.0
         isrestart=.false.

      elseif(ftype.eq.'RESTART') then
         write(*,*)"Reading restart file"
         if(iope) then
           write(lunit,*)'************RESTART file*****************'
         endif
         start=starti
         isrestart=.true.

c---------Set u,v,w,tke variables to true.  (We should read them
c--------- as values exist for them in the restart cdf file.)
c         tkeread=.true.

c---------Get id's of restart variables.
         call varnme(MDUM,name)
         ilen=lnblk(name,15)
         umid=ncvid(cdf,name(:ilen),rcode)
 
c         call varnme(MDUMM,name)
c         ilen=lnblk(name,15)
c         ummid=ncvid(cdf,name(:ilen),rcode)
c
         call varnme(MDWM,name)
         ilen=lnblk(name,15)
         wmid=ncvid(cdf,name(:ilen),rcode)

c         call varnme(MDWMM,name)
c         ilen=lnblk(name,15)
c         wmmid=ncvid(cdf,name(:ilen),rcode)

         call varnme(MDVM,name)
         ilen=lnblk(name,15)
         vmid=ncvid(cdf,name(:ilen),rcode)

c         call varnme(MDVMM,name)
c         ilen=lnblk(name,15)
c         vmmid=ncvid(cdf,name(:ilen),rcode)

         call varnme(MDPRES,name)
         ilen=lnblk(name,15)
         prid=ncvid(cdf,name(:ilen),rcode)

c         call varnme(MDKH,name)
c         ilen=lnblk(name,15)
c         khid=ncvid(cdf,name(:ilen),rcode)

c         call varnme(MDKM,name)
c         ilen=lnblk(name,15)
c         kmid=ncvid(cdf,name(:ilen),rcode)

c         call varnme(MDRHO,name)
c         ilen=lnblk(name,15)
c         rhopid=ncvid(cdf,name(:ilen),rcode)
      endif
 
c------Set density identification.
      star(1)=1
      star(2)=1
      star(3)=1
      star(4)=1

      coun(1)=itotal
      coun(2)=jtotal
      coun(3)=kdim
      coun(4)=1

      counm(1)=itotal
      counm(2)=jtotal
      counm(3)=1
      counm(4)=1

      counr(1)=1
      counr(2)=1
      counr(3)=kdim
      counr(4)=1
c
      call read3d(cdf,thid,th)
c
#if defined(SNGLSALT)
      call read3d(cdf,salid,sal)
#else
      call read3ddp(cdf,salid,sal)
#endif
c
#if defined(TRACER)
      do itrc=1,ntrace
        call read3d(cdf,trid(itrc),trid(itrc))
      enddo
#endif
c
      if(uread) then
        call read3d(cdf,uid,u)
      endif
c
      if(vread) then
        call read3d(cdf,vid,v)
      endif
c
      if(wread) then
        call read3d(cdf,wid,w)
      endif
c
      if(tkeread) then
        call read3d(cdf,tkid,tke)
      endif
c #if defined(PRINTDEBUG)
      write(*,*)"basic fields read in ",my_pe
c #endif

c
c------Reinitialize star.
      star(1)=1
      star(2)=1
      star(3)=1
      star(4)=1
      call ncvgt(cdf,rhoid,star,counr,phi1d,rcode)
      do k=1,kdim
        rhobm1(k) = phi1d(k)
      enddo
#if defined(SPONGE)
      call ncvgt(cdf,salbarid,star,counr,salbar(1),rcode)
      call ncvgt(cdf,thbarid,star,counr,thbar(1),rcode)
#endif
c------Read in the masks
      if (hfton) then
        call read2d(cdf,hf_tid,star,counm,hf_tmask)
        call imag2d(hf_tmask)
      endif
      if (lfton) then
        call read2d(cdf,lf_tid,star,counm,lf_tmask)
        call imag2d(lf_tmask)
      endif
      if (strton) then
        call read2d(cdf,str_tid,star,counm,str_tmask)
        call imag2d(str_tmask)
        call read2d(cdf,str_tbcid,star,counm,str_tbc)
        call imag2d(str_tbc)
      endif
      do k=0,kdim1
        raytau(k) = 0.0
      enddo
      if(useray) then
#if defined(RADTOP)
           write(lunit,*)"Using Rayleigh Damping for upper b.c."
           write(lunit,*)"Relaxation time at top is",taumax
           do n=1,10
             raytau(kdim-10+n) = (sin(1.5708*n/10)**2)/taumax
           enddo
           write(lunit,*)"Rayleigh Time Array ",(raytau(kdim-9+n),n=1,9)
#else
           write(lunit,*)"Radtop set, no Rayleigh damping on top"
#endif
#if defined(RADIATE)
           write(lunit,*)"Using Rayleigh Damping for lower b.c."
#if defined(PRINTDEBUG)
           write(*,*)"Using Rayleigh Damping for lower b.c."
#endif
           write(lunit,*)"Relaxation time at bot is",taumax
           do k=0,kdim1
             raytau(k) = 0.0
           enddo
           do n=1,8
             raytau(9-n) = (sin(1.5708*n/8)**2)/taumax
           enddo
           write(lunit,*)"Rayleigh Time Array ",(raytau(9-n),n=1,9)
#endif
       endif


c only if START
      if(ftype.eq.'START') then
c  
c---------Initialize variables, if data does not already exist (u,v,w,tke).
c
c
         write(*,*) "Initialize variables uvwtke "
c
         call MPI_BARRIER(comm,ier)

         do j=1,jdim
            do 145 k=0,kdim1
               do 140 i=1,idim
                  if (.not.uread) u(i,j,k)=0.0
                  if (.not.wread) w(i,j,k)=0.0
c                  if (.not.tkeread) tke(i,j,k)=0.0
                  if (.not.vread) v(i,j,k)=0.0
                  p(i,j,k)=0.0
                  pm(i,j,k)=0.0
                  pmm(i,j,k)=0.0
                  um(i,j,k)=u(i,j,k)
                  umm(i,j,k)=u(i,j,k)
                  vm(i,j,k)=v(i,j,k)
                  vmm(i,j,k)=v(i,j,k)
                  if(k.ne.0) then
                     wm(i,j,k)=w(i,j,k)
                     wmm(i,j,k)=w(i,j,k)
                  endif
 140           continue
 145        continue
         enddo

         call topbotb(th,MDTEMP,ns)
#if defined(SNGLSALT)
         call topbotb(sal,MDSAL,ns)
#else 
         call topbotdp(sal,MDSAL,ns)
#endif

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
         call imag(p)
         call imag(pm)
         call imag(pmm)
c
         write(*,*) "Imag points set "
c
         call MPI_BARRIER(comm,ier)

#if defined(FRONT)
c
c         if(iope) then
           thinmin = th(1,1,kdim)-0.25
           ml = kdim-40./dz
           do k=0,kdim1
c             if(k.gt.ml) then
c               xof = (500.*(k - ml)/(kdim-ml))/dx
c               if(iope) write(*,*)"xof = ",xof,k
c             endif
             do j=1,jdim
               do i=1,idim
                
c
                 ii = iblkstart(my_pe+1)+i
                 jj = jblkstart(my_pe+1)+j
c
c double front warm in the middle
c
               if(k.gt.ml) then
                 coef = 1.0
               else
                 coef = exp((k-ml)/5.)
               endif
               if(ii.lt.itotal/2) then
                 phi(i,j,k) = th(i,j,k)+
     *               0.06*tanh((real(ii-itotal/3.))/80.)*coef
c     *               0.04*tanh((float(i-idim/3)+xof)/100.)
c     *               0.2*tanh(float(i-idim/4)/100.)
               else
                 phi(i,j,k) = th(i,j,k)-
     *            0.06*tanh((real(ii-2.*itotal/3.))/80.)*coef
c     *             0.04*tanh((float(i-idim/2-(idim/2-idim/3))-xof)/100.)
c     *             0.2*tanh(float(i-idim/2-idim/4)/100.)
               endif
c               if(i+iblkstart(my_pe+1).gt.509.and.
c     *            i+iblkstart(my_pe+1).lt.514
c     *             .and.k.ge.kdim-2) then
c                 tke(i,j,k) = 1.0
c               else if(i.gt.330.and.i.lt.340
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               else
c                 tke(i,j,k) = 0.0
c               endif
               enddo
             enddo
           enddo
c
c turn off front
c
           do k=0,kdim1
           do j=1,jdim
           do i=1,idim
             if(phi(i,j,k).gt.thinmin) then
               th(i,j,k) = phi(i,j,k)
             endif
           enddo
           enddo
           enddo
c         call imag(tke)
#endif
         call imag(th)

         rho0=0.0
         do 105 k=1,kdim
            drhol=0.0
            do 106 j=1,jdim
               do 107 i=1,idim
                  s=sal(i,j,k)
                  drhol=drhol+iestbl(s,th(i,j,k),2*k)
 107           continue
 106        continue
            call MPI_ALLREDUCE(drhol,drho,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,comm, ierr)

            rhobm1(k)=drho/(jtotal*itotal)
            rhob(k)=1000.0+rhobm1(k)
            if(iope) then
              write(*,*)"k,rhob,t,s ",k,rhob(k),th(1,1,k),sal(1,1,k)
            endif
            rho0=rho0+rhob(k)
 105     continue
         rhob(0)=rhob(1)
         rhob(kdim1) = rhob(kdim)
         rho0=rho0/float(kdim)
         do j=1,jdim
            do 161 k=1,kdim
               do 162 i=1,idim
                  s=sal(i,j,k)
                  rhop(i,j,k)=iestbl(s,th(i,j,k),2*k)-rhobm1(k)
 162           continue
 161        continue
            
          do 163 i=1,idim
               rhop(i,j,0)=rhop(i,j,1)
               rhop(i,j,kdim1)=rhop(i,j,kdim)
 163        continue
         enddo
c         call denupdte
         call imag(rhop)
         do i=0,idim1
           do j=0,jdim1
             vm(i,j,0) = 0.0
           enddo
         enddo
#if defined(FRONT)
c
c apply thermal wind
c
         do k=1,kdim
           do j=1,jdim
             do i=1,idim
               ii = iblkstart(my_pe+1)+i
               jj = jblkstart(my_pe+1)+j
               if(ii.lt.itotal/2) then
                 vz = grav*
     *               0.06/cosh((real(ii-itotal/3.))/80.)**2
               else
                 vz = -grav*
     *               0.06/cosh((real(ii-2.*itotal/3))/80.)**2
               endif
               vz = 2.e-4*vz/(rhob(k)*tomegaz)
               vm(i,j,k) = vm(i,j,k-1)+vz*dz
               v(i,j,k) = vm(i,j,k)
               vmm(i,j,k) = vm(i,j,k)
             enddo
           enddo
         enddo 
#endif
         if(iope) then
           write(*,*)"vm ",(vm(i,1,1),i=1,idim)
           write(*,*)"th ",(th(i,1,1),i=1,idim)
           write(*,*)"sal ",(sal(i,1,1),i=1,idim)
           write(*,*)"rhop ",(rhop(i,1,1),i=1,idim)
         endif
         call imag(v)
         call imag(vm)
         call imag(vmm)
         do j=1,jdim
            do 180 i=1,idim
               p(i,j,0)=0.0
               do 170 k=1,kdim1
                  p(i,j,k)=p(i,j,k-1)-
     &                 0.5*grav*(rhop(i,j,k)+rhop(i,j,k-1))*dz/rhob(k)
 170           continue
 180        continue
         enddo
         call imag(th)
         call imag(rhop)
         call imag(p)
#if defined(SNGLSALT)
         call imag(sal)
#else
         call imagdp(sal)
#endif
#if defined(HBOT)
         if(iope) then
           h_mtn = 20.
           a_mtn = 500.
           do j=1,jtotal
             do i=1,idim
c               phi2(i,j) = 0.8+
c     *            0.4*sin(3.14159*float(i-1)/(idim/2.))*
c     *            cos(3.14159*float(j-1)/(jtotal/8.))
c                phi2(i,j) = 1./(1.+ (i-idim/2)**2/25.)
c                 phi2(i,j) = 4.0*dx*exp(dble(-((i-idim/2)/8.0)**2))
c
c isolated small bump
c
c                 phi2(i,j) = 6.*dx*exp(dble(-((j-jtotal/2)/16.0)**2))
c
c bigger wider bump
c
c                 phi2(i,j) = 14.*dx*exp(dble(-((j-jtotal/2)/28.0)**2))
c
c                 phi2(i,j) = 10.*dx*exp(dble(-((j-jtotal/2)/20.0)**2))
c witch of agnesi
c                 phi2(i,j) = 6.0*dx/(1.0+((j-jtotal/2)/20.)**2)
c     *               *(1.+0.2*sin(3.14159*float(i-1)/(idim/8)))
c
                 phi2(i,j) = (h_mtn*a_mtn**2)
     *                    /  (((j*dy-jtotal*dy*0.5)**2 + a_mtn**2))

c
c shelf break case
c
c                 phi2(i,j) = 5.*dx*exp(dble(-((j-jtotal/2)/100.0)**2))
c     *               *(1.+0.2*sin(3.14159*float(i-1)/(idim/8)))
c                 if(j.lt.jtotal/2) then
c                   phi2(i,j) = 5.*dx
c                 endif
c
c                phi2(i,j) = 0.0
c
c     *             *exp(dble(-((i-idim/2)/8.0)**2))
c
c do a square meltpond
c
c             phi2(i,j) = 4.*dz
c             if(i.le.4.or.i.ge.idim-3) then
c                phi2(i,j) = kdim*dz
c             endif
c             if(j.le.4.or.j.ge.jtotal-3) then
c                phi2(i,j) = kdim*dz
c             endif
cc
c do a circular meltpond
c
c             r = sqrt(float(i-idim/2)**2+float(j-jtotal/2)**2)
c             if(r.gt.(idim/2 - 4)) then
c             if(r.gt.315) then
c             if(r.gt.125) then
c               phi2(i,j) = kdim*dz
c             endif
c
c do an irregular meltpond
c
c               ang = atan(float(j-jtotal/2)/(float(i-idim/2)+1e-6))
c
c               rnew = float(idim/2) - 0.1*r*(1.5 - 0.5*
c     *           (1.5*cos(ang*25.)+
c     *           7.2*sin(ang*6.)+cos(ang*100.)+0.3*cos(ang*200.)))
c
c c4
c
c               rnew = float(idim/2) - 0.1*r*(1.5 - 0.5*
c     *           (0.0*cos(ang*25.)+
c     *           6.0*sin(ang*6.)))
c
c c3
c               rnew = float(idim/2) - 0.1*r*(1.5 - 0.5*
c     *           (2.0*cos(ang*25.)+
c     *           8.0*sin(ang*6.)))
cc
c c5
c
c               rnew = float(idim/2) - 0.1*r*(1.5 - 0.5*
c     *           (1.0*cos(ang*25.)+
c     *           7.0*sin(ang*6.)+ 1.0*cos(ang*70.)))
c
c               rnew = float(idim/2) - 0.1*r*(1.5 - 0.5*(cos(ang*25.)+
c     *           6.*sin(ang*8.)+1.8*cos(ang*35.)+1.2*cos(ang*200.)))
c
c               if(r.gt.(rnew-90)) then
c                 phi2(i,j) = kdim*dz
c               endif
c
c
             enddo
c             if(i.eq.2) then
c               write(*,*)(phi2(i,j),j=1,jtotal)
c               read(*,*) junk
c             endif
           enddo
         endif
#if defined(USEMPI)
         call distrib_2d(hbot,phi2)
#else
         do j=1,jdim
           do i=1,idim
             hbot(i,j) = phi2(i,j)
           enddo
         enddo
#endif
         call imag2d(hbot)
         call volbot(hbot)
         if(my_pe.eq.3) then
           do i=1,10
             write(*,*)"kdist ",kdist(i,1)
           enddo
         endif
#endif
c
c Initialize Ice field
c
#if defined(ICE)
         if(iope) then
           do j=1,jtotal
             do i=1,idim
c
c
c               phi2(i,j) = 0.4+
c     *            0.4*sin(3.14159*float(i-1)/(idim/4.))*
c     *            cos(3.14159*float(j-1)/(jtotal/32.))
c                 phi2(i,j) = htop*dx*exp(dble(-((j-jtotal/2)/hwid)**2))
c
                 phi2(i,j) = htop*dx/(1.0+((j-jtotal/2)/hwid)**2)
c     *                    *(1.+0.2*sin(3.1416*float(i-1)/(idim/4)))
c     *                    *cos(3.1416*float(i-1)/(idim/8)))
c     *                    +1.0 +
c     *            0.6*sin(3.14159*float(i-1)/(idim/4.))*
c     *            cos(3.14159*float(j-1)/(jtotal/4.))

c     *             *exp(dble(-((i-idim/2)/15.0)**2))
c                 phi2(i,j) = 0.1
c                 if(i.gt.13.and.i.lt.19) then
c                  phi2(i,j) = 0.2
c                 endif
             enddo
             if(i.eq.2) then
c               write(*,*)(phi2(i,j),j=1,jdim)
c               write(*,*)"Hice"
c               read(*,*) junk
             endif
           enddo
         endif
#if defined(USEMPI) 
         call distrib_2d(hice,phi2)
#endif
         do j=1,jdim
           do i=1,idim
#if !defined(USEMPI)
             hice(i,j) = phi2(i,j)
#endif
c
#if !defined(WARMICE)
             hice(i,j) = 1.0
c #else
c              hice(i,j) = 0.015
#endif
c             sice(i,j) = sal(i,j,kdim)
             sice(i,j) = 0.0
             tice(i,j) = th(i,j,kdim)
             vbr(i,j) = 0.0
c
c if hf mask is on, then intialize a 1.5 m thick ice slab at the surface
c for non lead ice
c
c             if(hfton) then
c               if(hf_tmask(i,j).eq.0) then
c                 hice(i,j) = 2.0
c                 sice(i,j) = 0.0
c                 tice(i,j) = -1.5
c             else
c                 hice(i,j) = 0.001
c                 sice(i,j) = 0.0
c                 tice(i,j) = -1.5
c               endif
c             endif
           enddo
         enddo
c
c section for prescribing lead in summer situation
c
#if defined(ICEWALL)
        write(*,*)"Initialize lead "
        stx = idim/2 -  120
        endx = idim/2 + 120
        do j=0,jtotal1
          do i=0,idim1
            phi2(i,j) = 0.0
            if((i.ge.stx).and.(i.le.endx)) then
             if((j.ge.stx).and.(j.le.endx)) then
                phi2(i,j) = 1.0
             endif
            endif
          enddo
        enddo
#if defined(USEMPI)
        call distrib_2d(hf_tmask,phi2)
        call distrib_2d(lf_tmask,phi2)
        call distrib_2d(str_tmask,phi2)
#else
        do j=0,jdim1
          do i=0,idim1
            hf_tmask(i,j) = phi2(i,j)
            lf_tmask(i,j) = phi2(i,j)
            str_tmask(i,j) = phi2(i,j)
          enddo
        enddo
#endif
         strton = .true.
         hfton = .true.
         lfton = .false.
         call imag2d(hf_tmask)
         call imag2d(lf_tmask)
         call imag2d(str_tmask)
         do j=1,jdim
           do i=1,idim
             if(hf_tmask(i,j).eq.1) then
               hice(i,j) = 0.0
             else
               hice(i,j) = 1.5
               sice(i,j) = 4.0
               tice(i,j) = -0.054*sice(i,j)
             endif
           enddo
         enddo
#endif
         call imag2d(hice)
         call imag2d(tice)
         call imag2d(sice)
         call voltop(hice)
c
c set water properties to ice values when in ice
c
c         do i=1,idim
c           do j=1,jdim
c             do k=kdist(i,j)+1,kdim1
c               th(i,j,k) = tice(i,j)
c               sal(i,j,k) = sice(i,j)
c             enddo
c           enddo
c         enddo
c
#if !defined(WARMICE)
c
c initialize frazil volume fraction to zero
c
         do k=0,kdim1
           do j=0,jdim1
             do i=0,idim1
               cfr(i,j,k) = 0.0
             enddo
           enddo
         enddo
#endif
#endif
      endif
      

      if(ftype.eq.'RESTART') then
c
c
#if defined(PRINTDEBUG)
       write(*,*)"Reading restart file"
#endif
c---------Read in variable values.
c
c
c------Reinitialize star.
         star(1)=1
         star(2)=1
         star(3)=1
         star(4)=1
c
         call read3d(cdf,umid,um)
c
c         call read3d(cdf,ummid,umm)
c
         call read3d(cdf,vmid,vm)
c
c         call read3d(cdf,vmmid,vmm)
c
         call read3d(cdf,wmid,wm)
c
c         call read3d(cdf,wmmid,wmm)
c
         call read3d(cdf,prid,p)
c
c         call read3d(cdf,tkid,tke)
c
c         call read3d(cdf,kmid,kmom)
c
c         call read3d(cdf,khid,kheat)
c
c         call read3d(cdf,rhopid,rhop)
c          do k=0,kdim1
c             do j=0,jdim1
c               do i=0,idim1
c               if(i+iblkstart(my_pe+1).gt.509.and.
c     *            i+iblkstart(my_pe+1).lt.514
c     *             .and.k.ge.kdim-2) then
c                 tke(i,j,k) = 1.0
c               else
c                 tke(i,j,k) = 0.0
c               endif
c
c               if(i.gt.2100.and.i.lt.2110
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               else if(i.gt.2230.and.i.lt.2240
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               endif
c               if(i.gt.330.and.i.lt.340
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               else if(i.gt.500.and.i.lt.510
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               endif
c               tke(i,j,k) = 0.0
c               enddo
c             enddo
c           enddo

#if defined(PRINTDEBUG)
         write(*,*)"Standard Fields are in"
#endif
#if defined(HBOT)
         call read2d(cdf,hbotid,hbot,star,counm,phi2)
         call imag2d(hbot)
         call volbot(hbot)
         if(my_pe.eq.3) then
           do i=1,10
             write(*,*)"kdist ",kdist(i,3)
           enddo
         endif
#endif
c
#if defined(ICE)
c
#if !defined(WARMICE)
         if(cfrread) then
           call read3d(cdf,cfrid,cfr,phi)
         endif
#endif
c
         call read2d(cdf,hiceid,star,counm,hice,phi2)
c
         call read2d(cdf,ticeid,star,counm,tice,phi2)
c
         call read2d(cdf,siceid,star,counm,sice,phi2)
c
         if(vbrread) then
          call read2d(cdf,vbrid,star,counm,vbr,phi2)
         endif
c
         if(int(start).eq.36000) then
#if defined(PRINTDEBUG)
	    write(*,*)"Initialize Lead"
#endif
#if !defined(WARMICE)
	    do i=1,idim
	      do j=1,jdim
                if(hf_tmask(i,j).ne.0)then
		  hice(i,j) = 0.01
		  tice(i,j) = -1.51
		  sice(i,j) = 0.0
                  if(.not.vbrread) then
                   vbr(i,j)  = 0.0
                  endif
		endif

	      enddo
	    enddo
#endif
	 endif
c
#if !defined(WARMICE)
         call imag(cfr)
#endif
         call imag2d(hice)
         call imag2d(sice)
         call imag2d(vbr)
         call imag2d(tice)
         call voltop(hice)
#endif
c
c
         do j=0,jdim+1
            do i=0,idim+1
               randn(i,j)=1.0
            enddo
         enddo
         do k=1,kdim
           do j=0,jdim1
             do i=0,idim1
#if !defined(LEAPF)
c                p(i,j,k) = 0.0
#endif
c                px(i,j,k) = 0.0
                pm(i,j,k) = p(i,j,k)
                pmm(i,j,k) = p(i,j,k)
                u(i,j,k) = um(i,j,k)
                umm(i,j,k) = um(i,j,k)
                v(i,j,k) = vm(i,j,k)
                vmm(i,j,k) = vm(i,j,k)
                w(i,j,k) = wm(i,j,k)
                wmm(i,j,k) = wm(i,j,k)
             enddo
           enddo
         enddo
         rho0=0.0
         do 82 i=1,kdim
            rhob(i)=1000.0+rhobm1(i)
            rho0=rho0+rhob(i)
            if(iope) then
              write(*,*)"k,rhob,t,s ",i,rhob(i)
            endif
 82      continue
         rhob(0)=rhob(1)
         rhob(kdim1) = rhob(kdim)
         rho0=rho0/float(kdim)
      
c---------Get other restart cdf global attributes.
         call ncagt(cdf,NCGLOBAL,'heatf',heatf,rcode)
         call ncagt(cdf,NCGLOBAL,'ncool',ncool,rcode)       
         call ncagtc(cdf,NCGLOBAL,'bccom_1',bccom1t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'bccom_2',bccom2t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'bccom_3',bccom3t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'sgscom_1',sgscom1t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'sgscom_2',sgscom2t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'sgscom_3',sgscom3t,80,rcode)
         call ncagtc(cdf,NCGLOBAL,'IC_FILE_NAME',ifilet,80,rcode)
         call ncpopt(0)
         call ncagtc(cdf,NCGLOBAL,'BC_FILE_NAME',bfilet,80,rcode)
         call ncpopt(NCVERBOS+NCFATAL)
         call ncagt(cdf,NCGLOBAL,'mdflgs',mdflgst,rcode)
         call ncagt(cdf,NCGLOBAL,'mdvals',mdvalst,rcode)

         call ncagt(cdf,NCGLOBAL,'notimes',notimest,rcode)
         call ncagt(cdf,NCGLOBAL,'nofld',nofldt,rcode)
         call ncagt(cdf,NCGLOBAL,'nostart',nostartt,rcode)
         call ncagt(cdf,NCGLOBAL,'nostop',nostopt,rcode)
         call ncagt(cdf,NCGLOBAL,'nodt',nodtt,rcode)
         call ncagt(cdf,NCGLOBAL,'nofld',nofldt,rcode)
         call ncagt(cdf,NCGLOBAL,'ndiag',ndiagt,rcode)
         call ncagt(cdf,NCGLOBAL,'diloc',diloct,rcode)
         call ncagt(cdf,NCGLOBAL,'djloc',djloct,rcode)
         call ncagt(cdf,NCGLOBAL,'dkloc',dkloct,rcode)
c
c call topbotb
c
c compute starting index
c
         ns=int(start/delt+0.5)
c
#if defined(PRINTDEBUG)
         write(*,*)"ic_read, calling topbotb ",my_pe
#endif
c
         call topbotb(um,MDUVEL,ns)
         call topbotb(umm,MDUVEL,ns)
         call topbotb(vm,MDVVEL,ns)
         call topbotb(vmm,MDVVEL,ns)
         call topbotb(wm,MDWVEL,ns)
         call topbotb(wmm,MDWVEL,ns)
         call topbotb(th,MDTEMP,ns)
#if defined(SNGLSALT)
         call topbotb(sal,MDSAL,ns)
#else 
         call topbotdp(sal,MDSAL,ns)
#endif
         call topbotb(p,MDPRES,ns)
         call topbotb(kmom,MDKH,ns)
         call topbotb(kheat,MDKH,ns)
         call topbotb(tke,MDTKE,ns)
#if defined(PRINTDEBUG)
         write(*,*)"ic_read, topbotb done ",my_pe
#endif
c
         call imag(um)
         call imag(umm)
         call imag(vm)
         call imag(vmm)
         call imag(wm)
         call imag(wmm)
         call imag(th)
#if defined(SNGLSALT)
         call imag(sal)
#else
         call imagdp(sal)
#endif
         call imag(p)
         call imag(kmom)
         call imag(kheat)
         call imag(tke)
#if defined(PRINTDEBUG)
         write(*,*)"ic_read, imag done ",my_pe
#endif

c
         if(iope) then
           write(lunit,200) itime,iver

           k=lnblk(titlei,80)
           write(lunit,201)titlei(1:k)

           k=lnblk(ifilet,80)
           write(lunit,206)ifilet(1:k)
           k=lnblk(bfilet,80)
           write(lunit,207)bfilet(1:k)
           write(lunit,202)bccom1t,bccom2t,bccom3t

           write(lunit,211)sgscom1t,sgscom2t,sgscom3t

           write(lunit,214)mdflgst
           write(lunit,215)mdvalst
         endif
c
c initialize old pe
c
#if defined(ANALY)
c
c compute density perturbation, only works for constant salinity case
c
#if defined(SNGLSALT)
c           call density(th,sal,den)
#else
c           call densitydp(th,sal,den)
#endif
c
c           do k=1,kdim
c             pel(k) = 0.0
c             zz = (k-1)*dz + 0.5*dz
c             do j=1,jdim
c               do i=1,idim
c                 pel(k) = pel(k)+(den(i,j,k)+rhobm1(k))*zz*9.81/rhob(k)
c               enddo
c             enddo
c           enddo
#if defined(USEMPI)
c
c do processor summation
c
c           call MPI_REDUCE(pel,peo,kdim,
c     &             MPI_DOUBLE_PRECISION,MPI_SUM,mastertask,
c     &             comm, ier)
c
           call MPI_BARRIER(comm,ier)
#else
           do k=1,kdim
             peo(k) = pel(k)
           enddo
#endif
c
c ice
c
#endif

 200     format(' Restart from a previous run: ',/a/a)
 201     format(a)
 202     format(/'Original B.C. Parameters: ',/2x,a80/2x,a80/2x,a80)
 203     format('Original B.C. Flags:'/5(4i20/))
 204     format('Original B.C. Values:'/5(1p4g20.8/))
 206     format('Original I.C. file: ',a)
 207     format('Original B.C. File: ',a)
 211     format(/'Original SGS Parameters:'/2x,a80/2x/a80/2x,a80)
 212     format(/'Original SGS Flags:'/5(4i20/))
 213     format(/'Original SGS Values:'/5(1p4g20.8/))
c not sure of format statements 214 and 215
 214     format(/'Original Model Flags:'/5(4i20/))
 215     format(/'Original Model Values:'/5(1p4g20.8/))
      endif

      call ncclos(cdf,rcode)
      
#if defined(USEMPI) && defined(PRINTDEBUG)
      write(*,*)' done reading ', ifile(1:60),' pe number ',my_pe
#endif

      return
      end
c
c
c The following are i/o routines that split off mpi
c
#if !defined(SNGLSALT)
      subroutine read3ddp(cdfin,cdfid,start,count,fld,phi)
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
      real*8 fld(0:idim1,0:jdim1,0:kdim1)
      real phi(0:idim1,0:jtotal1,0:kdim1)
      integer cdfin,cdfid,start(4),count(4),rcode
c
      integer i,j,k
#if defined(USEMPI)
      real*8 phi8(0:idim1,0:jtotal1,0:kdim1)

      if(iope) then
        call ncvgt(cdfin,cdfid,start,count,phi,rcode)
      endif
      do k=0,kdim1
         do j=0,jtotal1
            do i=0,idim1
               phi8(i,j,k) = phi(i,j,k)
            enddo
         enddo
      enddo
      call distrib_3ddp(fld,phi8)
#else
      call ncvgt(cdfin,cdfid,start,count,phi,rcode)
      do k=0,kdim1
         do j=0,jtotal1
            do i=0,idim1
               fld(i,j,k) = phi(i,j,k)
            enddo
         enddo
      enddo
#endif
      return
      end
#endif
c
c
      subroutine read2d(cdfin,cdfid,start,count,fld,phi)
c
c this routine inputs a two dimensional field including the
c image points (for use in restart). If using MPI, it does a
c distribution to all processors.
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
c
#include "param.inc"
#include "communicate.inc"
c
      real fld(0:idim1,0:jdim1)
      real fldt(idim,jdim)
      real phi(itotal,jtotal)
      integer cdfin,cdfid,start(4),count(4),rcode
c
#if defined(USEMPI)
      if(iope) then
        call ncvgt(cdfin,  cdfid, start, count,   phi,rcode)
      endif
      call distrib_2d(fldt,phi)
#else
      call ncvgt(cdfin,  cdfid, start, count,   fldt,rcode)
#endif
      do j=1,jdim
        do i=1,idim
          fld(i,j) = fldt(i,j)
        enddo
      enddo
      return
      end

