      program olem
c 
c     $Id: olem.F,v 1.22 1999/08/16 20:22:32 eds Exp $
c 
c Oceanic Large Eddy Model (OLEM).
c
c by Eric Skyllingstad and Donald Denbo
c
c The model is a non-hydrostatic, Boussinesq, incompressible Navier-
c Stokes primitive equation three-dimensional model.  The model is 
c designed to be used in process studies of the upper ocean and for
c oceanic deep convection.  The model design uses vertical slabs and
c image points surrounding the model domain to simplify model structure.
c 
c model units are in mks.
c
c ***********************************************************************
c ***********************************************************************
c 
c a vertical cut (i,k) of the grid is:
c
c          0    0    1    1                      idim idim idim1 idim1 
c
c kdim1    u    I    u    I    u    I    u    I    u    I    u    I
c
c kdim1         w    +----w---------w---------w---------w----+    w kdim1
c                  +.|.....................................+ |
c kdim     u    I  . u    I    u    I    u    I    u    I  . u    I kdim
c                  . |                                     . |
c kdim          w  . |    w         w         w         w  . |    w kdim
c                  . |                                     . |
c          u    I  . u    I    u    I    u    I    u    I  . u    I
c                  . |                                     . |
c               w  . |    w         w         w         w  . |    w
c                  . |                                     . |
c          u    I  . u    I    u    I    u    I    u    I  . u    I
c                  . |                                     . |
c               w  . |    w         w         w         w  . |    w
c                  . |                                     . |
c 1        u    I  . u    I    u    I    u    I    u    I  . u    I
c                  +.|.................................... + |
c 1             w    +----w---------w---------w---------w----+    w
c                                                          
c 0        u    I    u    I    u    I    u    I    u    I    u    I
c
c 0             w         w         w         w         w         w
c
c          0    0    1    1                      idim idim idim1 idim1 
c 
c
c ***********************************************************************
c ***********************************************************************
c
c ***********************************************************************
c ***********************************************************************
c
c a horizontal cut (i,j) of the grid is:
c
c          0    0    1    1                      idim idim idim1 idim1 
c
c jdim1    u    I    u    I    u    I    u    I    u    I    u    I
c
c jdim1         v    +----v---------v---------v---------v----+    v jdim1
c                  +.|.....................................+ |
c jdim     u    I  . u    I    u    I    u    I    u    I  . u    I jdim
c                  . |                                     . |
c jdim          v  . |    v         v         v         v  . |    v jdim
c                  . |                                     . |
c          u    I  . u    I    u    I    u    I    u    I  . u    I
c                  . |                                     . |
c               v  . |    v         v         v         v  . |    v
c                  . |                                     . |
c          u    I  . u    I    u    I    u    I    u    I  . u    I
c                  . |                                     . |
c               v  . |    v         v         v         v  . |    v
c                  . |                                     . |
c 1        u    I  . u    I    u    I    u    I    u    I  . u    I
c                  . |                                     . |
c 1             v  . +----v---------v---------v---------v--.-+    v
c                  +...................................... + 
c 0        u    I    u    I    u    I    u    I    u    I    u    I
c
c 0             v         v         v         v         v         v
c
c          0    0    1    1                      idim idim idim1 idim1 
c 
c ***********************************************************************
c ***********************************************************************
c
c ***********************************************************************
c ***********************************************************************
c 
c a vertical cut (j,k) of the grid is:
c
c          0    0    1    1                      jdim jdim jdim1 jdim1 
c
c kdim1    v    I    v    I    v    I    v    I    v    I    v    I
c
c kdim1         w    +----w---------w---------w---------w----+    w kdim1
c                  +.|.....................................+ |
c kdim     v    I  . v    I    v    I    v    I    v    I  . v    I kdim
c                  . |                                     . |
c kdim          w  . |    w         w         w         w  . |    w kdim
c                  . |                                     . |
c          v    I  . v    I    v    I    v    I    v    I  . v    I
c                  . |                                     . |
c               w  . |    w         w         w         w  . |    w
c                  . |                                     . |
c          v    I  . v    I    v    I    v    I    v    I  . v    I
c                  . |                                     . |
c               w  . |    w         w         w         w  . |    w
c                  . |                                     . |
c 1        v    I  . v    I    v    I    v    I    v    I  . v    I
c                  +.|.................................... + |
c 1             w    +----w---------w---------w---------w----+    w
c                                                          
c 0        v    I    v    I    v    I    v    I    v    I    v    I
c
c 0             w         w         w         w         w         w
c
c          0    0    1    1                      jdim jdim jdim1 jdim1 
c 
c
c ***********************************************************************
c ***********************************************************************
c
c the grid starts and ends on u, v or w fields, I is the location of
c pressure, salinity, temperature and tracers.  The pressure field
c forcing is updated with small time steps, the advection etc with large
c time steps.  'k' is the vertical coordinate, k=1 at the bottom.  'i'
c and 'j' are the horizontal coordinates, i=1 on the left and j=1 at the
c lower corner.  w(k=0) is not used.  The solid box in the actual model
c domain.  The dotted box includes those variables that are actually
c computed.  The i=0,idim1, j=0,jdim1 and k=0,kdim1 points are image
c points and are computed by "image" and "topbotb".
c 
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemi.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
#include "olemtke.inc"
#include "olemtim.inc"
#include "cyclic.inc"
#include "float.inc"
#include "refdate.inc"
#include "signal.inc"
#include "communicate.inc"
c 
      logical rstart
      character pfile*80,ifile*80,lfile*80,froot*80,bcfile*80
      character ffile*80,modid*30,hostname*50,dummy*120,penum*1
      character*120 bc1,bc2,bc3
      integer lnblk,k,ier,i,nerr,ilen,ns,n,rcode,index
#if defined(USEMPI)
      integer mpierr
#endif
      integer staunit
      integer mlen,numprocs
      real cputime,fctr,gpsec
c
#if defined(sun)
      integer ieeer,myhandler,ieee_handler,siger,signal
      external myhandler,sig_handler
#endif
c
c
c 
#include "moddef.inc"
#include "version.inc"
#if defined(Win32)
      call deltime(cputime)
#endif
c
c setup file lun's
c
c first call setupmpi to get the MPI system up and running
c
#if defined(USEMPI)
      call mpi_init(mpierr)
      call setupmpi()
#else
      iope = .true.
#endif
c 
      punit=1
      iunit=2
      ounit=3
      dunit=4
      runit=8
      lunit=7
      funit=9
      staunit = 12+my_pe
c
c setup constants
c
c salt diffusivity
c
      kms=7.4e-10
c
      grav=9.81
c
c earths rotation 2pi/(24*60*60)
c
      omega=7.2722052e-5
c
c specific heat (W sec kg-1 C-1) at T=0-10 degC, S=35 PSU and P=0 
c
      cp = 3986.0
c
c latent heat of vaporization (W sec kg-1) at T=0 degC
c
      lv = 2.5008e6
c
#if defined(sun)
c
c set up ieee handler
c
      ieeer = ieee_handler('set','common',myhandler)
c
c set up signal handler
c
      sig_restart=.false.
      sig_quit=.false.
      sig_quick=.false.
      sig_usr1=.false.
      sig_usr2=.false.
c
c      call signal(15, sig_handler)
c      call signal(24, sig_handler)
c      call signal(30, sig_handler)
c      call signal(31, sig_handler)
c      call signal(8, sig_handler)
c      call signal(SIGTRAP,XL_TRCE)
c
#endif
c 
c set model unique identifier
c 
c      ident=modid()
      ident = "olem"
c
c get the name of the machine 
c
      call hostnme(hostname)
c set up model cpp and flags
c
      modflags(1) = 
     *'MODEL_FLAGS1'
      modflags(2) = 
     *'MODEL_FLAGS2'
      modflags(3) =
     *'MODEL_FLAGS3'
      modflags(4) =
     *'MODEL_FLAGS4'
c
c write out ident and dimensioned size
c
      do i=1,80
        pfile(i:i) = ' '
        ifile(i:i) = ' '
        bcfile(i:i) = ' '
        ffile(i:i) = ' '
        froot(i:i) = ' '
      enddo

#if defined(USEMPI)
      call MPI_Comm_size(comm, numprocs, nerr)
#else
      numprocs = 1
#endif

      if(iope) then
        write (*,5) ident, idim, jtotal, kdim, nprocs
 5      format (/' Run ident = ', a/
     *         ' Grid size = (',i3,'x',i3,'x',i3,') with ',
     *         i3,' processors'/)
c
        write(*,6)
 6      format(' Compilation options: ')
        do i=1,4
           mlen = lnblk(modflags(i),80)
           if(mlen.ge.1) then
              write(*,8)modflags(i)(:mlen)
           endif
 8         format(1x,a)
        enddo
        write(*,9)
 9      format(/)
        if(numprocs .ne. nprocs) then
           write(*,7) numprocs, nprocs
 7         format(' Number of processors (',i3,') not equal to ',i3)
#if defined(USEMPI)
           call MPI_FINALIZE(nerr)
#endif
           call exit(1)
        endif
c
c open script file containing names of files
c
#if defined(SPLUME)
        open(unit=21,file='inputp.scp')
#else
        open(unit=21,file='inputi.scp')
#endif
c     
c 
c read file names and open proper files
c 
        write (*,10)
 10     format ('           Enter parameter file name: ',$)
        read (21,20) pfile
 20     format (a80)
        write (*,30)
 30     format ('      Enter initialization file name: ',$)
        read (21,20) ifile
        bcfile = ' '
        write(*,35)
 35     format (' Enter boundary conditions file name: ',$)
        read (21,20) bcfile
        ffile = ' '
        write(*,36)
 36     format ('               Enter float file name: ',$)
        read (21,20) ffile
        write (*,40)
 40     format ('         Enter output file root name: ',$)
        read (21,20) froot
        k=lnblk(froot,80)
        lfile=froot(1:k)//'.log'
        open (lunit,file=lfile,status='new')
c        open (12,file=froot(1:k)//'_sta.dat')
      endif
#if defined(USEMPI)
      call broadcast_string(pfile,80,my_pe)
      call broadcast_string(ifile,80,my_pe)
      call broadcast_string(bcfile,80,my_pe)
      call broadcast_string(ffile,80,my_pe)
      call broadcast_string(froot,80,my_pe)
#endif
c
#if defined(ICE)
      if(my_pe.lt.6)then
        k=lnblk(froot,80)
        write(penum,'(i1)')my_pe
        lfile=froot(1:k)//penum//'_sta.dat'
        open (staunit,file=lfile,status='new')
      endif
#endif
c
      havefloats=ffile.ne.' '
c 
c read parameter file
c 
#if defined(PRINTDEBUG)
      if(iope .and. nprocs .gt. 1) then
        write(*,42)"File names:",
     *        pfile,ifile,bcfile,ffile,froot,lfile
c MJH
 42     format( /,A,6(/, 4x, A),/,'End file names.',/,/ )
      endif
#endif
c
c read parameter file
c
      ilen = lnblk(pfile,80)
#if defined(PRINTDEBUG)
      if(nprocs.gt.1)
     * write(*,41)"Opening parameter file in process ",pfile,my_pe
c MJH
 41   format(A, /, 4X, 'File: ', A, /, 4X, 'Processor: ', I2)
#endif

      open (punit,file=pfile(1:ilen),status='old')
      write(*,*)"calling in_read"
      call in_read (ier)
      close (punit)
      if (ier.ne.0)then
        write(*,*) 'rdparm: error reading parameter input file'
#if defined(USEMPI)
        call MPI_ABORT(comm,1,mpierr)
#endif
      endif

      rstart = nrstrt.ne.0
      nerr=0
c
      if(iope) then
        ilen=lnblk(mver,80)
        write (lunit,60) ident,mver(:ilen)
 60     format ('Oceanic Large Eddy Model 3-D version',t50,a30/a//
     *        'Adams-Bashforth time stepping, pois3d pressure solver')
        do i=1,4
           mlen = lnblk(modflags(i),80)
           if(mlen.ge.1) then
              write(lunit,8) modflags(i)(:mlen)
           endif
        enddo
c
        write (lunit,61) hostname
 61     format (/'        Running on host: ',a)
c
        write (lunit,64) pfile,ifile,bcfile,ffile,froot
 64     format ('         Parameter file: ',a50/
     *        'Initial conditions file: ',a50/
     *       'Boundary conditions file: ',a50/
     *                     'Float file: ',a50/
     *        '       Output root name: ',a50)
c
        write(*,*)"calling wrparm"
        call wrparm
      endif
c 
c read initialization fields
c 
      write(*,*)"Reading ic file"
      call ic_read (ifile)
#if defined(PRINTDEBUG)
      if(iope) write(*,*)"ic is in"
#endif
c 
c set up required constants
c 
      strtf(1)=1
      strtf(2)=1
      strtf(3)=1
      strtf(4)=1

      countf(1)=idim
      countf(2)=1
      countf(3)=1
      countf(4)=1

      countuf(1)=idim1
      countuf(2)=1
      countuf(3)=1
      countuf(4)=1

      countvf(1)=idim
      countvf(2)=1
      countvf(3)=1
      countvf(4)=1

      countwf(1)=idim
      countwf(2)=1
      countwf(3)=1
      countwf(4)=1
c
      odz=1.0/dz
      ody=1.0/dy
      odx=1.0/dxx
      tdz=0.5/dz
      tdy=0.5/dy
      tdx=0.5/dx(1)
      fdz=1.0/(4.*dz)
      fdy=1.0/(4.*dy)
      fdx=1.0/(4.*dxx)
      tdzs=1.0/(dz*dz)
      tdys=1.0/(dy*dy)
      tdxs=1.0/(dxx*dxx)
      edz=0.25*tdz
      edy=0.25*tdy
      edx=0.25*tdx
c
      dtdx=delt/dxx
      dtdy=delt/dy
      dtdz=delt/dz
c
      maxu=dxx/(delt)
      maxv=dy/(delt)
      maxw=dz/(delt)
      maxt=40.
c



c
c set up rotation constants
c
      if(.not.rotate)then
         tomegaz=0.0
         tomegay=0.0 
      else if(.not.useomegy)then
         tomegay=0.0
      endif
c
      if(iope) then
        write (lunit,100) iver
 100    format (/a80)
        ilen=lnblk(ifile,80)
        write (lunit,110) ifile(1:ilen),itime
 110    format ('           File: ',a/
     1        '   File Created: ',a30)
        write (lunit,115)
 115    format (/'I.C. Comments: ')
        ilen=lnblk(com1,80)
        if (ilen.gt.0) write (lunit,120) com1(1:ilen)
 120    format (1x,a)
        ilen=lnblk(com2,80)
        if (ilen.gt.0) write (lunit,120) com2(1:ilen)
        ilen=lnblk(com3,80)
        if (ilen.gt.0) write (lunit,120) com3(1:ilen)
        if (ntrace.gt.0) then
           write (lunit,122)
 122       format (/'I.C. Tracer Info: ')
           do 128 i=1,ntrace
              write (dummy,'(a80)') trfile(i)
              ilen=lnblk(dummy,80)
              write (lunit,124) i,dummy(1:ilen)
 124          format (1x,'trcfile_',i1,': ',a)
              write (dummy,'(a120)') trccom(i)
              ilen=lnblk(dummy,120)
              write (lunit,126) i,dummy(1:ilen)
 126          format (2x,'trccom_',i1,': ',a)
 128       continue
        endif
        write (lunit,130)
 130    format()
c
      endif
c 
      nstart=int(start/delt+0.5)+1
      ns=nstart-1
c 
c boundary condition setup
c 
      ilen=lnblk(bcfile,80)
      write(*,*)"reading bc file"
      if(ilen.gt.0) then
         call bc_read(bcfile)
         bcon = .true.
      else
         if(iope) then
           write(lunit,'(a)')'B.C. file not entered, use defaults'
         endif
         bcon = .false.
         bctime=" "
         bcver=" "
      endif
c
c setup reference date
c
      datemode = (carddate .or. icdate) .and. usedate 
     *     .and. (bcdate .or. .not.bcon)
      if(datemode) then
         if(iope) then
           write(lunit,400)
 400       format(/' Date Mode is active')
         endif
         if(icdate .and. useicdate) then
            refdate(1) = icref(1)
            refdate(2) = icref(2)
         else
            refdate(1) = cardref(1)
            refdate(2) = cardref(2)
         endif
         call frmtdate(refdate, dummy)
         if(iope) then
           write(lunit,401)dummy
 401       format(4x, 'Reference Date = ',a20)
         endif
c
c may need to recompute start and ns!
c
      endif
c
#if defined(PRINTDEBUG)
      if(iope)write(*,*)"read in bcinit"
#endif
      call bcinit (ns,bc1,bc2,bc3,ier)
      if(.not.bcon) then
         bccom1=bc1
         bccom2=bc2
         bccom3=bc3
      endif
      if(ier.ne.0)then
        write(*,*) 'bcinit: error in B.C. initialization'
#if defined(USEMPI)
        call MPI_ABORT(comm,1,mpierr)
#endif
      endif
c
      if(iope)write(*,*)'call bcupdate'
      call bc_update(ns)
c
      if(iope)write(*,*)'call icinit'
      call icinit(ier)
      if(ier.ne.0)then 
        write(*,*) 'icinit: error in IC initialization'
#if defined(USEMPI)
        call MPI_ABORT(comm,1,mpierr)
#endif
      endif
c
      if(iope)write(*,*)'call oleminit'
      call olem_init(ier)
      if(ier.ne.0)then 
        write(*,*) 'sgsinit: error in SGS initialization'
#if defined(USEMPI)
        call MPI_ABORT(comm,1,mpierr)
#endif
      endif
c
c float initialization
c
      if(iope)write(*,*)'call flt routines'
      if(havefloats)then
         write(lunit,65)ffile
 65      format('Reading float setup form ',a)
c         call flt_read(ffile, ier)
         if(iope) then
           if(ier.ne.0)then
             write(*,*) 'floatinit: error in float initailization'
#if defined(USEMPI)
             call MPI_ABORT(comm,1,mpierr)
#endif
           endif
c           call flt_create(froot)
         endif
c         if(.not.isrestart)call flt_write(ns)
      endif
c 
c write output file header
c 
      stop=ntime*delt
c
#if defined(PRINTDEBUG)
      if(iope) write(*,*)"fld create"
#endif
      if(iope) then
        do 300 i=1,notimes
           call fld_create(ifile, bcfile, froot, i)
 300    continue
      endif
c
c Write the initial conditions
c
c      if(.not.isrestart)call fld_write (ns,nerr)
c
c setup diagnostic output
c
      if(iope)write(*,*)'call dag routines'
      if(ndiag.ne.0) then
         if(iope) then
           call dag_create(froot)
         endif
         if(iope)write(*,*)'dump dag file'
         if(.not.isrestart) call dag_write(ns,nerr)
      endif
c
c MJH MODIFY
c      if(iope) then
c        do i=1,slbtimes
c           call slb_create(ifile,bcfile,froot,i)
c        enddo
c      endif
c MJH MODIFY
c Changed to calculate s(X)loc for all processors under MPI
      do i=1,slbtimes
         if(siloc(i) .eq. -1) then
            index = int(sxloc(i)/dx(i) + 0.5)
            if(index .lt. 1) index = 1
            if(index .gt. idim) index = idim
            siloc(i) = index
         endif
         if(sjloc(i) .eq. -1) then
            index = int(syloc(i)/dy + 0.5)
            if(index .lt. 1) index = 1
            if(index .gt. jtotal) index = idim
            sjloc(i) = index
         endif
         if(skloc(i) .eq. -1) then
            index = kdim - int(szloc(i)/dz + 0.5)
            if(index .lt. 1) index = 1
            if(index .gt. kdim) index = kdim
            skloc(i) = index
         endif
         if(iope) then
            write(*,*)"create slb file"
            call slb_create(ifile,bcfile,froot,i)
         endif
      enddo
c END MJH MODIFY
      if(iope) write(*,*)"write slb file"
      if(.not.isrestart)call slb_write(ns,nerr)
c
c get cputime for initialization
c
      call deltime(cputime)
      if(iope) then
        write(lunit,62)cputime
 62    format(///'Model initialization required ',f10.2,'  CPU seconds')
      endif
c
      nerr=0
      veltim=0.0
      prestim=0.0
      advtim=0.0
      updtetim=0.0
      floattim=0.0
      writetim=0.0
c
c      if(iope) then
c #if defined(Win32)
c        close(lunit)
c        open(lunit,file=lfile,status='old',access='append')
c #else
c        call flush(lunit)
c #endif
c      endif
#if defined(USEMPI)
c
c let everything catch up
c
      call MPI_BARRIER(comm,mpierr)
#endif
c
c 
c**** start main time loop of model
c 
c
      do 190 n=nstart,ntime
         relstep = n - nstart + 1
c 
         if(iope) write(*,*) "Time Step ",n
         timstp=n
         call loop3d (n,nerr)
#if defined(USEMPI)
         if(nerr.eq.1) then
           call MPI_ABORT(comm,1,mpierr)
         endif
         call MPI_BARRIER(comm,mpierr)
#endif
c 
         if ((rstart .and. (mod(relstep,nrstrt).eq.0)) .or.
     *        sig_restart .or.
     *        ((n.eq.ntime).and.enddumpon))then
c            call energy(n)
            call restart(n,froot,ifile,bcfile)
#if defined(AVEFLUX)
c
c limit x axis to itotal/4.
c
c            if(iblkstart(my_pe+1).lt.itotal/3.) then
              call wrflux(n,froot,ifile,bcfile)
c            endif
#endif
            if(sig_quick)goto 200
            if(sig_usr1)then
               sig_usr1=.false.
#if defined(sun)
               siger = signal(SIGUSR1, sig_handler, -1)
#endif
            endif
            sig_restart=.false.
            if(sig_usr2)sig_quit=.true.
         endif
c
         if(nerr.eq.1) then
            call restart(n,froot,ifile,bcfile)
         endif
c 
         if (ndiag.ne.0 .and. mod(relstep,ndiag).eq.0)then
c            call energy(n)
            call dag_write(n,nerr)
         endif
c
c uses relstep internally to output step
c
         call slb_write(n,nerr)
c
         if (havefloats) then
           if(fltdt.ne.0 .and. mod(relstep,fltdt).eq.0)then
c              call flt_write(n)
           endif
         endif
c 
c         call fld_write (n,nerr)
c
         call deltime(cputime)
         writetim=writetim+cputime
c
         if (nerr.eq.1 .or. sig_quit) go to 200
c
 190  continue
 200  continue
c
c get cputime for model execution
c
      cputime=veltim+prestim+advtim+updtetim+floattim+writetim
      fctr=1.0e6/(float((n-nstart+1))*float(idim*kdim*jdim))
      if( cputime .gt. 0.0 ) then
         gpsec = float(n-nstart+1)*float(idim*jtotal*kdim)/cputime
      else
         gpsec = 0.0
      endif
      if(iope) then
         write(*,*)"TOTAL CPUTIME ",cputime
         if( cputime .gt. 0.0) then 
            write(lunit,63)veltim,veltim*100./cputime,
     *           prestim,prestim*100./cputime,advtim,
     *           advtim*100./cputime,
     *           updtetim,updtetim*100./cputime,floattim,
     *           floattim*100./cputime,writetim,writetim*100./cputime,
     *           cputime,cputime*fctr,gpsec
            write(*,63)veltim,veltim*100./cputime,
     *           prestim,prestim*100./cputime,advtim,
     *           advtim*100./cputime,
     *           updtetim,updtetim*100./cputime,floattim,
     *           floattim*100./cputime,writetim,writetim*100./cputime,
     *           cputime,cputime*fctr,gpsec
         else
            write(lunit,63)veltim, 0.0,
     *           prestim,0.0, advtim, 0.0, 
     *           updtetim, 0.0, floattim, 
     *           0.0, writetim, 0.0,
     *           cputime,cputime*fctr,gpsec
            write(*,63)veltim, 0.0,
     *           prestim,0.0, advtim, 0.0, 
     *           updtetim, 0.0, floattim, 
     *           0.0, writetim, 0.0,
     *           cputime,cputime*fctr,gpsec
         endif
 63     format(/'         Velocity calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '         Pressure calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '        Advection calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '           Update calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '            Float calculation ',f10.2,'  CPU seconds',
     *     f10.2,' percent'/
     *        '               Writing Output ',f10.2,'  CPU seconds',
     *     f10.2,' percent'//
     *        '     Model execution required ',f10.2,'  CPU seconds'//
     *        '      per time step per point ',f11.3,
     *     ' CPU microseconds'//
     *        '            Model performance ',f11.0,
     *        '  total grid points per second')
c
        close(lunit)
        write(*,*)' closing netcdf...'
c
        call ncclos(timid,rcode)
c
        do i=1,notimes
           call ncclos(fldid(i),rcode)
        enddo
c
        do i=1,slbtimes
           call ncclos(slbid(i),rcode)
        enddo
c
        if(havefloats)call ncclos(fltid, rcode)
      endif
c
#if defined(USEMPI)
      call MPI_FINALIZE(nerr)
#endif
c
#if defined(PRINTDEBUG)
      if(iope) then
         write(*,*)' MODEL COMPLETE.'
      endif
#endif
      end



