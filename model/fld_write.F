      subroutine fld_write(n,nerr)
c
c   $Id: fld_write.F,v 1.9 1999/03/19 02:36:14 dwd Exp $
c
c 2.0:  Works with compiler switches for mpi or f77 version.

#if defined(USEMPI)
       include 'mpif.h'
#endif
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
#include "cyclic.inc"
#include "communicate.inc"
c 
#include "moddef.inc"
c 
      real tmp(0:idim1,0:jdim1,0:kdim1),zi,q,qs,beta,ws,ts
      real sigthta,kh,s,ustar,tau,rainfall
      integer time,strtt(4),countt(4)
      integer l,i,n,j,k,nerr,ii,kk,jj,rcode,jt
      logical first(5)
c
c temporary collected arrays
c
      real tmp2d(0:idim1,0:jtotal1)
c
      data first/.true.,.true.,.true.,.true.,.true./
      data countt/1,1,1,1/
      data strtt/1,1,1,1/
c
      j=1
C-----BEGIN MAIN LOOP
      do i=1,notimes
         if(((n.ge.nostart(i)).and.(n.le.nostop(i))
     *        .and.(mod(n-nostart(i),nodt(i)).eq.0)) .or.
     *        nerr.eq.1) then
c     
c     write output fields
c     
#if defined(PRINTDEBUG)
            write(*,*) "*****fld_write activated, n= ",n," proc= ",my_pe
#endif

            time=n*delt
            if(first(i))then
               first(i)=.false.
               tindexf(i) = 1
               strtf(4) = 1
               strtt(4) = 1
               if(iope) then
                  call ncvpt(fldid(i), tfid(i), 1, 1, time,rcode)
               endif
            else
               tindexf(i) = tindexf(i) + 1
               strtf(4)=tindexf(i)
               strtt(4)=tindexf(i)
               if(iope) then
                  call ncvpt(fldid(i), tfid(i), tindexf(i), 1, time,
     $                 rcode)
               endif
            endif

C-----BEGIN SECONDARY LOOP
            do jt=1,nofld(i)
               if (ofield(i,jt).eq.MDUVEL) then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   MDUVEL"
#endif
                  call write_fld_var( fldid(i), ufid(i), um, 
     *                 countuf, jtotal, kdim )
               else if (ofield(i,jt).eq.MDVVEL) then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   MDVVEL"
#endif
                  call write_fld_var( fldid(i), vfid(i), vm, 
     *                 countvf, jtotal1, kdim )
               else if (ofield(i,jt).eq.MDWVEL)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   MDWVEL"
#endif
                  call write_fld_var( fldid(i), wfid(i), wm, 
     *                 countwf, jtotal, kdim1 )
               else if (ofield(i,jt).eq.MDPRES)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   PRES"
#endif
                  call write_fld_var( fldid(i), pfid(i), p,
     *                 countf, jtotal, kdim )
               else if (ofield(i,jt).eq.MDTEMP)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   TEMP"
#endif
                  call write_fld_var( fldid(i), thfid(i), th,
     *                 countf, jtotal, kdim )
               else if (ofield(i,jt).eq.MDSAL)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   SAL"
#endif

#if !defined(SNGLSALT)
c MJH Do operation in parallel first: optimizes for clean code/parallel
                  do ii=1,idim
                     do j=1,jdim
                        do k=1,kdim
                           tmp(ii,j,k)=sal(ii,j,k)
                        enddo
                     enddo
                  enddo
                  call write_fld_var( fldid(i), salfid(i), tmp,
     *                 countf, jtotal, kdim )
#else
                  call write_fld_var( fldid(i), salfid(i), sal,
     *                 countf, jtotal, kdim )
#endif
               else if (ofield(i,jt).eq.MDTKE)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   TKE"
#endif
c MJH Do operation in parallel first: optimizes for clean code/parallel
                  do ii=1,idim
                     do j=1,jdim
                        do k=1,kdim
                           tmp(ii,j,k)=tke(ii,j,k)*tke(ii,j,k)
                        enddo
                     enddo
                  enddo
                  call write_fld_var( fldid(i), tkefid(i), tmp,
     *                 countf, jtotal, kdim )
               else if(ofield(i,jt).eq.MDKH)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   KH"
#endif
                  call write_fld_var( fldid(i), khfid(i), kheat,
     *                 countf, jtotal, kdim )
               else if(ofield(i,jt).eq.MDKM)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   KM"
#endif
                  call write_fld_var( fldid(i), kmfid(i), kmom,
     *                 countf, jtotal, kdim )
               else if (ofield(i,jt).eq.MDRHO)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   RHO"
#endif
                  call write_fld_var( fldid(i), rho_pfid(i), rhop,
     *                 countf, jtotal, kdim )
               else if (ofield(i,jt).eq.MDRHOT)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   RHOT"
#endif
c MJH Do operation in parallel first: optimizes for clean code/parallel
                  do ii=0,idim1
                     do jj=1,jdim
                        do kk=1,kdim
                           tmp(ii,jj,kk)=rhobm1(kk)+rhop(ii,jj,kk)
                        enddo
                     enddo
                  enddo
c
                  call write_fld_var( fldid(i), rhofid(i), tmp,
     *                 countf, jtotal, kdim )
               else if(ofield(i,jt).eq.MDHFSGS)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   HFSGS"
#endif

c MJH Do operation in parallel first: optimizes for clean code/parallel
                  do kk=1,kdim1
                     do jj=1,jdim
                        do ii=1,idim
                           kh=(kheat(ii,jj,kk)+kheat(ii,jj,kk-1))*0.5
                           tmp(ii,jj,kk)=-kh*(th(ii,jj,kk)-
     *                          th(ii,jj,kk-1))*odz
                        enddo
                     enddo
                  enddo
                  call write_fld_var( fldid(i), hf_sgsid(i), tmp,
     *                 countwf, jtotal, kdim1 )
#if defined(TRACER)
               else if (ofield(i,jt).ge.MDTRCR) then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   TRCR"
#endif
                  do l=1,ntrace
                     if(ofield(i,jt).eq.MDTRCR+l)then
                        call write_fld_var( fldid(i), 
     *                       rhofid(i), trc(0,0,0,1),
     *                       countf, jtotal, kdim )
                     end if
                  enddo
#endif
               else
                  write (lunit,*) ' field type ',ofield(i,jt),
     *                 ' not recognized'
               end if
C-----END SECONDARY LOOP - jt=1,nofld(i)
            enddo
c
c allways output the following!
c
c
c zi
c
            call zmixed(zi)
            if(iope) then 
               call ncvpt(fldid(i),zifid(i),strtt,countt,zi,rcode)
            endif
c     
c     compute w_star, u_star, t_star, q
c     
            q=0.0

#if defined(USEMPI)
            call collect_2d( randn, tmp2d )
            do jj=1,jtotal
               do ii=1,idim
                  q=q+heatf*tmp2d(ii,jj)
               enddo
            enddo
#else
            do jj=1,jdim
               do ii=1,idim
                  q=q+heatf*randn(ii,jj)
               enddo
            enddo
#endif
            q=q/float(idim*jtotal)
c     
            qs=-q/(cp*rho0)
c     
            s=sal(1,1,kdim)
            beta=-(sigthta(s,th(1,1,kdim)+0.001,0.0) - 
     *           sigthta(s,th(1,1,kdim),0.0))/(0.001*rho0)
c     
            If(beta.le.0 .or. qs.le.0)then
               ws=1.
               ts=1.
            else
               ws=(beta*grav*zi*qs)**0.3333333
               ts=qs/ws
            endif
            tau = sqrt(ustr**2 + vstr**2)
            ustar = sqrt(tau/rho0)
            if(rainon)rainfall = 3600.*rainf*1000.0
            if(iope) then 
c     
c     w_star
c     
               call ncvpt(fldid(i),w_starfid(i),strtt,countt,ws,
     *              rcode)
c     
c     u_star
c     
               call ncvpt(fldid(i),u_starfid(i),strtt,countt,ustar,
     *              rcode)
c     
c     t_star
c     
               call ncvpt(fldid(i),t_starfid(i),strtt,countt,ts,
     *              rcode)
c     
c     q
c     
               call ncvpt(fldid(i),qfid(i),strtt,countt,q,rcode)
c     
c     hf_top
c     
               if(heaton) call ncvpt(fldid(i),hf_topfid(i),strtt,
     *              countt,heatf,rcode)
c     
c     ustr_t
c     
               if(uston) call ncvpt(fldid(i),ustr_tfid(i),strtt,
     *              countt,ustr,rcode)
c     
c     vstr_t
c     
               if(vston) call ncvpt(fldid(i),vstr_tfid(i),strtt,
     *              countt,vstr,rcode)
c     
c     swf_top
c     
               if(swheaton) call ncvpt(fldid(i),swf_topfid(i),strtt,
     *              countt,swheatf,rcode)
c     
c     rain
c     
               if(rainon) call ncvpt(fldid(i),rainfid(i),strtt,
     *              countt,rainfall,rcode)
c     
c     lhf_top
c     
               if(lheaton) call ncvpt(fldid(i),lhf_topfid(i),
     *              strtt,countt,latent,rcode)
c     
#if defined(STOKES)
c     
c     wave_l
c     
               if(wavelon) call ncvpt(fldid(i),wave_lfid(i),strtt,
     *              countt,wavelen,rcode)
c     
c     wave_h
c     
               if(wavehon) call ncvpt(fldid(i),wave_hfid(i),strtt,
     *              countt,waveamp,rcode)
c     
c     w_angle
c     
               if(angleon) call ncvpt(fldid(i),w_anglefid(i),strtt,
     *              countt,waveang,rcode)
#endif
c     Endif(iope)
            endif
c     
c     sync the netCDF file
c     
            if(iope) then
               call ncsnc(fldid(i), rcode)
            endif
         endif
         if(nerr.eq.1)return
C-----END MAIN LOOP - i=1,notimes
      enddo
      return
      end
      
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C-----Output subroutine for field variables------------------------
c
c
      subroutine write_fld_var( ncid, varid, vararr, loc_count, 
     *     loc_jmax, loc_kmax )
c
c
#include "param.inc"
#include "olemcdf.inc"
#include "communicate.inc"
      integer ncid, varid
      real vararr(0:idim1,0:jdim1,0:kdim1)
      integer loc_count(4), loc_jmax, loc_kmax
      integer j,k
      integer rcode
#if defined(USEMPI)
      real tmp(0:idim1,0:jtotal1,0:kdim1)

      call collect_3d( vararr, tmp )
#endif
      if(iope) then 
         do j=1,loc_jmax
            strtf(2)=j
            do k=1,loc_kmax
               strtf(3)=k
#if defined(USEMPI)
               call ncvpt(ncid,varid,strtf,loc_count,
     *              tmp(1,j,k),rcode)
#else
               call ncvpt(ncid,varid,strtf,loc_count,
     *              vararr(1,j,k),rcode)
#endif
            enddo
         enddo
      endif
      return
      end




