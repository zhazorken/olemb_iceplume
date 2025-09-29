      subroutine slb_write (n,nerr)
c
c   $Id: slb_write.F,v 1.12 1999/08/16 20:22:32 eds Exp $
c
#include "param.inc"
#include "olemf.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "olemcdf.inc"
#include "netcdf.inc"
#include "cyclic.inc"
#include "finitv.inc"
#include "communicate.inc"
c #if defined(ICE)
#include "ice.inc"
c #endif
c 
#include "moddef.inc"
c 
      real temp(0:idim1,0:jdim1,0:kdim1)
      real fmeltl(jdim,kdim),fmeltc(jdim,kdim)
c
      integer time,strtt(4)
      integer countx(4),county(4),countxw(4),countyw(4)
      integer countz(4)
      integer i,n,j,k,nerr,rcode,jt,l,si,sj,sk,cdfout,type,ko
      logical first
c
      data first/.true./
      data strtt/1,1,1,1/
c
      j=1
      countx(1) = itotal
      countx(2) = 1
      countx(3) = kdim
      countx(4) = 1
c
      county(1) = 1
      county(2) = jtotal
      county(3) = kdim
      county(4) = 1
c
      countxw(1) = itotal 
      countxw(2) = 1  
      countxw(3) = kdim1
      countxw(4) = 1  
c
      countyw(1) = 1
      countyw(2) = jtotal
      countyw(3) = kdim1
      countyw(4) = 1
c
      countz(1) = itotal
      countz(2) = jtotal
      countz(3) = 1
      countz(4) = 1
      if(ndiag.ne.0) then
      do k=1,kdim
        do j=1,jdim
           do i=1,idim
c             dcor(i,j,k)=dcor(i,j,k)/float(ndiag)
c             dpres(i,j,k)=dpres(i,j,k)/float(ndiag)
c             dvadvv(i,j,k)=dvadvv(i,j,k)/float(ndiag)
           enddo
        enddo
      enddo
      endif
c
c      call MPI_ALLREDUCE(fmelt,fmeltc,jdim*kdim,
c     *                      MPI_REAL,MPI_SUM,comm,nerr)


c 
c write output fields
c 
      do ko=1,slbtimes
c
         cdfout=slbid(ko)
         type=slbtype(ko)
         si=siloc(ko)
         sj=sjloc(ko)
         sk=skloc(ko)
c
         if(nslab(ko).ne.0 .and. mod(relstep,nslab(ko)).eq.0)then
#if defined(PRINTDEBUG)
            write(*,*) "*****slb_write activated, n= ",n," proc= ",my_pe
#endif
            time=n*delt
            if(iope) then
              if(first)then
                 first=.false.
                 tindexs(ko) = 1
                 strtt(4) = 1
                 call ncvpt(cdfout, tsid(ko), 1, 1, time,rcode)
              else
                 tindexs(ko) = tindexs(ko) + 1
                 strtt(4)=tindexs(ko)
                 call ncvpt(cdfout, tsid(ko), tindexs(ko), 1, time,
     *                rcode)
              endif
            endif
c
            do jt=1,noslb(ko)
c
c first output the surface ice slab
c
#if defined(ICE)
               if (oslab(ko,jt).eq.MDHICE) then
                 call slab2d(cdfout,hice,hifidz(ko),strtt,countz)
               else if(oslab(ko,jt).eq.MDFICE) then
                 call slab2d(cdfout,iceh,hfifidz(ko),strtt,countz)
               else if(oslab(ko,jt).eq.MDSFICE) then
                 call slab2d(cdfout,icesal,sfifidz(ko),strtt,countz)
c
               else if (oslab(ko,jt).eq.MDUVEL) then
#else
               if (oslab(ko,jt).eq.MDUVEL) then
#endif
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:MDUVEL ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                     call vertslab(cdfout,umm,ufidx(ko),ufidy(ko),strtt,
     *                 countx,county,si,sj) 
                  endif
                  call horzslab(cdfout,umm,ufidz(ko),strtt,countz,sk)

#if defined(HBOT)
               else if (oslab(ko,jt).eq.MDHBOT) then
                 call slab2d(cdfout,umass,botfidz(ko),strtt,countz)
                 do i=0,idim1
                   do j=0,jdim1
                     umass = 0.0
                   enddo
                 enddo
#endif
               else if (oslab(ko,jt).eq.MDUF) then
                  call horzslab(cdfout,uf,uffidz(ko),strtt,countz,sk)
               else if (oslab(ko,jt).eq.MDVF) then
                  call horzslab(cdfout,vf,vffidz(ko),strtt,countz,sk)
               else if (oslab(ko,jt).eq.MDWF) then
                  call horzslab(cdfout,wf,wffidz(ko),strtt,countz,sk)
               else if (oslab(ko,jt).eq.MDRHOF) then
                  call horzslab(cdfout,pf,rffidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDVVEL) then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:MDVVEL  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                     call vertslab(cdfout,vm,vfidx(ko),vfidy(ko),strtt,
     *                 countx,county,si,sj)
                  endif
                  call horzslab(cdfout,vm,vfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDWVEL)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:MDWVEL  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                    call vertslabw(cdfout,wm,wfidx(ko),wfidy(ko),strtt,
     *                 countxw,countyw,si,sj)
                  endif
                  call horzslab(cdfout,wm,wfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDVADV)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:vadv  ", my_pe
#endif
c                 do k=1,kdim
c                    do j=1,jdim
c                      do i=1,idim
c                        dvadvv(i,j,k)=dvadvv(i,j,k)/float(nslab(ko))
c                      enddo
c                    enddo
c                  enddo

                 if(type.eq.SLAB_ALL)then
                    call vertslab_a(cdfout,dvadvv,vafidx(ko),vafidy(ko),
     *                 strtt,countx,county,si,sj)
                 endif
                 call horzslab(cdfout,dvadvv,vafidz(ko),strtt,countz,sk)
c
               else if (oslab(ko,jt).eq.MDDCOR)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:dcor  ", my_pe
#endif
c                  do k=1,kdim
c                    do j=1,jdim
c                      do i=1,idim
c                        dcor(i,j,k)=dcor(i,j,k)/float(nslab(ko))
c                      enddo
c                    enddo
c                  enddo
                  if(type.eq.SLAB_ALL)then
                    call vertslab_a(cdfout,dcor,dcfidx(ko),dcfidy(ko),
     *                 strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,dcor,dcfidz(ko),strtt,countz,sk)
c
               else if (oslab(ko,jt).eq.MDDPRES)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:dpres  ", my_pe
#endif
c                 do k=1,kdim
c                    do j=1,jdim
c                      do i=1,idim
c                        dpres(i,j,k)=dpres(i,j,k)/float(nslab(ko))
c                      enddo
c                    enddo
c                  enddo

                  if(type.eq.SLAB_ALL)then
                    call vertslab_a(cdfout,dpres,dpfidx(ko),dpfidy(ko),
     *                 strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,dpres,dpfidz(ko),strtt,countz,sk)
c
               else if (oslab(ko,jt).eq.MDSGF)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:sgf  ", my_pe
#endif
c
c using to dump out vertical slab average
c
                  if(type.eq.SLAB_ALL)then
                    call vertslab_a(cdfout,sgf,sgfidx(ko),sgfidy(ko),
     *                 strtt,countx,county,si,sj)
                  endif
c                  call horzslaba(cdfout,sgf,sgfidz(ko),strtt,countz,sk)
                  call horzslaba(cdfout,sgf,sgfidz(ko),strtt,countz,sk)
c

               else if (oslab(ko,jt).eq.MDPRES)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:PRES  ", my_pe
#endif
                do k=1,kdim
                    do j=1,jdim
                      do i=1,idim
                        temp(i,j,k)=p(i,j,k)/rhob(k)
                      enddo
                    enddo
                  enddo

                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,temp,pfidx(ko),pfidy(ko),strtt,
     *                 countx,county,si,sj)
                  endif
                  call horzslab(cdfout,temp,pfidz(ko),strtt,countz,sk)
               else if (oslab(ko,jt).eq.MDTEMP)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:TEMP  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,th,thfidx(ko),thfidy(ko),strtt,
     *                 countx,county,si,sj)
                  endif
                  call horzslab(cdfout,th,thfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDSAL)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:SAL  ", my_pe
#endif
                  do k=0,kdim1
                    do j=0,jdim1
                      do i=0,idim1
                        temp(i,j,k) = sal(i,j,k)
                      enddo
                    enddo
                  enddo

                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,temp,salfidx(ko),salfidy(ko),
     *                   strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,temp,salfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDTKE)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:TKE  ", my_pe
#endif

c                  do i=1,idim
c                  enddo
                  if(type.eq.SLAB_ALL)then
                   call vertslab_a(cdfout,tke,tkefidx(ko),tkefidy(ko),
     *                   strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,tke,tkefidz(ko),strtt,countz,sk)
c                  do k=0,kdim1
c                    do j=0,jdim1
c                      do i=0,idim1
c                        tke(i,j,k) = 0.0
c                      enddo
c                    enddo
c                  enddo
               else if (oslab(ko,jt).eq.MDRHO)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:RHO  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,rhop,rho_pfidx(ko),
     *                   rho_pfidy(ko),strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,rhop,rho_pfidz(ko),strtt,
     *                 countz,sk)

               else if (oslab(ko,jt).eq.MDKH)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:KH  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,kheat,khfidx(ko),khfidy(ko),
     *                   strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,kheat,khfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDKM)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:KM  ", my_pe
#endif
                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,kmom,kmfidx(ko),kmfidy(ko),
     *                   strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,kmom,kmfidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDRHOT)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:RHOT  ", my_pe
#endif
                     do k=0,kdim1
                       do j=0,jdim1
                         do i=0,idim1
                           temp(i,j,k) = rhop(i,j,k)+rhobm1(k)
                         enddo
                       enddo
                     enddo

                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,temp,rhofidx(ko),rhofidy(ko),
     *                   strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,temp,rhofidz(ko),strtt,countz,sk)

               else if (oslab(ko,jt).eq.MDHFSGS)then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:HFSGS  ", my_pe
#endif
                  do k=1,kdim
                    do i=1,idim
                      do j=1,jdim
                        temp(i,j,k)=(kheat(i,j,k)+kheat(i,j,k-1))*
     *                       0.5*(th(i,j,k)-th(i,j,k-1))*odz
                      enddo
                    enddo
                  enddo

                  if(type.eq.SLAB_ALL)then
                    call vertslab(cdfout,bax,hf_sgsfidx(ko),
     *                   hf_sgsfidy(ko),strtt,countx,county,si,sj)
                  endif
                  call horzslab(cdfout,bax,hf_sgsfidz(ko),strtt,
     *                 countz,sk)
#if defined(TRACER)
               else if (oslab(ko,jt).ge.MDTRCR) then
#if defined(PRINTDEBUG)
                  write(*,*) "*****   slb:TRCR  ", my_pe
#endif
                  do l=1,ntrace
                     do k=0,kdim1
                       do j=0,jdim1
                         do i=0,idim1
                           temp(i,j,k) = trc(i,j,k,l)
                         enddo
                       enddo
                     enddo
                     if(oslab(ko,jt).eq.MDTRCR)then
                        if(type.eq.SLAB_ALL)then
                           call vertslab(cdfout,temp,trcrfidx(ko,l),
     *                           trcrfidy(ko,l),strtt,countx,county,
     *                           si,sj)
                        endif
                        call horzslab(cdfout,temp,trcrfidz(ko,l),strtt,
     *                       countz,sk)
                     endif
                  enddo
#endif
               else
                  if(iope) then
                    write (lunit,*) ' field type ',oslab(ko,jt),
     *                 ' not recognized'
                  endif
               end if
c
            enddo
c
c sync the netCDF file
c
            if(iope) then
              call ncsnc(cdfout, rcode)
            endif
            do k=1,kdim
              do j=1,jdim
                 do i=1,idim
c                   dcor(i,j,k)= 0.0
c                   dpres(i,j,k)= 0.0
c                   dvadvv(i,j,k)= 0.0
                 enddo
              enddo
            enddo

#if defined(PRINTDEBUG)
            write(*,*) "*****slb_write completed, n= ",n," proc= ",my_pe
#endif
         endif
         if(nerr.eq.1)return
      enddo
      return
      end
c
c routines for writing out vertical slabs and horizontal slabs
c
c
      subroutine vertslab(cdf,phi,fidx,fidy,strtt,countx,county,si,sj)
c
#include "param.inc"
      integer cdf,fidx,fidy,strtt(4),countx(4),county(4),si,sj
      integer rcode,i,j,k
      real phi(0:idim1,0:jdim1,0:kdim1)
      real tmpxz(itotal,kdim1),tmpyz(jtotal,kdim1)
#include "communicate.inc"
c
#if defined(USEMPI)
c      call collect_yz_ave(phi,tmpyz,si)
      call collect_yz(phi,tmpyz,si)
      call collect_xz(phi,tmpxz,sj)
c      call collect_xz_ave(phi,tmpxz,sj)
#else
      do k=1,kdim
        do i=1,idim
          tmpxz(i,k)=phi(i,sj,k)
        enddo
        do j=1,jdim
          tmpyz(j,k)=phi(si,j,k)
        enddo
      enddo
#endif
      if(iope) then
        call ncvpt(cdf,fidx,strtt,countx,tmpxz,
     *                    rcode)
        call ncvpt(cdf,fidy,strtt,county,tmpyz,
     *                    rcode)
      endif
      return
      end
c
c routines for writing out vertical slabs and horizontal slabs
c
c
      subroutine vertslab_a(cdf,phi,fidx,fidy,strtt,countx,county,si,sj)
c
#include "param.inc"
      integer cdf,fidx,fidy,strtt(4),countx(4),county(4),si,sj
      integer rcode,i,j,k
      real phi(0:idim1,0:jdim1,0:kdim1)
      real tmpxz(itotal,kdim1),tmpyz(jtotal,kdim1)
#include "communicate.inc"
c
#if defined(USEMPI)
      call collect_yz_ave(phi,tmpyz,si)
c      call collect_yz(phi,tmpyz,si)
c      call collect_xz(phi,tmpxz,sj)
      call collect_xz_ave(phi,tmpxz,sj)
#else
      do k=1,kdim
        do i=1,idim
          tmpxz(i,k)=phi(i,sj,k)
        enddo
        do j=1,jdim
          tmpyz(j,k)=phi(si,j,k)
        enddo
      enddo
#endif
      if(iope) then
        call ncvpt(cdf,fidx,strtt,countx,tmpxz,
     *                    rcode)
        call ncvpt(cdf,fidy,strtt,county,tmpyz,
     *                    rcode)
      endif
      return
      end

c
      subroutine vertslabw(cdf,phi,fidx,fidy,strtt,countx,county,si,sj)
c
#include "param.inc"
      integer cdf,fidx,fidy,strtt(4),countx(4),county(4),si,sj
      integer rcode,i,j,k
      real phi(0:idim1,0:jdim1,0:kdim1)
      real tmpxz(itotal,kdim1),tmpyz(jtotal,kdim1)
#include "communicate.inc"
c
#if defined(USEMPI)
      call collect_yz(phi,tmpyz,si)
c      call collect_xz_ave(phi,tmpxz,sj)
      call collect_xz(phi,tmpxz,sj)
#else
      do k=1,kdim
        do i=1,idim
          tmpxz(i,k)=phi(i,sj,k)
        enddo
        do j=1,jdim
          tmpyz(j,k)=phi(si,j,k)
        enddo
      enddo
#endif
      if(iope) then

        do j=1,jtotal
          tmpyz(j,kdim1) = 0.0
        enddo
        do i=1,itotal
          tmpxz(i,kdim1) = 0.0
        enddo

        call ncvpt(cdf,fidx,strtt,countx,tmpxz,
     *                    rcode)
        call ncvpt(cdf,fidy,strtt,county,tmpyz,
     *                    rcode)
      endif
      return

      end
c
c
      subroutine horzslab(cdf,phi,fid,strtt,countz,sk)
c
#if defined(USEMPI)
      include "mpif.h"
#endif
c
#include "param.inc"

      integer cdf,fid,strtt(4),countz(4),sk,rcode,i,j
      real phi(0:idim1,0:jdim1,0:kdim1),tmpxy(itotal,jtotal)
#include "communicate.inc"
c
#if defined(USEMPI)
      call collect_xy(phi,tmpxy,sk)
#else
      do j=1,jdim
        do i=1,idim
          tmpxy(i,j)=phi(i,j,sk)
        enddo
      enddo
#endif
      if(iope) then
        call ncvpt(cdf,fid,strtt,countz,tmpxy,rcode)
      endif
c
      return
      end
      subroutine horzslaba(cdf,phi,fid,strtt,countz,sk)
c
#if defined(USEMPI)
      include "mpif.h"
#endif
c
#include "param.inc"

      integer cdf,fid,strtt(4),countz(4),sk,rcode,i,j,ktem
      real phi(0:idim1,0:jdim1,0:kdim1),tmpxy(itotal,jtotal)
      real phit(0:idim1,0:jdim1,0:kdim1)

#include "communicate.inc"
c
#if defined(USEMPI)
      do j=1,jdim
        do i=1,idim
          phit(i,j,sk) = 0.0
          do k=kdim-61,kdim
          if(phi(i,j,k).gt.0.0) then
            phit(i,j,sk)=phi(i,j,k)+phit(i,j,sk)
          endif
          enddo
        enddo
      enddo

      call collect_xy(phit,tmpxy,sk)
#else
c
c compute vertical average
c
#endif
      if(iope) then
        call ncvpt(cdf,fid,strtt,countz,tmpxy,rcode)
      endif
c
      return
      end


      subroutine slab2d(cdf,phi,fid,strtt,countz)
c
#if defined(USEMPI)
      include "mpif.h"
#endif
c
#include "param.inc"

      integer cdf,fid,strtt(4),countz(4),sk,rcode,i,j
      real phi(0:idim1,0:jdim1),tmpxy(itotal,jtotal)
      real tmpout(idim,jtotal)
#include "communicate.inc"
c
#if defined(USEMPI)
      call collect_2d(phi,tmpxy)
      do j=1,jtotal
        do i=1,idim
          tmpout(i,j)=tmpxy(i,j)
        enddo
      enddo
#else
      do j=1,jdim
        do i=1,idim
          tmpout(i,j)=phi(i,j)
        enddo
      enddo
#endif
      if(iope) then
        call ncvpt(cdf,fid,strtt,countz,tmpout,rcode)
      endif

      return
      end
