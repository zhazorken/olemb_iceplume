      integer function lnblk (str,ilen)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
      character str*(*)
      integer ilen,i
c 
      do 10 i=ilen,1,-1
          if (str(i:i).ne.' ') go to 20
 10   continue
      i=0
 20   lnblk=i
      return
      end

c-----------------------------------------------------------------------
      subroutine varcde (str,icode)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
      character*15 str
      integer icode,i,k
 
c     ivar  =  variable name (input)
c     icode =  variable code number (output)

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
#include "varinit.inc"
       integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(str.eq.mdvnme(i))goto 11
         if(str(1:4).eq.'trcr')goto 12
 10   continue
      icode=-1
      return
 11   icode=mdvcde(i)
      return
 12   read(str(5:5),'(I1)')k
      icode=MDTRCR+k
      return
      end

c-----------------------------------------------------------------------
      subroutine varnme(icode,str)
      character*15 str
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if(icode.ge.MDTRCR)then
         write(str,999)icode-MDTRCR
 999     format('trcr',i1)
      else
         str=' '
      endif
      return
 11   str=mdvnme(i)
      return
      end

c-----------------------------------------------------------------------
c  Assigns 4 attributes to the specified variable:  units, long_name,
c    FORTRAN_format, and missing_value.
c  Called by winit3, a subroutine of ~golem/src/bin/winit3.F,
c    and other programs...
c  Subroutines called: varunt, varstr, varfmt, which are in this file.
c  Functions called: lnblk from ~golem/src/include/lnblk.F
c  NetCDF subroutines are also called.
c-----------------------------------------------------------------------
      subroutine setvar(fld,cdfid,varid)

      integer fld,cdfid,varid
      integer ilen,lnblk,rcode
      character str*50,unt*20,fmt*10

#include "netcdf.inc"

      call varunt(fld,unt)
      ilen=max0(lnblk(unt,20),1)
      call ncaptc(cdfid,varid,'units',NCCHAR,ilen,unt,rcode)

      call varstr(fld,str)
      ilen=max0(lnblk(str,50),1)
      call ncaptc(cdfid,varid,'long_name',NCCHAR,ilen,str,rcode)

      call varfmt(fld,fmt)
      ilen=max0(lnblk(fmt,10),1)
      call ncaptc(cdfid,varid,'FORTRAN_format',NCCHAR,ilen,fmt,rcode)

      call ncapt(cdfid,varid,'missing_value',NCFLOAT,1,1.e35,rcode)

      return
      end

c-----------------------------------------------------------------------
      subroutine varstr(icode,str)
      character*50 str
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if(icode.ge.MDTRCR)then
         write(str,999)icode-MDTRCR
 999     format('Tracer ',i1)
      else
         str=' '
      endif
      return
 11   str=mdvstr(i)
      return
      end

c-----------------------------------------------------------------------
      subroutine varunt(icode,unt)
      character*20 unt
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      unt=' '
      return
 11   unt=mdvunt(i)
      return
      end

c-----------------------------------------------------------------------
c  modified for tracers by jrb on 7-20-92
c-----------------------------------------------------------------------
      subroutine varfmt(icode,fmt)
      character*10 fmt
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
      integer psize
      parameter (psize=300)

      do 10 i=1,mdvlen
         if(icode.eq.mdvcde(i))goto 11
 10   continue
      if (icode.ge.MDTRCR) then
         fmt='t6.2'
      else
         fmt='e10.2'
      endif
      return
 11   fmt=mdvfmt(i)
      return
      end

      real function ran2 (idum)
c
c   $Id: ic_create.F,v 4.8 1998/10/02 15:44:09 eds Exp $
c
c        returns a uniform random deviate between 0.0 and 1.0.  Set
c        IDUM to any negative value to initialize or reinitialize
c        the sequence.
c
      save
      integer m,ia,ic
      real rm 
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      integer ir(97),idum,iff,j,mod,iy
      data iff/0/
      if (idum.lt.0.or.iff.eq.0) then
          iff=1
          idum=mod(ic-idum,m)
          do 10 j=1,97
              idum=mod(ia*idum+ic,m)
              ir(j)=idum
 10       continue
          idum=mod(ia*idum+ic,m)
          iy=idum
      end if
      j=1+(97*iy)/m
      if (j.gt.97.or.j.lt.1)then
         write(*,*)'Problem with random number'
         return
      endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
