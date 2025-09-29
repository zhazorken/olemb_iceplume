      subroutine varcde (str,icode)
c
c   $Id: varcde.F,v 1.4 1998/07/27 20:12:36 dwd Exp $
c
      character*15 str
      integer icode,i,k
 
c     ivar  =  variable name (input)
c     icode =  variable code number (output)

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"
#include "varinit.inc"
 
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
c    FORTRAN_format, PPLUS_label, and missing_value.
c  Subroutines called: varunt, varstr, varfmt, which are in this file.
c  Functions called: lnblk from ~golem/src/include/lnblk.F
c  NetCDF subroutines are also called.
c-----------------------------------------------------------------------
      subroutine setvar(fld,cdfid,varid)

      integer fld,cdfid,varid
      integer ilen,lnblk,rcode
      character str*50,unt*20,fmt*10,ppl*80

#include "param.inc"
#include "netcdf.inc"
#include "moddef.inc"
#include "olemtke.inc"

      call varunt(fld,unt)
      ilen=max0(lnblk(unt,20),1)
      if(ilen.gt.0)call ncaptc(cdfid,varid,'units',NCCHAR,
     *     ilen,unt,rcode)

      if(fld .eq. MDTKE .and. .not.usedeardorff) then
         str = 'subgrid kinetic-energy dissipation'
      else
         call varstr(fld,str)
      endif

      ilen=max0(lnblk(str,50),1)
      if(ilen.gt.0)call ncaptc(cdfid,varid,'long_name',NCCHAR,
     *     ilen,str,rcode)

      call varfmt(fld,fmt)
      ilen=max0(lnblk(fmt,10),1)
      if(ilen.gt.0)call ncaptc(cdfid,varid,'FORTRAN_format',
     *     NCCHAR,ilen,fmt,rcode)

	  call varppl(fld,ppl)
      ilen=max0(lnblk(ppl,80),1)
      if(ilen.gt.0)call ncaptc(cdfid,varid,'PPLUS_label',
     *     NCCHAR,ilen,ppl,rcode)
      
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
c
	  subroutine varppl(icode,ppl)
      character*80 ppl
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"

	  do i=1,mdvlen
		 if(icode.eq.mdvcde(i)) then
		    ppl = mdvppl(i)
            return
         endif
      enddo
      ppl = ""
      return
      end
c-----------------------------------------------------------------------
      subroutine varunt(icode,unt)
      character*20 unt
      integer icode,i

#include "param.inc"
#include "moddef.inc"
#include "modvar.inc"

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
