	subroutine mks(str,val,flgs,unts)
c
c   $Id: mks.F,v 1.2 1998/08/04 18:12:50 dwd Exp $
c
	character str*(*)
	real val,fact
	logical flgs(3)
	integer unts(3),i,nchar,iex,isign
c
c returns a value in mks units
c    where:
c            M**im K**ik S**is
c
c    syntax:
c		str := fpn(units)
c		fpn := floating point number
c	       units:= length, mass or time
c
c    example:
c		4.3 (m)   (4.3 meters)
c		1.7e2 (km) (1700 kilometers)
c		1000 (kgm-3) (1000 kilograms per cubic meter)
c
c    M   is type 1
c    K   is type 2
c    S   is type 3
c
c    built in definitions:
c
c        N = kg m s-2
c        W = kg m2 s-3
c
	character sstr*80,ustr*80,frmt*30,tmp*80
	integer iparno,iparnc,slen,ulen,ispc,ilen
c
	character unit(9)*3
	integer lunit(9),tunit(9),nunits,type
	real funit(9)
c
        character def(2)*1,result(2)*20
        integer   ldef(2),lresult(2),ndefs
c
	data  unit/'min','m', 'km','cm','kg',  'g','s',   'd',  'h'/
	data lunit/    3,  1,    2,   2,   2,    1,  1,     1,    1/
	data tunit/    3,  1,    1,   1,   2,    2,  3,     3,    3/
	data funit/  60.,1.0,1000.,0.01, 1.0,0.001,1.0,86400.,3600./
	data nunits/9/
c
        data     def/       'N',        'W'/
        data    ldef/         1,          1/
        data  result/'kg m s-2','kg m2 s-3'/
        data lresult/         8,          9/
        data ndefs/2/
c
	flgs(1)=.false.
	flgs(2)=.false.
	flgs(3)=.false.
c
	unts(1)=0
	unts(2)=0
	unts(3)=0
	ilen=len(str)
c
	iparno=index(str,'(')
        if(iparno.eq.0)iparno=ilen
	iparnc=index(str,')')
	ispc=index(str,' ')
	if(ispc.ne.0 .and. ispc.lt.iparno)then
	    slen=ispc-1
	else
	    slen=iparno-1
	endif
c
	if(iparno.eq.ilen)then
	    write(frmt,999)slen
999	    format('(f',i3.3,'.0)')
	    read(str(:slen),frmt)val
	    return
	else
	    sstr=str(1:slen)
	    ustr=str(iparno+1:iparnc-1)
            ulen=iparnc-1-iparno
	    call squish(ustr,1,ulen)
	    write(frmt,999)slen
	    read(sstr,frmt)val
	endif
c
c find a substitute definitions
c
 20     do 30 i=1,ndefs
           if(ustr(1:ldef(i)).eq.def(i)(:ldef(i))) then
              ustr=ustr(ldef(i)+1:ulen)//result(i)(:lresult(i))
              ulen=ulen-ldef(i)+lresult(i)
              call squish(ustr,1,ulen)
           endif
 30     continue
c
c determine units
c
	do 10 i=1,nunits
	    if(ustr(1:lunit(i)).eq.unit(i)(:lunit(i)))goto 11
10	continue
c
c error unit not recognized
c
	return
c
11	fact=funit(i)
	type=tunit(i)
	nchar=lunit(i)+1
c
c get sign of exponent
c
	isign=1
	if(ustr(nchar:nchar).eq.'-')then
	    isign=-1
	    nchar=nchar+1
	else if(ustr(nchar:nchar).eq.'+')then
	    isign=1
	    nchar=nchar+1
	endif
c
c get exponent
c
	iex=1
       	if(lge(ustr(nchar:nchar),'0') .and. 
     *		lle(ustr(nchar:nchar),'9'))then
	    read(ustr(nchar:nchar),998)iex
998	    format(i1)
	    nchar=nchar+1
	endif
c
	fact=fact**(isign*iex)
        unts(type)=unts(type)+iex*isign
        flgs(type)=unts(type).ne.0
c
	ulen=ulen-nchar+1
	tmp=ustr(nchar:)
	ustr=tmp
	val=val*fact
	if(ulen.gt.0)goto 20
	return
	end
c
      SUBROUTINE SQUISH(JSTRNG,JFIRST,NPOS)
	integer jfirst,npos
C** 
C**    @(#)squish.f	1.2    9/18/89
C**
C**
C***********************************************************************
C**
C**                 PLOT+ Scientific Graphics System
C**
C***********************************************************************
C**
C**
C
C	Removes blank characters starting at position JFIRST
C	for the next NPOS characters and packs the characters
C	back into the original string.  Good for getting rid
C	of extra blanks.
C
C	JSTRNG -- character, input and output string
C	JFIRST -- INT, position to begin processing
C	NPOS -- INT, number of characters to process
C
      CHARACTER*1 JSTRNG*120,KLPAR,KRPAR,JBLNK,JCH
      LOGICAL REMOVE
	integer next,nparen,jpos,last

      DATA KLPAR/'('/,KRPAR/')'/,JBLNK/' '/

      IF(NPOS.EQ.0)RETURN
      LAST=IABS(NPOS)+JFIRST-1
      REMOVE=.FALSE.
      IF(NPOS.LT.0)REMOVE=.TRUE.
      NEXT=JFIRST
      NPAREN=0

      DO 500 JPOS=JFIRST,LAST
	JCH=JSTRNG(JPOS:JPOS)
	JSTRNG(JPOS:JPOS)=JBLNK
      IF(JCH.NE.KLPAR)GO TO 100
      NPAREN=NPAREN+1
      GO TO 150

  100 IF(NPAREN.EQ.0)GO TO 200
      IF(JCH.NE.KRPAR)GO TO 300
      NPAREN=NPAREN-1
  150 IF(REMOVE)JCH=JBLNK
      GO TO 300

  200 IF(JCH.EQ.JBLNK)GO TO 500
  300 JSTRNG(NEXT:NEXT)=JCH
      NEXT=NEXT+1
  500 CONTINUE

      RETURN
      END
c
	logical function mkstst(flgs,unts,fm,fk,fs,im,ik,is)
	logical flgs(3),fm,fk,fs
	integer unts(3),im,ik,is
c
c returns .true. if there is an error
c
c returns .true. if flgs and fm,fk and fs are different
c                   unless flgs are all false
c
	logical tstm,tstk,tsts
c
	tstm=flgs(1).eqv.fm
	tstk=flgs(2).eqv.fk
	tsts=flgs(3).eqv.fs
c
	if((.not.flgs(1)).and.(.not.flgs(2)).and.(.not.flgs(3)))then
	    mkstst=.false.
	    return
	endif
c
	if(tstm.and.tstk.and.tsts)then
	    mkstst=.false.
	    if(tstm .and. im.ne.unts(1))mkstst=.true.
	    if(tstk .and. ik.ne.unts(2))mkstst=.true.
	    if(tsts .and. is.ne.unts(3))mkstst=.true.
	else
	    mkstst=.true.
	endif
	return
	end
