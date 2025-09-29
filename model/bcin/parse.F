	subroutine parse(str,value)
	character str*(*),value*(*)
c
c   $Id: parse.F,v 4.2 1998/06/18 14:08:01 dwd Exp $
c
c get next value from string  value is terminated by a comma
c
	character tmp*160
	integer icma,ilen,i
c
c remove leading spaces
c
	ilen=len(str)
	do 10 i=1,ilen
	    if(str(i:i).ne.' ')goto 11
10	continue
	value=' '
	return
11	tmp=str(i:)
	str=tmp
c
	icma=index(str,',')
	if(icma.eq.0)then
	    value=str
	    str=' '
	else
	    value=str(:icma-1)
	    tmp=str(icma+1:)
	    str=tmp
	endif
	return
	end
