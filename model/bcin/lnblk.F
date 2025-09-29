      integer function lnblk (str,ilen)
c
c   $Id: lnblk.F,v 4.1 1995/04/12 21:33:46 eds Exp $
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
