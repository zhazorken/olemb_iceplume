      integer function lnblk (str,ilen)
c
c   $Id: lnblk.F,v 1.1.1.1 1998/07/20 22:03:48 eds Exp $
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
