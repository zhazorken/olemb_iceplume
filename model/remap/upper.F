      subroutine upper(string,nc)
      integer nc
      character*(*) string
c
c $Id: upper.F,v 1.1.1.1 1998/07/20 22:03:50 eds Exp $
c
      integer i,cdiff,la,lz,chr
c
      la=ichar('a')
      lz=ichar('z')
      cdiff=ichar('A')-la
c
      do i=1,nc
         chr=ichar(string(i:i))
         if(chr.ge.la .and. chr.le.lz) then
            string(i:i)=char(chr+cdiff)
         endif
      enddo
c
      return
      end
