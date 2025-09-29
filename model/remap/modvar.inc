c
c   $Id: modvar.inc,v 1.2 1998/07/21 23:58:27 eds Exp $
c
      integer mdvcde(vardim),mdvlen
      character*15 mdvnme(vardim)
      character*50 mdvstr(vardim)
      character*20 mdvunt(vardim)
      character*10 mdvfmt(vardim)
	  character*80 mdvppl(vardim)
c
      common/varnum/mdvcde,mdvlen
      common/varchr/mdvnme,mdvstr,mdvunt,mdvfmt,mdvppl

