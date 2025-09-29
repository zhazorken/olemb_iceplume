c
c   $Id: refdate.inc,v 1.1.1.1 1998/07/20 22:03:48 eds Exp $
c
c  refdate - reference date for model run
c  cardref - reference date read from input
c    bcref - reference date from bc file
c    icref - reference date from ic or restart file
c
c datemode - if true, use date system
c carddate - if true, date read from input
c   icdate - if true, date in ic or restart file
c   bcdate - if true, date in bc file
c
c   usedate - if true, use date system (if possible)
c useicdate - use date from ic file (overriding card input)
c
      integer refdate(2),cardref(2),bcref(2),icref(2)
      logical datemode,usedate,useicdate
      logical carddate,icdate,bcdate
c
      common/refdat/refdate,cardref,bcref,icref,datemode,
     *     usedate,useicdate,carddate,icdate,bcdate
