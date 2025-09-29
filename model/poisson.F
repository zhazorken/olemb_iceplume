#if !defined(NOPOISSON)
      subroutine pois3d (lperod,l,c1,mperod,m,c2,nperod,n,a,b,c,ldimf,
     1                   mdimf,f,ierror,w)
c
c   $Id: poisson.F,v 1.2 1998/07/21 21:34:19 eds Exp $
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                        f i s h p a k                          *
c     *                                                               *
c     *                                                               *
c     *     a package of fortran subprograms for the solution of      *
c     *                                                               *
c     *      separable elliptic partial differential equations        *
c     *                                                               *
c     *                  (version 3.1 , october 1980)                 *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *        john adams, paul swarztrauber and roland sweet         *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the national center for atmospheric research          *
c     *                                                               *
c     *                boulder, colorado  (80307)  u.s.a.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the national science foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * *
c
c     subroutine pois3d solves the linear system of equations
c
c       c1*(x(i-1,j,k)-2.*x(i,j,k)+x(i+1,j,k))
c     + c2*(x(i,j-1,k)-2.*x(i,j,k)+x(i,j+1,k))
c     + a(k)*x(i,j,k-1)+b(k)*x(i,j,k)+c(k)*x(i,j,k+1) = f(i,j,k)
c
c     for  i=1,2,...,l , j=1,2,...,m , and k=1,2,...,n .
c
c     the indices k-1 and k+1 are evaluated modulo n, i.e.
c     x(i,j,0) = x(i,j,n) and x(i,j,n+1) = x(i,j,1). the unknowns
c     x(0,j,k), x(l+1,j,k), x(i,0,k), and x(i,m+1,k) are assumed to take
c     on certain prescribed values described below.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c     lperod   indicates the values that x(0,j,k) and x(l+1,j,k) are
c              assumed to have.
c
c              = 0  if x(0,j,k) = x(l,j,k) and x(l+1,j,k) = x(1,j,k).
c              = 1  if x(0,j,k) = x(l+1,j,k) = 0.
c              = 2  if x(0,j,k) = 0  and x(l+1,j,k) = x(l-1,j,k).
c              = 3  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = x(l-1,j,k).
c              = 4  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = 0.
c
c     l        the number of unknowns in the i-direction. l must be at
c              least 3.
c
c     c1       the real constant that appears in the above equation.
c
c     mperod   indicates the values that x(i,0,k) and x(i,m+1,k) are
c              assumed to have.
c
c              = 0  if x(i,0,k) = x(i,m,k) and x(i,m+1,k) = x(i,1,k).
c              = 1  if x(i,0,k) = x(i,m+1,k) = 0.
c              = 2  if x(i,0,k) = 0 and x(i,m+1,k) = x(i,m-1,k).
c              = 3  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = x(i,m-1,k).
c              = 4  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = 0.
c
c     m        the number of unknowns in the j-direction. m must be at
c              least 3.
c
c     c2       the real constant which appears in the above equation.
c
c     nperod   = 0  if a(1) and c(n) are not zero.
c              = 1  if a(1) = c(n) = 0.
c
c     n        the number of unknowns in the k-direction. n must be at
c              least 3.
c
c
c     a,b,c    one-dimensional arrays of length n that specify the
c              coefficients in the linear equations given above.
c
c              if nperod = 0 the array elements must not depend upon the
c              index k, but must be constant.  specifically,the
c              subroutine checks the following condition
c
c                          a(k) = c(1)
c                          c(k) = c(1)
c                          b(k) = b(1)
c
c                  for k=1,2,...,n.
c
c     ldimf    the row (or first) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  ldimf must be at least l.
c
c     mdimf    the column (or second) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  mdimf must be at least m.
c
c     f        a three-dimensional array that specifies the values of
c              the right side of the linear system of equations given
c              above.  f must be dimensioned at least l x m x n.
c
c     w        a one-dimensional array that must be provided by the
c              user for work space.  the length of w must be at least
c              30 + l + m + 2*n + max(l,m,n) +
c              7*(int((l+1)/2) + int((m+1)/2)).
c
c
c            * * * * * *   on output   * * * * * *
c
c     f        contains the solution x.
c
c     ierror   an error flag that indicates invalid input parameters.
c              except for number zero, a solution is not attempted.
c              = 0  no error
c              = 1  if lperod .lt. 0 or .gt. 4
c              = 2  if l .lt. 3
c              = 3  if mperod .lt. 0 or .gt. 4
c              = 4  if m .lt. 3
c              = 5  if nperod .lt. 0 or .gt. 1
c              = 6  if n .lt. 3
c              = 7  if ldimf .lt. l
c              = 8  if mdimf .lt. m
c              = 9  if a(k) .ne. c(1) or c(k) .ne. c(1) or b(i) .ne.b(1)
c                      for some k=1,2,...,n.
c              = 10 if nperod = 1 and a(1) .ne. 0 or c(n) .ne. 0
c
c              since this is the only means of indicating a possibly
c              incorrect call to pois3d, the user should test ierror
c              after the call.
c
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(n),b(n),c(n),f(ldimf,mdimf,n),
c     arguments      w(see argument list)
c
c     latest         december 1, 1978
c     revision
c
c     subprograms    pois3d,pos3d1,trid,rffti,rfftf,rfftf1,rfftb,
c     required       rfftb1,costi,cost,sinti,sint,cosqi,cosqf,cosqf1
c                    cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,cffti1,
c                    cfftb,cfftb1,passb2,passb3,passb4,passb,cfftf,
c                    cfftf1,passf1,passf2,passf3,passf4,passf,pimach,
c
c     special        none
c     conditions
c
c     common         value
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in july,1977
c
c     algorithm      this subroutine solves three-dimensional block
c                    tridiagonal linear systems arising from finite
c                    difference approximations to three-dimensional
c                    poisson equations using the fourier transform
c                    package sclrfftpak written by paul swarztrauber.
c
c     space          6561(decimal) = 14641(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine pois3d is roughly proportional
c                    to l*m*n*(log2(l)+log2(m)+5), but also depends on
c                    input parameters lperod and mperod.  some typical
c                    values are listed in the table below when nperod=0.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(k) = c(k) = -0.5*b(k) = 1,       k=1,2,...,n
c
c                    and, when nperod = 1
c
c                       a(1) = c(n) = 0
c                       a(n) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine pois was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                    e = max(abs(z(i,j,k)-x(i,j,k)))/max(abs(x(i,j,k)))
c
c                    where the two maxima are taken over i=1,2,...,l,
c                    j=1,2,...,m and k=1,2,...,n, was computed.  the
c                    value of e is given in the table below for some
c                    typical values of l,m and n.
c
c
c                       l(=m=n)   lperod    mperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         16        0         0         272     1.e-13
c                         15        1         1         287     4.e-13
c                         17        3         3         338     2.e-13
c                         32        0         0        1755     2.e-13
c                         31        1         1        1894     2.e-12
c                         33        3         3        2042     7.e-13
c
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos,sin,atan
c     resident
c     routines
c
c     reference      none
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      integer lperod,l,mperod,m,nperod,ldimf,mdimf,ierror,lp,mp,np,k
      integer iwd,iwt,nhm1,nh,nodd,i,j,iwyrt,iwx,iwy,iwbb,n,nhpk,nhmk
      real c1,c2,a,b,c,f,w,save
c
      dimension       a(*)       ,b(*)       ,c(*)       ,
     1                f(ldimf,mdimf,*)       ,w(*)       ,save(6)
      lp = lperod+1
      mp = mperod+1
      np = nperod+1
c
c     check for invalid input.
c
      ierror = 0
      if (lp.lt.1 .or. lp.gt.5) ierror = 1
      if (l .lt. 3) ierror = 2
      if (mp.lt.1 .or. mp.gt.5) ierror = 3
      if (m .lt. 3) ierror = 4
      if (np.lt.1 .or. np.gt.2) ierror = 5
      if (n .lt. 3) ierror = 6
      if (ldimf .lt. l) ierror = 7
      if (mdimf .lt. m) ierror = 8
      if (np .ne. 1) go to 103
      do 101 k=1,n
         if (a(k) .ne. c(1)) go to 102
         if (c(k) .ne. c(1)) go to 102
         if (b(k) .ne. b(1)) go to 102
  101 continue
      go to 104
  102 ierror = 9
  103 if (nperod.eq.1 .and. (a(1).ne.0. .or. c(n).ne.0.)) ierror = 10
  104 if (ierror .ne. 0) go to 122
      iwyrt = l+1
      iwt = iwyrt+m
      iwd = iwt+max0(l,m,n)+1
      iwbb = iwd+n
      iwx = iwbb+n
      iwy = iwx+7*((l+1)/2)+15
      go to (105,114),np
c
c     reorder unknowns when nperod = 0.
c
  105 nh = (n+1)/2
      nhm1 = nh-1
      nodd = 1
      if (2*nh .eq. n) nodd = 2
      do 111 j=1,m
         do 110 i=1,l
            do 106 k=1,nhm1
               nhpk = nh+k
               nhmk = nh-k
               w(k) = f(i,j,nhmk)-f(i,j,nhpk)
               w(nhpk) = f(i,j,nhmk)+f(i,j,nhpk)
  106       continue
            w(nh) = 2.*f(i,j,nh)
            go to (108,107),nodd
  107       w(n) = 2.*f(i,j,n)
  108       do 109 k=1,n
               f(i,j,k) = w(k)
  109       continue
  110    continue
  111 continue
      save(1) = c(nhm1)
      save(2) = a(nh)
      save(3) = c(nh)
      save(4) = b(nhm1)
      save(5) = b(n)
      save(6) = a(n)
      c(nhm1) = 0.
      a(nh) = 0.
      c(nh) = 2.*c(nh)
      go to (112,113),nodd
  112 b(nhm1) = b(nhm1)-a(nh-1)
      b(n) = b(n)+a(n)
      go to 114
  113 a(n) = c(nh)
  114 continue
      call pos3d1 (lp,l,mp,m,n,a,b,c,ldimf,mdimf,f,w,w(iwyrt),w(iwt),
     1             w(iwd),w(iwx),w(iwy),c1,c2,w(iwbb))
      go to (115,122),np
  115 do 121 j=1,m
         do 120 i=1,l
            do 116 k=1,nhm1
               nhmk = nh-k
               nhpk = nh+k
               w(nhmk) = .5*(f(i,j,nhpk)+f(i,j,k))
               w(nhpk) = .5*(f(i,j,nhpk)-f(i,j,k))
  116       continue
            w(nh) = .5*f(i,j,nh)
            go to (118,117),nodd
  117       w(n) = .5*f(i,j,n)
  118       do 119 k=1,n
               f(i,j,k) = w(k)
  119       continue
  120    continue
  121 continue
      c(nhm1) = save(1)
      a(nh) = save(2)
      c(nh) = save(3)
      b(nhm1) = save(4)
      b(n) = save(5)
      a(n) = save(6)
  122 continue
      return
      end
      subroutine pos3d1 (lp,l,mp,m,n,a,b,c,ldimf,mdimf,f,xrt,yrt,t,d,
     1                   wx,wy,c1,c2,bb)
c
      integer lp,l,mp,m,n,ldimf,mdimf,lr,mr,nr,i,j,ifwrd,is,k
      integer lrdel,mrdel
      real a,b,c,f,xrt,yrt,t,d,wx,wy,c1,c2,bb,pi,dum,di,scalx,dj,scaly
      real pimach,dx,dy
c
      dimension       a(*)       ,b(*)       ,c(*)       ,
     1                f(ldimf,mdimf,*)       ,xrt(*)     ,yrt(*)     ,
     2                t(*)       ,d(*)       ,wx(*)      ,wy(*)      ,
     3                bb(*)
      pi = pimach(dum)
      lr = l
      mr = m
      nr = n
c
c     generate transform roots
c
      lrdel = ((lp-1)*(lp-3)*(lp-5))/3
      scalx = lr+lrdel
      dx = pi/(2.*scalx)
      go to (108,103,101,102,101),lp
  101 di = 0.5
      scalx = 2.*scalx
      go to 104
  102 di = 1.0
      go to 104
  103 di = 0.0
  104 do 105 i=1,lr
         xrt(i) = -4.*c1*(sin((float(i)-di)*dx))**2
  105 continue
      scalx = 2.*scalx
      go to (112,106,110,107,111),lp
  106 call sinti (lr,wx)
      go to 112
  107 call costi (lr,wx)
      go to 112
  108 xrt(1) = 0.
      xrt(lr) = -4.*c1
      do 109 i=3,lr,2
         xrt(i-1) = -4.*c1*(sin(float((i-1))*dx))**2
         xrt(i) = xrt(i-1)
  109 continue
      call rffti (lr,wx)
      go to 112
  110 call sinqi (lr,wx)
      go to 112
  111 call cosqi (lr,wx)
  112 continue
      mrdel = ((mp-1)*(mp-3)*(mp-5))/3
      scaly = mr+mrdel
      dy = pi/(2.*scaly)
      go to (120,115,113,114,113),mp
  113 dj = 0.5
      scaly = 2.*scaly
      go to 116
  114 dj = 1.0
      go to 116
  115 dj = 0.0
  116 do 117 j=1,mr
         yrt(j) = -4.*c2*(sin((float(j)-dj)*dy))**2
  117 continue
      scaly = 2.*scaly
      go to (124,118,122,119,123),mp
  118 call sinti (mr,wy)
      go to 124
  119 call costi (mr,wy)
      go to 124
  120 yrt(1) = 0.
      yrt(mr) = -4.*c2
      do 121 j=3,mr,2
         yrt(j-1) = -4.*c2*(sin(float((j-1))*dy))**2
         yrt(j) = yrt(j-1)
  121 continue
      call rffti (mr,wy)
      go to 124
  122 call sinqi (mr,wy)
      go to 124
  123 call cosqi (mr,wy)
  124 continue
      ifwrd = 1
      is = 1
  125 continue
c
c     transform x
c
      do 141 k=1,nr
         do 140 j=1,mr
            do 126 i=1,lr
               t(i) = f(i,j,k)
  126       continue
            go to (127,130,131,134,135),lp
  127       go to (128,129),ifwrd
  128       call rfftf (lr,t,wx)
            go to 138
  129       call rfftb (lr,t,wx)
            go to 138
  130       call sint (lr,t,wx)
            go to 138
  131       go to (132,133),ifwrd
  132       call sinqf (lr,t,wx)
            go to 138
  133       call sinqb (lr,t,wx)
            go to 138
  134       call cost (lr,t,wx)
            go to 138
  135       go to (136,137),ifwrd
  136       call cosqf (lr,t,wx)
            go to 138
  137       call cosqb (lr,t,wx)
  138       continue
            do 139 i=1,lr
               f(i,j,k) = t(i)
  139       continue
  140    continue
  141 continue
      go to (142,164),ifwrd
c
c     transform y
c
  142 continue
      do 158 k=1,nr
         do 157 i=1,lr
            do 143 j=1,mr
               t(j) = f(i,j,k)
  143       continue
            go to (144,147,148,151,152),mp
  144       go to (145,146),ifwrd
  145       call rfftf (mr,t,wy)
            go to 155
  146       call rfftb (mr,t,wy)
            go to 155
  147       call sint (mr,t,wy)
            go to 155
  148       go to (149,150),ifwrd
  149       call sinqf (mr,t,wy)
            go to 155
  150       call sinqb (mr,t,wy)
            go to 155
  151       call cost (mr,t,wy)
            go to 155
  152       go to (153,154),ifwrd
  153       call cosqf (mr,t,wy)
            go to 155
  154       call cosqb (mr,t,wy)
  155       continue
            do 156 j=1,mr
               f(i,j,k) = t(j)
  156       continue
  157    continue
  158 continue
      go to (159,125),ifwrd
  159 continue
c
c     solve tridiagonal systems in z
c
      do 163 j=1,mr
         do 162 i=1,lr
            do 160 k=1,nr
               bb(k) = b(k)+xrt(i)+yrt(j)
               t(k) = f(i,j,k)
  160       continue
            call trid (nr,a,bb,c,t,d)
            do 161 k=1,nr
               f(i,j,k) = t(k)
  161       continue
  162    continue
  163 continue
      ifwrd = 2
      is = -1
      go to 142
  164 continue
      do 167 k=1,nr
         do 166 j=1,mr
            do 165 i=1,lr
               f(i,j,k) = f(i,j,k)/(scalx*scaly)
  165       continue
  166    continue
  167 continue
      return
      end
      subroutine trid (mr,a,b,c,y,d)
c
      integer mr,m,i,mm1,ip
      real a,b,c,y,d,z
c
      dimension       a(*)       ,b(*)       ,c(*)       ,y(*)       ,
     1                d(*)
      m = mr
      mm1 = m-1
      z = 1./b(1)
      d(1) = c(1)*z
      y(1) = y(1)*z
      do 101 i=2,mm1
         z = 1./(b(i)-a(i)*d(i-1))
         d(i) = c(i)*z
         y(i) = (y(i)-a(i)*y(i-1))*z
  101 continue
      z = b(m)-a(m)*d(mm1)
      if (z .ne. 0.) go to 102
      y(m) = 0.
      go to 103
  102 y(m) = (y(m)-a(m)*y(mm1))/z
  103 continue
      do 104 ip=1,mm1
         i = m-ip
         y(i) = y(i)-d(i)*y(i+1)
  104 continue
      return
      end
      subroutine cosqb (n,x,wsave)
c
      integer n
      real x,wsave,tsqrt2,x1
c
      dimension       x(*)       ,wsave(*)
      data tsqrt2 /2.82842712474619/
      if (n-2) 101,102,103
  101 x(1) = 4.*x(1)
      return
  102 x1 = 4.*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
  103 call cosqb1 (n,x,wsave,wsave(n+1))
      return
      end
      subroutine cosqb1 (n,x,w,xh)
c
      integer n,i,ns2,k,modn,np2,kc
      real x,w,xh,xim1
c
      dimension       x(*)       ,w(*)       ,xh(*)
      ns2 = (n+1)/2
      np2 = n+2
      do 101 i=3,n,2
         xim1 = x(i-1)+x(i)
         x(i) = x(i)-x(i-1)
         x(i-1) = xim1
  101 continue
      x(1) = x(1)+x(1)
      modn = mod(n,2)
      if (modn .eq. 0) x(n) = x(n)+x(n)
      call rfftb (n,x,xh)
      do 102 k=2,ns2
         kc = np2-k
         xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
         xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do 103 k=2,ns2
         kc = np2-k
         x(k) = xh(k)+xh(kc)
         x(kc) = xh(k)-xh(kc)
  103 continue
      x(1) = x(1)+x(1)
      return
      end
      subroutine cosqf (n,x,wsave)
c
      integer n
      real x,wsave,sqrt2,tsqx
c
      dimension       x(*)       ,wsave(*)
      data sqrt2 /1.4142135623731/
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
  103 call cosqf1 (n,x,wsave,wsave(n+1))
      return
      end
      subroutine cosqf1 (n,x,w,xh)
c
      integer n,ns2,k,modn,i,np2,kc
      real x,w,xh,xim1
c
      dimension       x(*)       ,w(*)       ,xh(*)
      ns2 = (n+1)/2
      np2 = n+2
      do 101 k=2,ns2
         kc = np2-k
         xh(k) = x(k)+x(kc)
         xh(kc) = x(k)-x(kc)
  101 continue
      modn = mod(n,2)
      if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
      do 102 k=2,ns2
         kc = np2-k
         x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
         x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
      call rfftf (n,x,xh)
      do 103 i=3,n,2
         xim1 = x(i-1)-x(i)
         x(i) = x(i-1)+x(i)
         x(i-1) = xim1
  103 continue
      return
      end
      subroutine cosqi (n,wsave)
c
      integer n,k
      real wsave,pih,dt,fk
c
      dimension       wsave(*)
      data pih /1.57079632679491/
      dt = pih/float(n)
      fk = 0.
      do 101 k=1,n
         fk = fk+1.
         wsave(k) = cos(fk*dt)
  101 continue
      call rffti (n,wsave(n+1))
      return
      end
      subroutine cost (n,x,wsave)
c
      integer n,nm1,k,ns2,i,modn,np1,kc
      real x,wsave,x1h,x1p3,tx2,t1,t2,xim2,c1,xi
c
      dimension       x(*)       ,wsave(*)
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
      x(2) = x(1)-x(2)
      x(1) = x1h
      return
  102 if (n .gt. 3) go to 103
      x1p3 = x(1)+x(3)
      tx2 = x(2)+x(2)
      x(2) = x(1)-x(3)
      x(1) = x1p3+tx2
      x(3) = x1p3-tx2
      return
  103 c1 = x(1)-x(n)
      x(1) = x(1)+x(n)
      do 104 k=2,ns2
         kc = np1-k
         t1 = x(k)+x(kc)
         t2 = x(k)-x(kc)
         c1 = c1+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(k) = t1-t2
         x(kc) = t1+t2
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
      call rfftf (nm1,x,wsave(n+1))
      xim2 = x(2)
      x(2) = c1
      do 105 i=4,n,2
         xi = x(i)
         x(i) = x(i-2)-x(i-1)
         x(i-1) = xim2
         xim2 = xi
  105 continue
      if (modn .ne. 0) x(n) = xim2
  106 return
      end
      subroutine costi (n,wsave)
c
      integer n,ns2,kc,nm1,np1,k
      real wsave,pi,fk,dt
c
      dimension       wsave(*)
      data pi /3.14159265358979/
      if (n .le. 3) return
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/float(nm1)
      fk = 0.
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.
         wsave(k) = 2.*sin(fk*dt)
         wsave(kc) = 2.*cos(fk*dt)
  101 continue
      call rffti (nm1,wsave(n+1))
      return
      end
      real function pimach (dum)
c
c     this subprogram supplies the value of the constant pi correct to
c     machine precision where
c
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
      real dum
c
      pimach = 3.14159265358979
      return
      end
      subroutine rfftb (n,r,wsave)
c
      integer n
      real r,wsave
c
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftb1 (n,c,ch,wa,ifac)
c
      integer n,ifac,nf,na,l1,iw,k1,ip,ido,l2,ix2,ix3,ix4,idl1,i
      real c,ch,wa
c
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call radb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na .ne. 0) go to 113
         call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
      return
      end
      subroutine rfftf (n,r,wsave)
c
      integer n
      real r,wsave
c
      dimension       r(*)       ,wsave(*)
      if (n .eq. 1) return
      call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rfftf1 (n,c,ch,wa,ifac)
c
      integer n,ifac,nf,na,l2,iw,k1,ip,kh,idl1,ido,l1,ix2,ix3,ix4,i
      real c,ch,wa
c
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
      subroutine rffti (n,wsave)
c
      integer n
      real wsave
c
      dimension       wsave(*)
      if (n .eq. 1) return
      call rffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end
      subroutine rffti1 (n,wa,ifac)
c
      integer n,ifac,ntryh,nl,nf,j,ntry,nr,nq,i,ib,is,l1,nfm1,k1
      integer ip,ld,ido,l2,ipm,ii
      real wa,tpi,argld,argh,fi,arg
c
      dimension       wa(*)      ,ifac(*)    ,ntryh(4)
      data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      tpi = 6.28318530717959
      argh = tpi/float(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = float(ld)*argh
            fi = 0.
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      return
      end
      subroutine sinqb (n,x,wsave)
c
      integer n,ns2,k,kc
      real x,wsave,xhold
c
      dimension       x(*)       ,wsave(*)
      if (n .gt. 1) go to 101
      x(1) = 4.*x(1)
      return
  101 ns2 = n/2
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      call cosqb (n,x,wsave)
      do 103 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  103 continue
      return
      end
      subroutine sinqf (n,x,wsave)
c
      integer n,ns2,k,kc
      real x,wsave,xhold
c
      dimension       x(*)       ,wsave(*)
      if (n .eq. 1) return
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  101 continue
      call cosqf (n,x,wsave)
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      return
      end
      subroutine sinqi (n,wsave)
c
      integer n
      real wsave
c
      dimension       wsave(*)
      call cosqi (n,wsave)
      return
      end
      subroutine sint (n,x,wsave)
c
      integer n,np1,iw1,iw2,iw3
      real x,wsave
c
      dimension       x(*)       ,wsave(*)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      iw3 = iw2+np1
      call sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
      return
      end
      subroutine sint1(n,war,was,xh,x,ifac)
c
      integer n,ifac,i,np1,ns2,k,kc,modn
      real war,was,xh,x,sqrt3,xhold,t1,t2
c
      dimension war(*),was(*),x(*),xh(*),ifac(*)
      data sqrt3 /1.73205080756888/
      do 100 i=1,n
      xh(i) = war(i)
      war(i) = x(i)
  100 continue
      if (n-2) 101,102,103
  101 xh(1) = xh(1)+xh(1)
      go to 106
  102 xhold = sqrt3*(xh(1)+xh(2))
      xh(2) = sqrt3*(xh(1)-xh(2))
      xh(1) = xhold
      go to 106
  103 np1 = n+1
      ns2 = n/2
      x(1) = 0.
      do 104 k=1,ns2
         kc = np1-k
	 t1 = xh(k)-xh(kc)
	 t2 = was(k)*(xh(k)+xh(kc))
	 x(k+1) = t1+t2
	 x(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+2) = 4.*xh(ns2+1)
      call rfftf1 (np1,x,xh,war,ifac)
      xh(1) = .5*x(1)
      do 105 i=3,n,2
	 xh(i-1) = -x(i)
	 xh(i) = xh(i-2)+x(i-1)
  105 continue
      if (modn .ne. 0) go to 106
      xh(n) = -x(n+1)
  106 do 107 i=1,n
      x(i) = war(i)
      war(i) = xh(i)
  107 continue
      return
      end
      subroutine sinti (n,wsave)
c
      integer n,ns2,np1,k
      real wsave,pi,dt
c
      dimension       wsave(*)
      data pi /3.14159265358979/
      if (n .le. 1) return
      ns2 = n/2
      np1 = n+1
      dt = pi/float(np1)
      do 101 k=1,ns2
         wsave(k) = 2.*sin(k*dt)
  101 continue
      call rffti (np1,wsave(ns2+1))
      return
      end
      subroutine radb2 (ido,l1,cc,ch,wa1)
c
      integer ido,l1,k,i,ic,idp2
      real cc,ch,wa1,tr2,ti2
c
      dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
      end
      subroutine radb3 (ido,l1,cc,ch,wa1,wa2)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,wa2,taur,taui
      real dr3,di3,tr2,ti2
      real cr2,ci2,cr3,ci3,dr2,di2
c
      dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
      return
      end
      subroutine radb4 (ido,l1,cc,ch,wa1,wa2,wa3)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,wa2,wa3,ti4,tr2,tr3,cr2,sqrt2,tr1,ti1
      real cr3,ci4,ti2,ti3,tr4,ci3,cr4,ci2
c
      dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data sqrt2 /1.414213562373095/
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
      end
      subroutine radb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,wa2,wa3,wa4,tr11,ti11,tr12,ti12,ti5,ti4,tr2,tr3,cr2
      real ci5,cr3,ci4,ti2,ti3,tr5,tr4,di3,ci3,cr4,di5,ci2,cr5,dr2
      real di2,dr4,di4,dr3,dr5
c
      dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
      return
      end
      subroutine radbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
      integer ido,ip,l1,idl1,idp2,nbd,ipp2,ipph,k,i,j,j2,jc,ic,l
      integer ik,is,idij,lc
      real cc,c1,c2,ch,ch2,wa,tpi,arg,dcp,dsp,ar1,ai1,ar1h,dc2,ds2
      real ar2,ai2,ar2h
c
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
  103 do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 ik=1,idl1
         do 121 j=2,ipph
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 k=1,l1
            do 129 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
      do 135 k=1,l1
         do 134 j=2,ip
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 return
      end
      subroutine radf2 (ido,l1,cc,ch,wa1)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,tr2,ti2
c
      dimension       ch(ido,2,l1)           ,cc(ido,l1,2)           ,
     1                wa1(*)
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
      end
      subroutine radf3 (ido,l1,cc,ch,wa1,wa2)
c
      integer ido,l1,k,i,idp2,ic
      real cc,ch,wa1,wa2,taur,taui,cr2,dr2,di2,dr3,di3,ci2,tr2
      real tr3,ti2,ti3
c
      dimension       ch(ido,3,l1)           ,cc(ido,l1,3)           ,
     1                wa1(*)     ,wa2(*)
      data taur,taui /-.5,.866025403784439/
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
      return
      end
      subroutine radf4 (ido,l1,cc,ch,wa1,wa2,wa3)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,wa2,wa3,hsqt2,tr1,tr2,cr2,ci2,cr3,ci3,cr4,ci4
      real ti1,ti2,ti4,tr3,tr4,ti3
c
      dimension       cc(ido,l1,4)           ,ch(ido,4,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)
      data hsqt2 /.7071067811865475/
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
      end
      subroutine radf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
c
      integer ido,l1,k,idp2,i,ic
      real cc,ch,wa1,wa2,wa3,wa4,tr11,ti11,tr12,ti12,cr2,ci5,cr3,ci4
      real dr2,di2,dr3,di3,dr4,di4,dr5,di5,ci2,ci3,tr2,tr5,tr3,tr4
      real ti3,ti4,cr5,cr4,ti2,ti5
c
      dimension       cc(ido,l1,5)           ,ch(ido,5,l1)           ,
     1                wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
      data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
      return
      end
      subroutine radfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
      integer ido,ip,l1,idl1,ipph,ipp2,idp2,nbd,ik,j,k,is,idij,i,jc
      integer l,j2,ic,lc
      real cc,c1,c2,ch,ch2,wa,tpi,arg,dcp,dsp,ar1,ai1,ar1h,dc2,ds2
      real ar2,ai2,ar2h
c
      dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          ,
     1                c1(ido,l1,ip)          ,c2(idl1,ip),
     2                ch2(idl1,ip)           ,wa(*)
      data tpi/6.28318530717959/
      arg = tpi/float(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 k=1,l1
            do 116 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.
      ai1 = 0.
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 k=1,l1
         do 133 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
      return
      end
#endif
