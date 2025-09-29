c
c   $Id: ies80.F,v 1.2 1998/07/24 16:59:40 dwd Exp $
c
c International Equation of State of Sea Water, 1980
c ref: Millero and Poisson (1981) and Unesco Technical Paper
c in Marine Science Number 36 (Unesco, 1981).
c
c IES80 is computed in the given form
c
c rho(s,t,p)=rho(s,t,0)/[1-p/K(s,t,p)]
c
c where, K(s,t,p) is the secant bulk modulus,
c        s -- salinity in PSS-78
c        t -- temperature in degrees C
c        p -- pressure (bars)
c
c        rho -- density (kg m-3)
c
c values for checking:
c
c   S   T(C)  p(bars)    rho(s,t,p)      K(s,t,p)
c
c   0     5        0      999.96675    20337.80375
c               1000     1044.12802    23643.52599
c        25        0      997.04796    22100.72106
c               1000     1037.90204    25405.09717
c  35     5        0     1027.67547    22185.93358
c               1000     1069.48914    25577.49819
c        25        0     1023.34306    23726.34949
c               1000     1062.53817    27108.94504
c
c
      real function sig0(s,t)
c
c     rho0(s,t)=1000+sig0(s,t)
c
      real s,t
c
      real rt,rts,rts32
c
      rt=((((6.536332e-9*t-1.120083e-6)*t+1.001685e-4)*t-9.09529e-3)*t
     *     +6.793952e-2)*t
      rts=(((5.3875e-9*t-8.2467e-7)*t+7.6438e-5)*t-4.0899e-3)*t
     *     +8.24493e-1
      rts32=(-1.6546e-6*t+1.0227e-4)*t-5.72466e-3
c
      sig0=-1.57406e-1+rt+(4.8314e-4*s+rts)*s+rts32*sqrt(s*s*s)
      return
      end
c
      real function sbm(s,t,p)
c
c     K(s,t,p)=20000+sbm(s,t,p)
c
      real s,t,p
c
      real rt,rtp,rtp2,rts,rts32,rtps,rtp2s
c
      rt=(((-5.155288e-5*t+1.360477e-2)*t-2.327105)*t+148.4206)*t
      rtp=((-5.77905e-7*t+1.16092e-4)*t+1.43713e-3)*t+3.239908
      rtp2=(5.2787e-8*t-6.12293e-6)*t+8.50935e-5
      rts=((-6.167e-5*t+1.09987e-2)*t-0.603459)*t+54.6746
      rts32=(-5.3009e-4*t+1.6483e-2)*t+7.944e-2
      rtps=(-1.6078e-6*t-1.0981e-5)*t+2.2838e-3
      rtp2s=(9.1697e-10*t+2.0816e-8)*t-9.9348e-7
c
      sbm=-347.79+rt
     *     +(rtp2*p+rtp)*p
     *     +(rts+(rtp2s*p+rtps)*p)*s
     *     +(1.91075e-4*p+rts32)*sqrt(s*s*s)
      return
      end
c
      real function sig(s,t,p)
c
c     rho(s,t,p)=1000+sig(s,t,p)
c
c     s   salinity in ppt
c     t   temperature in C
c     p   pressure in bars  (approx.  depth*0.1)
c
#if defined(cray)
      real k,pd
#else
      real*8 k,pd
#endif
      real s,t,p,sbm,sig0
c
      k=2.0d4+sbm(s,t,p)
      pd=p
      sig=(1.0d3+sig0(s,t))/(1.-pd/k)-1.0d3
      return
      end
c
      real function sigthta(s,thta,p)
c
c     sig for potential temperature
c
c     salinity in PSS-78
c     potential temperature (C)
c     pressure bars  (approx. depth*0.1)
c
      real s,thta,p,ipts,t,sig
c
      t=ipts(s,thta,p*10.)
c
      sigthta=sig(s,t,p)
      return
      end
c
      subroutine mkies(rtmin,rtmax,rsmin,rsmax,zl,nrz)
c
c     potential temperature
c     salinity in ppt
c     depth in m
c
#include "param.inc"
#include "olemp.inc"
#include "dentable.inc"
c
      integer nrz
      real rtmin,rtmax,rsmin,rsmax,zl(nrz)
c
      integer k,j,i,index
      real p,s,t,thta
c
      real rhobase,s0,t0,drdt,d2rdt2,drds,d2rds2,d2rdtds
c
      common/rholin/rhobase,s0,t0,drdt,d2rdt2,drds,d2rds2,d2rdtds
c
      real rh,ipts,sbm
      real sigthta
c
      data dtemp,dsal,tempsize,salsize/dt,dss,tsize,ssize/
c
      tmin=rtmin
      smin=rsmin
      nt=min(tsize,int((rtmax-tmin)/dtemp)+1)
      ns=min(ssize,int((rsmax-smin)/dsal)+1)
      nz=min(2*kdim1,nrz)
      tmax=tmin+dtemp*(nt-1)
      smax=smin+dsal*(ns-1)
c
      if(appden) call rholset(censal, centhta)
c
      do 10 k=1,nz
         zlev(k)=zl(k)
c
c very crude conversion of depth in meters to bars
c
         p=zlev(k)*0.1
         do 20 j=1,ns
            s=smin+(j-1)*dsal
            do 30 i=1,nt
               thta=tmin+(i-1)*dtemp
               index = i+(j-1)*tsize+(k-1)*tssize
#if defined(MARSH_DEN)
               rhotbl(index)=0.2*(centhta-thta)
#else
               index = i + (j-1)*tsize + (k-1)*tssize
               if(appden)then
                  if(potden)then
                     t=thta
                  else
                     t=ipts(s,thta,p*10.0)
                  endif
                  rh=rhobase+drdt*(t-t0)+drds*(s-s0)
                  if(.not.linden)then
                     rh=rh + 0.5*(d2rdt2*(t-t0)**2
     *                     + d2rdtds*(t-t0)*(s-s0)+d2rds2*(s-s0)**2)
                  endif
                  if(potden)then
                     rhotbl(index) = rh
                  else
                     rhotbl(index) = (1.0d3+rh)/(1.0d0-p/(2.0d4+
     *                    sbm(s,t,p)))-1.0d3
                  endif
               else
                  rhotbl(index) = sigthta(s,thta,p)
               endif
#endif
 30         continue
 20      continue
 10   continue
      return
      end
c
      real function iesfast(s,t,k)
c
c does not do the range check!!!!
c
c     s   salinity in ppt
c     t   potential temp
c     k   level from mkies
c
#include "param.inc"
#include "dentable.inc"
c
      real s,t
      integer k
c
      integer i,j,index
      real y1,y2,y3,y4,ti,si
      real tmpt,tmps
c
      tmpt=(t-tmin)/dtemp
      tmps=(s-smin)/dsal
c
      i=ifix(tmpt)+1
      j=ifix(tmps)+1

      index=i+(j-1)*tsize+(k-1)*tssize
      y1=rhotbl(index)
      y2=rhotbl(index+1)
      y3=rhotbl(index+1+tsize)
      y4=rhotbl(index+tsize)
      ti=tmpt-(i-1)
      si=tmps-(j-1)
      iesfast=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
      return
      end
c
      real function iestbl(s,t,k)
c
c     s   salinity in ppt
c     t   potential temp
c     k   level from mkies
c
#include "param.inc"
#include "dentable.inc"
c
      real s,t
      integer k
c
      integer i,j,index
      real sigthta,y1,y2,y3,y4,ti,si
      real tmpt,tmps
c
      tmpt=(t-tmin)/dtemp
      tmps=(s-smin)/dsal
c
      i=ifix(tmpt)+1
      j=ifix(tmps)+1
c
      if(i.ge.nt .or. i.lt.1 .or. j.ge.ns .or. j.lt.1 .or. k.gt.nz
     *     .or. k.lt.1)then
c
c includes crude conversion from depth in meters to bars
c
         iestbl=sigthta(s,t,zlev(k)*0.1)
      else
         index=i+(j-1)*tsize+(k-1)*tssize
	 y1=rhotbl(index)
	 y2=rhotbl(index+1)
	 y3=rhotbl(index+1+tsize)
	 y4=rhotbl(index+tsize)
         ti=tmpt-(i-1)
         si=tmps-(j-1)
         iestbl=y1+(y4-y1)*si+(y2-y1)*ti+(y1-y2+y3-y4)*si*ti
      endif
      return
      end
c
      real function tf (s,p)
c
c************************************
c function to compute the freezing point of seawater
c
c--reference:  unesco tech. papers in the marine science no. 28
c		1978 eighth report jpots
c		annex 6 freezing point of seawater, f.j.millero,29-35
c
c--units:
c	pressure	p	decibars
c	salinity	s	pss-78
c	temperature	tf	degrees celsius
c	freezing pt.
c
c--checkvalue:
c	tf = -2.588567 deg c, for s=40, p=500 db.
c
c--from fofonoff and millard,1983,draft. may 84, ljm.
c
c***************************************
c
      real s,p
c
      tf = (-0.0575 + 1.710523e-3*sqrt(abs(s))-2.154996e-4*s)
     1     *s-7.53e-4*p
      return
      end
c
      real function ipts(s,thta,p)
      real s,thta,p
c
c Use Newton's method to find ipts temperature from
c potential temperature (referenced to surface)
c
c     pressure     p    decibars
c     salinity     s    pss-78
c     temperutare  thta potential temperature (pr=0)
c
      real eps,delta,pr,ta,tb,t,tnew,dfdt,theta
c
      eps=1.e-5
      delta=1.e-3
      pr=0.0
c
      t=thta
c
 100  ta=theta(s,t,p,pr)
      tb=theta(s,t+delta,p,pr)
      dfdt=(tb-ta)/delta
      tnew=t-(ta-thta)/dfdt
      if(abs(t-tnew) .le. eps)goto 200
      t=tnew
      goto 100
c
 200  ipts=tnew
      return
      end
c
      real function atg(s,t,p0) 
      real s,t,p0
c
c atg fcn *** oct 20 1980 ****
c **************************** 
c adiabatic temperature gradient deg c per decibar
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c units:      
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        psu (ipss-78)
c       adiabatic       atg      deg. c/decibar
c check: atg=3.255976e-4 c/dbar for s=40 psu,t=40 deg c,p0=10000 decibar
c   scale pressure to bars
c
      real p,ds
c
      p=p0/10.
c**************************
      ds = s - 35.0 
      atg = (((-2.1687e-13*t+1.8676e-11)*t-4.6206e-10)*p
     *     +((2.7759e-10*t-1.1351e-8)*ds+((-5.4481e-12*t
     *     +8.733e-10)*t-6.7795e-8)*t+1.8741e-6))*p 
     *     +(-4.2393e-7*t+1.8932e-5)*ds
     *     +((6.6228e-9*t-6.836e-7)*t+8.5258e-5)*t+3.5803e-4
c     scale atg to per decibar
      atg=0.1*atg
c****************
      return  
      end     
c
      real function theta(s,t0,p0,pr)
      real s,t0,p0,pr
c
c theta fcn ****** oct 20 1980 ******
c ***********************************
c to compute local potential temperature at pr
c using bryden 1973 polynomial for adiabatic lapse rate
c and runge-kutta 4-th order integration algorithm.
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c fofonoff,n.,1977,deep-sea res.,24,489-491
c units:      
c       pressure        p0       decibars
c       temperature     t0       deg celsius (ipts-68)
c       salinity        s        psu (ipss-78)
c       reference prs   pr       decibars
c       potential tmp.  theta    deg celsius 
c checkvalue: theta= 36.89072 c,s=40 psu,t0=40 deg c,
c p0=10000 decibars,pr=0 decibars
c             
c      set-up intermediate temperature and pressure variables
c
      real p,t,xk,h,atg,q
c
      p=p0
      t=t0
c**************
      h = pr - p
      xk = h*atg(s,t,p) 
      t = t + 0.5*xk
      q = xk  
      p = p + 0.5*h 
      xk = h*atg(s,t,p) 
      t = t + 0.29289322*(xk-q) 
      q = 0.58578644*xk + 0.121320344*q 
      xk = h*atg(s,t,p) 
      t = t + 1.707106781*(xk-q)
      q = 3.414213562*xk - 4.121320344*q
      p = p + 0.5*h 
      xk = h*atg(s,t,p) 
      theta = t + (xk-2.0*q)/6.0
      return  
      end     
c
      subroutine rholset(scen,tcen)
      real scen,tcen
c
      real rhobase,s0,t0,drdt,d2rdt2,drds,d2rds2,d2rdtds
c
      common/rholin/rhobase,s0,t0,drdt,d2rdt2,drds,d2rds2,d2rdtds
c
      real a0,a1,a2,a3,a4,a5
      real b0,b1,b2,b3,b4
      real c0,c1,c2
      real d2
      real bt,ct,bdt,cdt,s32
c
      a0=-1.57406e-1
      a1=6.793952e-2
      a2=-9.09529e-3
      a3=1.001685e-4
      a4=-1.120083e-6
      a5=6.536332e-9
c
      b0=8.24493e-1
      b1=-4.0899e-3
      b2=7.6438e-5
      b3=-8.2467e-7
      b4=5.3875e-9
c
      c0=-5.72466e-3
      c1=1.0227e-4
      c2=-1.6546e-6
c
      d2=4.8314e-4
c
      t0=tcen
      s0=scen
c
      bt=(((b4*t0+b3)*t0+b2)*t0+b1)*t0+b0
      ct=(c2*t0+c1)*t0+c0
c
      bdt=((4.*b4*t0+3.*b3)*t0+2.*b2)*t0+b1
      cdt=c1+2.*c2*t0
c
      s32=sqrt(s0*s0*s0)
c
      rhobase=a0+((((a5*t0+a4)*t0+a3)*t0+a2)*t0+a1)*t0 +
     *     bt*s0 + ct*s32 + d2*s0*s0
c
      drdt=a1+(((5.*a5*t0+4.*a4)*t0+3.*a3)*t0+2.*a2)*t0 +
     *     bdt*s0 + cdt*s32
c
      d2rdt2=2.*a2+((20.*a5*t0+12.*a4)*t0+6.*a3)*t0 + 
     *     ((12.*b4*t0+6.*b3)*t0+2.*b2)*s0 + 2.*c2*s32
c
      drds=bt + 1.5*ct*sqrt(s0) + 2.*d2*s0
c
      d2rds2= 0.75*ct/sqrt(s0) + 2.*d2
c
      d2rdtds = bdt + 1.5*cdt*sqrt(s0)
c
      return
      end

