c
c    $Id: dentable.inc,v 1.1.1.1 1998/07/20 22:03:48 eds Exp $
c
	real dt,dss
	integer tsize,ssize,tssize
c
        parameter (dt=0.02,dss=0.05,tsize=300,ssize=100,
     *             tssize=tsize*ssize)
c       parameter (dt=0.2,dss=0.5,tsize=200,ssize=80,tssize=tsize*ssize)
c
	real dtemp,dsal
	integer tempsize,salsize
c
	common/iessize/dtemp,dsal,tempsize,salsize
c
	real tmin,tmax,smin,smax
	real rhotbl(tsize*ssize*2*kdim1),zlev(2*kdim1)
	integer nt,ns,nz
c
	common/iestable/tmin,tmax,nt,smin,smax,ns,zlev,nz,
     *	rhotbl
