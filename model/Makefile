#
#  $Id: Makefile.in,v 1.10 1998/10/23 16:08:54 dwd Exp $
#
# Makefile for olem models
#
# Model conditional compiliation flags
#
# NEWSOLVE   use mpi solver (mpi version only)
# RADIATE    use radiation conditions
# SNGLSALT   use single precsion salt (cray version)
# NOPOISSON  use other code for pois3d and poistg (cray version)
# ALLOC_DISK allocate all the disk space at the start of the run
# MARSH_DEN  use rho=rho0(1+alpha(t0-t)) equation of state
# ANALY      create analysis fields
# DEBUG      include code for writing maximum pressure
# TKESHEAR   remove domain average shear from tke generation calculation
# CONSGRAD   constant gradient bottom B.C. for u and v
# SUNOS      define for SUNOS 4.1.x compilation for proper signal definition
# NDEBUG     debug turned off (Win32)
# PRINTDEBUG turns on lots of print statements
#
# Choose at ONE of the following
# 
# LEVEQ      use leveq77.F for advection
# LEER       use leer.F for advection
#
# These are options for Surface Wave Paramterization
#
# STOKES     use stokes drift vortex force
# STOKESF    use stokes drift coriolis force (STOKES must be defined) if using uvwflux
# STOKES1    define stokes drift using random perturbation
#
# cray       Cray f77 environment
# sun        Sun compiler f77
# Win32      Microsoft Fortran Powerstation 1.0 (f77)
# 
# CPP=/usr/mpi/intel/openmpi-1.2.5/bin/mpicc -E
CPP=cc
FC=mpif90
# FC=ftn
#
# set model compilation options (lines should be <= 60 chars)
#
# MODFLAGS1= -DUVWFLUX -DLEER
MODFLAGS1= -DUVWFLUX -DLEER -DOPENBC -DCONRES -DHTOP -DMELT
MODFLAGS2= -DSNGLSALT -DANALY -DUSEMPI -DVARDX -DAVEFLUX
# MODFLAGS2= -DSNGLSALT -DANALY -DUSEMPI -DVARDX
# MODFLAGS3= -DLEAPF -DTHREEEQ -DDISCHARGE
MODFLAGS3= -DLEAPF -DSCHULZ -DDISCHARGE -DCHANNEL
# MODFLAGS3= -DLEAPF -DTHREEEQ -DCHANNEL
# MODFLAGS3 = -DHBOT
# MODFLAGS3= -DHBOT -DOPENBC -DOPENR -DRECIRC
# MODFLAGS3= -DHBOT -DOPENBC -DOPENR
#
# set fortran compilation options
#
#
# INC_MPI    := /usr/mpi/intel/openmpi-1.2.5/include
# LIB_MPI    := /usr/mpi/intel/openmpi-1.2.5/lib64
# INC_NETCDF := /home/skylling/netcdf/include
# LIB_NETCDF := /home/skylling/netcdf/lib
#
#
FFLAGS=  -O3 -traceback -g  -init=arrays -init=zero
# FFLAGS=  -O3 -check bounds -traceback -g  -init=arrays -init=zero

# FFLAGS= -O3 -traceback -fpe0 -g
# FFLAGS= -O3 -traceback -g
# FFLAGS= -O3 -march=corei7 -axAVX -traceback
#
# MODFLAGS4=$(FFLAGS)
#
ROOT_DIR=/usr/local
#
# LIBS= -L$(ROOT_DIR)/lib -lnetcdf
LIBS= -lnetcdff
#
# Common sun f77 switches:
# -u         Undefined variables cause fatal error.
# -g         Debug.
# -fast      Fast execution. Optimize for speed.
# -parallel  Automatic and explicit loop parallelizations.
# -stackvar  Local variables are stored on the stack
# -O4        Auto inline of functions.  Highest optimization.
# -sb        Produce table information for SourceBrowser
# -C         includes code for checking array bounds
# -Usym      Undefine preprocessor symbol
# -Dsym      Define preprocessor symbol
#
LNFLAGS=
CPPFLAGS= $(MODFLAGS1) $(MODFLAGS2) $(MODFLAGS3)

SOURCES = bc_read.F bc_update.F bcinit.F dag_create.F \
          dag_write.F date.F deardorff.F deltime.F density.F extrap.F\
          denupdte.F filter.F fld_create.F  conresdx8.F\
          fmetais.F LegendreUAVA.F\
          hostnme.F ic_read.F icinit.F ies80.F imag.F imag8.F\
          image.F in_read.F indexx.F initrad.F intg.F \
          leernew.F leerdan.F lnblk.F \
          loop3dpr.F mks.F myhandler.F olem_init.F obcbox.F\
          parse.F ran2.F restart.F readUV.F\
          sig_handler.F slb_create.F write3dflx.F\
          slb_write.F smag.F topbotb.F topbotdp.F unix.F vol2slp.F wrflux.F\
          upper.F uvwenst.F varcde.F wallice.F wrparm.F zmixed.F


MPI_SOURCES = global.F iompi.F setupmpi.F shifts.F


ALL_SOURCES = $(SOURCES) $(MPI_SOURCES)


OBJECTS= $(ALL_SOURCES:.F=.o)

olem: info $(OBJECTS) olem.o
	$(FC) $(CPPFLAGS) $(FFLAGS) $(LIBS) -o olem230_v_d75w105x2_s0 *.o
#$(OBJECTS):
#	rm -f $*.f
#	$(CPP) $(CPPFLAGS) $*.F | grep -v '^#' | grep -v '^$$' > $*.f
#	$(FC) $(FFLAGS) -c $*.f

olem.o: olem.F
	@sed -e "s;MODEL_FLAGS1;$(MODFLAGS1);" olem.F | \
	sed -e "s;MODEL_FLAGS2;$(MODFLAGS2);" | \
	sed -e "s;MODEL_FLAGS3;$(MODFLAGS3);" | \
	sed -e "s;MODEL_FLAGS4;$(MODFLAGS4);" > foo.F
#	$(CPP) $(CPPFLAGS) foo.F | grep -v '^#' | grep -v '^$$' > foo.f
	$(FC) $(CPPFLAGS) $(FFLAGS)  -c foo.F
	@mv foo.o olem.o



PHAKE:

clean:
	-@rm -f .sb/Locked/\*
	-@rm -f .sb/NewRoot/\*
	-@rm -f .sb/OldRoot/\*
	-@rm -f .sb/Refd/\*
	-@rm -f .sb/Index
	rm -f foo.F *~ *.o *.for *.f *.s *.lst .make.state .nse_depinfo

check:
	cvscheck

info:
	@head -3 param.inc
	@echo " "
	@echo "  model flags: $(MODFLAGS1)"
	@echo "  model flags: $(MODFLAGS2)"
	@echo "  model flags: $(MODFLAGS3)"
	@echo "fortran flags: $(MODFLAGS4)"
	@echo " "
