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
# STOKESF    use stokes drift coriolis force (STOKES must be defined)
# STOKES1    define stokes drift using random perturbation
#
# cray       Cray f77 environment
# sun        Sun compiler f77
# Win32      Microsoft Fortran Powerstation 1.0 (f77)
# 


CPP=/usr/bin/cpp
FC=mpif77
#
# set model compilation options (lines should be <= 60 chars)
#
MODFLAGS1= -DUVWFLUX -DLEER -DCONRES -DPRINTDEBUG
MODFLAGS2= -DSNGLSALT -DANALY -DUSEMPI
MODFLAGS3= -DFRONT -DLEAPF -DSTOKES
# MODFLAGS3 = -DHBOT
# MODFLAGS3= -DHBOT -DOPENBC -DOPENR -DRECIRC
# MODFLAGS3= -DHBOT -DOPENBC -DOPENR
#
# set fortran compilation options
#
# FFLAGS= -u -C -g -qstrict -bmaxdata:2000000000 -bmaxstack:256000000 -q32
# FFLAGS= -u -O2 -qstrict -bmaxdata:2000000000 -bmaxstack:256000000 -q32
FFLAGS = -u -O2 -traceback
# FFLAGS = -u -fpstkchk -g -traceback
MODFLAGS4=$(FFLAGS)
#
ROOT_DIR=/mnt/home/skylling/local
# LIBS= $(ROOT_DIR)/libnetcdf.a
# LIBS= -L$(ROOT_DIR)/lib32/r4i4 -lnetcdf
LIBS= -L$(ROOT_DIR)/lib -lnetcdf
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
CPPFLAGS= $(MODFLAGS1) $(MODFLAGS2) $(MODFLAGS3) -I$(ROOT_DIR)/include -traditional-cpp


SOURCES =  \
          date.F deltime.F density.F \
          denupdte.F \
          dag_write.F filter.F\
          hostnme.F ic_read.F icinit.F ies80.F imag2.F \
          image.F in_read.F indexx.F initrad.F intg.F \
          lnblk.F imag2f.F extrap.F outfile.F\
          mks.F myhandler.F olem_init.F \
          parse.F ran2.F \
          sig_handler.F \
          topbotb.F topbotdp.F unix.F \
          upper.F varcde.F wrparm.F zmixed.F


MPI_SOURCES = global.F iompi.F setupmpi.F shifts.F


ALL_SOURCES = $(SOURCES) $(MPI_SOURCES)


OBJECTS= $(ALL_SOURCES:.F=.o)

remap: info $(OBJECTS) remap.o
	$(FC) $(LNFLAGS) $(FFLAGS)  -o remap1250 *.o libfftw.a $(LIBS)
$(OBJECTS):
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $*.F | grep -v '^#' | grep -v '^$$' > $*.f
	$(FC) $(FFLAGS) -c $*.f

remap.o: remap.F
	@sed -e "s;MODEL_FLAGS1;$(MODFLAGS1);" remap.F | \
	sed -e "s;MODEL_FLAGS2;$(MODFLAGS2);" | \
	sed -e "s;MODEL_FLAGS3;$(MODFLAGS3);" | \
	sed -e "s;MODEL_FLAGS4;$(MODFLAGS4);" > foo.F
	$(CPP) $(CPPFLAGS) foo.F | grep -v '^#' | grep -v '^$$' > foo.f
	$(FC) $(FFLAGS)  -c foo.f
	@mv foo.o remap.o



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
