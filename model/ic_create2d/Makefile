#
#  $Id: Makefile.in,v 4.7 1999/08/16 20:22:32 eds Exp $
#
# Makefile for ic_create
#
# Model conditional compiliation flags
#
#
# sun        Sun compiler f77
# cray       Cray f77 environment
# cf90       Cray fortran90 compiler
# cmf        Thinking Machines Inc. connection machine compiler
#

CPP=/usr/ccs/lib/cpp
FC=ifort
# FC=gfortran
MODFLAGS= -DSNGLSALT
MODCPP=
ROOT_DIR=/users/skylling/local
LIBS=  -lnetcdff
#
# Common sun f77 switches:
# -u         Undefined variables cause fatal error.
# -g         Debug.
# -fast      Fast execution. Optimize for speed.
# -parallel  Automatic and explicit loop parallelizations.
# -O4        Auto inline of functions.  Highest optimization.
# -sb        Produce table information for SourceBrowser
# -C         includes code for checking array bounds
# -Usym      Undefine preprocessor symbol
# -Dsym      Define preprocessor symbol
#
FFLAGS= -O3
LNFLAGS=
# FFLAGS= -u $(MODFLAGS) $(MODCPP) -I/home/vtmx2/include -I..
CPPFLAGS= $(MODFLAGS) $(MODCPP) -I$(ROOT_DIR)/include -I..

SOURCES = ies80.F

MAIN = ic_create.F

OBJECTS= $(SOURCES:.F=.o)

.KEEP_STATE:

ic_create: $(OBJECTS) ic_create.o
	$(FC)  -o ic_create ic_create.o $(OBJECTS) $(LIBS)




PHAKE:


clean:
	@rm -f .sb/Locked/\*
	@rm -f .sb/NewRoot/\*
	@rm -f .sb/OldRoot/\*
	@rm -f .sb/Refd/\*
	@rm -f .sb/Index
	rm -f *.o *~ *.fcm *.f90 *.l *.lst *.T *.pe.dp
	rm -f .make.state .nse_depinfo .make.depend*

