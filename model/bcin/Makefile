#
#  $Id: Makefile.in,v 4.3 1998/07/07 20:42:34 dwd Exp $
#
# Makefile for bcin
#
# Model conditional compiliation flags
#
# sun        Sun compiler f77
# cray       Cray f77 environment
# cf90       Cray fortran90 compiler
# cmf        Thinking Machines Inc. connection machine compiler
#

CPP=cpp
FC=ftn
MODFLAGS=
MODCPP=
% LIBS= -L/users/skylling/local/lib -lnetcdf
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
FFLAGS= -u
LNFLAGS=
CPPFLAGS= $(MODFLAGS) $(MODCPP) -I/users/skylling/local/include -I..




SOURCES = lnblk.F parse.F varcde.F 

MAIN = bcin.F

OBJECTS= $(SOURCES:.F=.o)

.KEEP_STATE:

bcin: $(OBJECTS) bcin.o
	$(FC) $(FFLAGS) -o bcin bcin.o $(OBJECTS)



PHAKE:


clean:
	@rm -f .sb/Locked/\*
	@rm -f .sb/NewRoot/\*
	@rm -f .sb/OldRoot/\*
	@rm -f .sb/Refd/\*
	@rm -f .sb/Index
	rm -f *.o *~ *.fcm *.f90 *.l *.lst *.T *.pe.dp
	rm -f .make.state .nse_depinfo .make.depend*

