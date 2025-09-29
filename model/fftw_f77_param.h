c     This file contains PARAMETER statements for various constants
c     that can be passed to FFTW routines.  You should include
c     this file in any FORTRAN program that calls the fftw_f77
c     routines (either directly or with an #include statement
c     if you use the C preprocessor).

      INTEGER FFTW_FORWARD,FFTW_BACKWARD
      PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      INTEGER FFTW_IN_PLACE,FFTW_USE_WISDOM
      PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
