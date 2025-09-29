      subroutine imag2p(phi)
c
c    $Id: imag.F,v 1.6 1998/10/21 21:56:38 dwd Exp $
c
c this routine sets the lateral image points
c
#if defined(USEMPI)
      include 'mpif.h'
#endif
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "cyclic.inc"
#include "olemtke.inc"
#include "moddef.inc"
#include "communicate.inc"
#if defined(ICE)
#include "ice.inc"
#endif
#include "finitv.inc"
c
      real phi(0:idim1,0:jdim1,0:kdim1)

      
      if((blockstart(my_pe+1)+jdim).ge.jtotal) then
        j = jdim1
        do i=0,idim1
          do k=0,kdim1
c
c wall or zero gradient boundary condition
c
            phi(i,j,k) = phi(i,j-1,k)
          enddo
        enddo
      endif
c
      call MPI_BARRIER(comm,ierr)

c
      if(blockstart(my_pe+1).eq.1) then
c 
        do i=0,idim1
          do k=0,kdim1
c
c do wall bc
c
           phi(i,0,k) = phi(i,1,k)

          enddo
        enddo
c
      endif
      call MPI_BARRIER(comm,ierr)
      return
      end
