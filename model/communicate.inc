c
c  $Id: communicate.inc,v 1.3 1998/10/21 21:56:34 dwd Exp $
c
      integer my_pe, mastertask,ns_type_real,mpitag_nshift,mpitag_sshift
      integer mpitag_eshift,mpitag_wshift
c
      integer nbr_north,nbr_south,mpitag_io,nscont_type_real
      integer mpi_local_real,mpi_global_real,mpi_local_xz
      integer ns_type_double,mpi_local_double,mpi_global_double
      integer mpi_local_double_xz,nscont_type_double
      integer nst_type_real
      integer nbr_east,nbr_west
      integer iblkstart(nprocs),iblkcord(nprocs)
      integer jblkstart(nprocs),jblkcord(nprocs)
      integer iplan,jplan,rfplan,rbplan,comm
      logical iope
      common /communicate/ my_pe, mastertask, ns_type_real,
     &       mpitag_nshift,mpitag_sshift,jblkstart,
     &       nbr_north,nbr_south,mpitag_io,nscont_type_real,
     &       nbr_east,nbr_west,iblkstart,mpitag_eshift,
     &       mpitag_wshift,
     &       mpi_local_real,mpi_global_real,iplan,jplan,comm,
     &       rfplan,rbplan,mpi_local_xz,jblkcord,iblkcord,
     &       ns_type_double,mpi_local_double,mpi_global_double,
     &       mpi_local_double_xz,nscont_type_double,nst_type_real
      common /comlog/iope
