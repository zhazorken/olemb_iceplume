c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c     begin file communicate.h
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


      integer mpitag_io
      parameter(mpitag_io = 5)

      integer nbr_west, nbr_east, nbr_north, nbr_south
      integer block_start(nprocs)
      integer comm_2d, mpi_ew_type,
     &                 mpi_interior_int,
     &                 mpi_block_real,
     &                 mpi_interior_dble,
     &                 mpi_interior_int_global,
     &                 mpi_global_real,
     &                 mpi_interior_dble_global

      common /communicate/ nbr_west, nbr_east, nbr_north, nbr_south,
     &                     comm_2d, block_start,
     &                     mpi_ew_type,
     &                     mpi_interior_int,
     &                     mpi_block_real,
     &                     mpi_interior_dble,
     &                     mpi_interior_int_global,
     &                     mpi_global_real,
     &                     mpi_interior_dble_global

c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c     end file communicate.h
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
