C Initialize U and V with depth profiles
C read in from netCDF file with header of teh form
C
C      netcdf ??????? {
C       dimensions:
C        Z = ## ;
C       variables:
C        float Z(Z) ;
C        float U(Z) ;
C        float V(Z) ;
C      }
C
C
      subroutine readUV(FILE_NAME,U,V,deltaZ)
       include 'param.inc'
       real, intent(in) :: DeltaZ
       character*(*), intent(in) :: FILE_NAME
       real, dimension(kdim), intent(inout) :: U, V
       real :: Zd,Zi,Ze
       real, dimension(2) :: Uin, Vin, Zin, Zerrs
       integer, dimension(2) :: idxclose
       integer :: NVAR
       parameter (NVAR=3)
       integer :: varids(NVAR)
       integer :: retval, dimid, ncid, NZ
       integer :: i,j,k,l,m
C       Reading U,V data from FILE_NAME
C       Open the data file
       retval = nf_open(FILE_NAME,0,ncid)
       call handle_err(retval)
C       File open
C       Get the dimension ID for Z
       retval = nf_inq_dimid(ncid,'Z',dimid)
       call handle_err(retval)
C       Get length of Z
       retval = nf_inq_dimlen(ncid,dimid,NZ)
       call handle_err(retval)
C       Get the variable IDs for Z, U and V
       retval = nf_inq_varid(ncid,'Z',varids(1))
       call handle_err(retval)

       retval = nf_inq_varid(ncid,'U',varids(2))
       call handle_err(retval)

       retval = nf_inq_varid(ncid,'V',varids(3))
       call handle_err(retval)

C
       do i=1,kdim
        Zi = (0.5-i)*DeltaZ
C        print *,Zi
C        Find the closest pair of points to Zi
        retval =NF_GET_VAR1_real(NCID, varids(1), 1, Zd)
        call handle_err(retval)

        Zin(1) = Zd
        Zerrs(1) = (Zin(1)-Zi)
        idxclose(1) = 1
        Zin(2) = Zd
        Zerrs(2) = (Zin(2)-Zi)
        idxclose(2) = 1


        do j=2,NZ
C         Ensure that |Zerrs(1)| is the larger
         if(abs(Zerrs(1)).lt.abs(Zerrs(2))) then
C         Copy into temp variables
          Ze = Zerrs(1)
          Zd = Zin(1)
          k = idxclose(1)
C          Move larger value into 1st position
          Zerrs(1) = Zerrs(2)
          Zin(1) = Zin(2)
          idxclose(1) = idxclose(2)
C          Move smaller value into 2st position
          Zerrs(2) = Ze
          Zin(2) = Zd
          idxclose(2) = k
         end if

         retval =NF_GET_VAR1_real(NCID, varids(1),j,Zd)
         call handle_err(retval)

         Ze = (Zd-Zi)

         if(abs(Ze).le.abs(Zerrs(1))) then
          Zerrs(1) = Ze
          Zin(1) = Zd
          idxclose(1) = j
        end if
C         Ensure that |Zerrs(1)| is the larger
         if(abs(Zerrs(1)).lt.abs(Zerrs(2))) then
C         Copy into temp variables
          Ze = Zerrs(1)
          Zd = Zin(1)
          k = idxclose(1)
C          Move larger value into 1st position
          Zerrs(1) = Zerrs(2)
          Zin(1) = Zin(2)
          idxclose(1) = idxclose(2)
C          Move smaller value into 2st position
          Zerrs(2) = Ze
          Zin(2) = Zd
          idxclose(2) = k
         end if


        end do
C        Interpolate
C        For cases not interior to the given domain
C        Use the value at the closest point
        if(Zerrs(1)*Zerrs(2).ge.0) then
         retval =NF_GET_VAR1_real (NCID, varids(2),idxclose(2),U(i))
         call handle_err(retval)
         retval =NF_GET_VAR1_real (NCID, varids(3),idxclose(2),V(i))
         call handle_err(retval)
C
C       For cases interior to the given domain
C       Use simple linear interpolation between the closes 2 points
        else
         retval =NF_GET_VAR1_real (NCID, varids(2),idxclose(1),Uin(1))
         call handle_err(retval)
         retval =NF_GET_VAR1_real (NCID, varids(3),idxclose(1),Vin(1))
         call handle_err(retval)
         retval =NF_GET_VAR1_real (NCID, varids(2),idxclose(2),Uin(2))
         call handle_err(retval)
         retval =NF_GET_VAR1_real (NCID, varids(3),idxclose(2),Vin(2))
         call handle_err(retval)
         U(i) = (Uin(1)*Zerrs(2)-Uin(2)*Zerrs(1))/(Zerrs(2)-Zerrs(1))
         V(i) = (Vin(1)*Zerrs(2)-Vin(2)*Zerrs(1))/(Zerrs(2)-Zerrs(1))
        end if
C        print *,i,Zi,Zin(1),Zin(2),U(i),V(i)
       end do
C       Close the data file
       retval = nf_close(ncid)
      end subroutine

      subroutine handle_err(errcode)
       implicit none
       include 'netcdf.inc'
       integer errcode
       if (errcode .ne. nf_noerr) then
        print *, 'Error: ', nf_strerror(errcode)
        stop 2
       end if
      end subroutine
