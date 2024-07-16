!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
#ifdef USEHF5
USE HDF5
#endif
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
include "mpif.h"
#endif

! Characeter
character(len=99) :: n_file
character(len=99) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for GPU !

#ifdef GPU
!$ACC UPDATE HOST(prim, eps)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI stuff !

#ifdef MPI
CALL MPI_COMM_SIZE(new_comm, mpi_size, ierror)
CALL MPI_COMM_RANK(new_comm, mpi_rank, ierror)
#endif 

#ifdef USEHF5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open file 

! write to character !
write(n_file,'(I10)') n_iter

! assign !
filename = './outfile/star_weno_data_'// trim(adjustl(n_file)) //'.hdf5'

! create interface !
call h5open_f(error)

!---------------------------------------------------------------------------------------------------!

! open the file !
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
#ifdef MPI
CALL h5pset_fapl_mpio_f(plist_id, new_comm, MPI_INFO_NULL, error)
#endif

! create the file #
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print grid !
IF(NXTOT > 1) CALL hdf5_print_grid(x_dir)
IF(NYTOT > 1) CALL hdf5_print_grid(y_dir)
IF(NZTOT > 1) CALL hdf5_print_grid(z_dir)

! Primitive variables !
CALL hdf5_print_prim

! Epsilon !
CALL hdf5_print_eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close file !
CALL h5fclose_f(file_id, error)

! close the file !
CALL h5pclose_f(plist_id, error)

!---------------------------------------------------------------------------------------------------!

! close interface !
CALL h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

END SUBROUTINE 

#ifdef USEHF5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output scalar variables to the hdf5 file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_scalar(name, vars)
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Input !
character(len=*), INTENT(IN) :: name
real*8, INTENT(IN) :: vars

! integer !
integer :: space_rank

! data dimension !
integer(HSIZE_T) :: dims(1)

! define DIMENSION !
space_rank = 1
dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dims,space_id,error)

! create dataset !
call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,space_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vars,dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(space_id,error)

!--------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output grids to the hdf5 file 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_grid(dir_in)
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: dir_in

! integer !
integer :: space_rank

! data dimension !
integer(HSIZE_T) :: dims(1)

! Count and offset !
INTEGER(HSIZE_T), DIMENSION(1) :: count
INTEGER(HSSIZE_T), DIMENSION(1) :: offset 

! local !
character(len=99) :: name

!------------------------------------------------------------------------------------!

! Set the rank of the data set !
space_rank = 1

SELECT CASE (dir_in)
!*****************************************!
! Set dimensions, count and off set !
case(x_dir)
  name = "xf"
  dims(1) = nxtot+1
  count = (/nx+1/)
  offset = (/starts(1)/)
case(y_dir)
  name = "yf"
  dims(1) = nytot+1
  count = (/ny+1/)
  offset = (/starts(2)/)
case(z_dir)
  name = "zf"
  dims(1) = nztot+1
  count = (/nz+1/)
  offset = (/starts(3)/)
!*****************************************!
end select 

! Create data space !
CALL h5screate_simple_f(space_rank, dims, space_id, error)

! Create data set !
CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, error)

! Create memory space !
CALL h5screate_simple_f(space_rank, count, mem_id, error)

! Select hyperslab in the file
CALL h5dget_space_f(dset_id, space_id, error)
CALL h5sselect_hyperslab_f (space_id, H5S_SELECT_SET_F, offset, count, error) 

SELECT CASE (dir_in)
!*****************************************!
! Write to hyper slab !
case(x_dir)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xf(0:nx), count, error, file_space_id = space_id, mem_space_id = mem_id)
case(y_dir)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yf(0:ny), count, error, file_space_id = space_id, mem_space_id = mem_id)
case(z_dir)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, zf(0:nz), count, error, file_space_id = space_id, mem_space_id = mem_id)
!*****************************************!
end select

! Close all stuff
CALL h5dclose_f(dset_id, error) 
CALL h5sclose_f(space_id, error)
CALL h5sclose_f(mem_id, error)

!------------------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output primitive variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_eps
USE HDF5
USE DEFINITION
IMPLICIT NONE

! integer !
integer :: space_rank

! data dimension !
integer(HSIZE_T) :: dims(3)

! Count and offset !
INTEGER(HSIZE_T), DIMENSION(3) :: count
INTEGER(HSSIZE_T), DIMENSION(3) :: offset 

! local !
character(len=99) :: name

!------------------------------------------------------------------------------------!

! Set the rank of the data set !
space_rank = 3

! Set dimensions, count and off set !
name = "eps"
dims(1) = nxtot
dims(2) = nytot
dims(3) = nztot
count = (/nx,ny,nz/)
offset = (/starts(1),starts(2),starts(3)/) 

! Create data space !
CALL h5screate_simple_f(space_rank, dims, space_id, error)

! Create data set !
CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, error)

! Create memory space !
CALL h5screate_simple_f(space_rank, count, mem_id, error)

! Select hyperslab in the file
CALL h5dget_space_f(dset_id, space_id, error)
CALL h5sselect_hyperslab_f (space_id, H5S_SELECT_SET_F, offset, count, error) 

!*****************************************!
! Write to hyper slab !
CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eps(1:nx,1:ny,1:nz), count, error, file_space_id = space_id, mem_space_id = mem_id)

! Close all stuff
CALL h5dclose_f(dset_id, error) 
CALL h5sclose_f(space_id, error)
CALL h5sclose_f(mem_id, error)

!------------------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output primitive variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_prim
USE HDF5
USE DEFINITION
IMPLICIT NONE

! integer !
integer :: space_rank

! data dimension !
integer(HSIZE_T) :: dims(4)

! Count and offset !
INTEGER(HSIZE_T), DIMENSION(4) :: count
INTEGER(HSSIZE_T), DIMENSION(4) :: offset 

! local !
character(len=99) :: name

!------------------------------------------------------------------------------------!

! Set the rank of the data set !
space_rank = 4

! Set dimensions, count and off set !
name = "prim"
dims(1) = no_of_eq
dims(2) = nxtot
dims(3) = nytot
dims(4) = nztot
count = (/no_of_eq,nx,ny,nz/)
offset = (/0,starts(1),starts(2),starts(3)/) 

! Create data space !
CALL h5screate_simple_f(space_rank, dims, space_id, error)

! Create data set !
CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, space_id, dset_id, error)

! Create memory space !
CALL h5screate_simple_f(space_rank, count, mem_id, error)

! Select hyperslab in the file
CALL h5dget_space_f(dset_id, space_id, error)
CALL h5sselect_hyperslab_f (space_id, H5S_SELECT_SET_F, offset, count, error) 

!*****************************************!
! Write to hyper slab !
CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, prim(:,1:nx,1:ny,1:nz), count, error, file_space_id = space_id, mem_space_id = mem_id)

! Close all stuff
CALL h5dclose_f(dset_id, error) 
CALL h5sclose_f(space_id, error)
CALL h5sclose_f(mem_id, error)

!------------------------------------------------------------------------------------!

END SUBROUTINE
#endif