!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output the hydrodynamic variable profiles
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE print_hydroprofile
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Characeter
character(len=99) :: n_file
character(len=99) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for GPU !

#ifdef GPU
!$ACC UPDATE HOST(prim2(imin2:imax2,:,:,:), epsilon2, bcell)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open file 

! write to character !
write(n_file,'(I10)') n_iter

! assign !
filename = './outfile/star_weno_data_'// trim(adjustl(n_file)) //'.hdf5'

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print important variables !

CALL hdf5_print_1darray('xf', xf(0:NX), NX+1)
CALL hdf5_print_1darray('yf', yf(0:NY), NY+1)
CALL hdf5_print_4darray('prim', prim(:,1:NX,1:NY,1:NZ), no_of_eq, NX, NY, NZ)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE 

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

! define DIMENSION !
space_rank = 1
dist_dims(1) = 1

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vars,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!--------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output 1D array to the hdf5 file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_1darray(name, vars, ngrid)
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Input !
character(len=*), INTENT(IN) :: name
integer, INTENT(IN) :: ngrid
real*8, intent(in) :: vars(ngrid)

! integer !
integer :: space_rank

! define DIMENSION !
space_rank = 1
dist_dims(1) = ngrid

! open dataspace !
call h5screate_simple_f(space_rank,dist_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vars,dist_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!--------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output 3D array to the hdf5 file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_3darray(name, vars, nx_in, ny_in, nz_in)
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Input !
character(len=*), INTENT(IN) :: name
integer, INTENT(IN) :: nx_in, ny_in, nz_in
real*8, intent(in) :: vars(nx_in,ny_in,nz_in)

! integer !
integer :: space_rank

! define DIMENSION !
space_rank = 3
eps_dims(1) = nx_in
eps_dims(2) = ny_in
eps_dims(3) = nz_in

! open dataspace !
call h5screate_simple_f(space_rank,eps_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vars,eps_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!--------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine output 4D array to the hdf5 file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE hdf5_print_4darray(name, vars, n0_in, nx_in, ny_in, nz_in)
USE HDF5
USE DEFINITION
IMPLICIT NONE

! Input !
character(len=*), INTENT(IN) :: name
integer, INTENT(IN) :: n0_in, nx_in, ny_in, nz_in
real*8, intent(in) :: vars(n0_in,nx_in,ny_in,nz_in)

! integer !
integer :: space_rank

! define DIMENSION !
space_rank = 4
data_dims(1) = n0_in
data_dims(2) = nx_in
data_dims(3) = ny_in
data_dims(4) = nz_in

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vars,data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!--------------------------------------------------------------------------!

END SUBROUTINE