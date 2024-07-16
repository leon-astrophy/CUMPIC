!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assing MPI rank to GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ASSIGN_DEVICE
USE OPENACC
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
#ifdef GPU
include "mpif.h" 
!*******************************************************************!

! Get rank !
call MPI_COMM_RANK(new_comm, mpi_rank, ierror)

! Set device !
call acc_set_device_num(mpi_rank, acc_get_device_type())

!*******************************************************************!
#endif
#endif

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine handles data transfer between host and devices 
! AT THE BEGINNING of the simulations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE POPULATE_DEVICE
USE DEFINITION
IMPLICIT NONE

!-----------------------------------------------------------------------------------------------!

! Now populate all necessary, and reuseable arrays to the graphic cards !
!$ACC enter DATA COPYIN(bfac_xin, bfac_yin, bfac_zin, bfac_xout, bfac_yout, bfac_zout, & 
!$ACC prim, primL, primR, cons, consL, consR, flux, fluxL, fluxR, eps, l_rk, u_old, sc, &
!$ACC xf, yf, zf, bcell, ecell, ecorn, face_type, bcell_type, bface_type, scalar_type, neighbors, &
!$ACC starts, stops)

!-----------------------------------------------------------------------------------------------!

! Update device for variables defined in acc routine seq !
!$ACC UPDATE DEVICE(no_of_eq, imin, imax) 
!$ACC UPDATE DEVICE(irho, ivx, ivy, ivz, itau, ibx, iby, ibz)
!$ACC UPDATE DEVICE(xF, yF, zF)
!$ACC UPDATE DEVICE(primL, primR)
!$ACC UPDATE DEVICE(consL, consR)
!$ACC UPDATE DEVICE(prim, cons, eps, bcell)
!$ACC UPDATE DEVICE(fluxL, fluxR, flux, sc)

!-----------------------------------------------------------------------------------------------!

! Populate arrays according to model !
CALL CUSTOM_POPULATE

!-----------------------------------------------------------------------------------------------!

END SUBROUTINE
