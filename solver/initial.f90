!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up everything you need for running simulations 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIAL_MODEL
USE DEFINITION
IMPLICIT NONE

! Include MPI only if MPI is used !
#ifdef MPI
include "mpif.h" 
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find MPI rank
#ifdef MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
#else
mpi_rank = 0
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First build up all database, EOS table and arrays !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'We shall now setup everything for the simulations'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for building arrays !

! Hydrodynamic arrays !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Build hydro arrays'
END IF
CALL BUILD_HYDRO
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Done building hydro arrays'
  WRITE(*,*)
END IF

! Grid face arrays !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Build grid variables'
END IF
CALL GETGRID
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Done building grid variables'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for setting initial conditions !

! Set initial conservative/primitive variables !
cons = 0.0D0
prim = 0.0D0

! Build initial models !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Now build the initial model'
END IF
CALL GET_MODEL
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Finished building initial model'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine do initial updates after setting up the inital model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initial_update
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
include "mpif.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find MPI rank
#ifdef MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
#else
mpi_rank = 0
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare everything ... !

! Check density !
CALL CUSTOM_CHECKRHO

! Convert from primitive to conservative variables !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Build conservative variables'
END IF
CALL FROMRVETOU
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Done building initial conservative variables'
  WRITE(*,*)
END IF

! set boundary conditions !
CALL BOUNDARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the items needed for SPATIAL

! Update pressure, etc ... !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Do Update'
END IF
CALL UPDATE (0)
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Done initial update'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print out !
IF(mpi_rank == 0) THEN
  WRITE(*,*) 'Finish initial...'
  WRITE(*,*)
END IF

END SUBROUTINE