!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! mpi.f90
! contains all essential functions for mpi parallelization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------
!
! This subroutine initialize MPI parallelization, reference: iharm3d
!
!----------------------------------------------------------------------

#ifdef MPI 
!***************************************************************************************************************************!

SUBROUTINE SETUP_MPI
USE DEFINITION
IMPLICIT NONE

include "mpif.h"

! Integer !
INTEGER :: i, j, k, l
INTEGER :: blocks, stride
INTEGER :: reduced_nvars

! For cartesian topology !
INTEGER, DIMENSION(1:3) :: dims
INTEGER, DIMENSION(1:3) :: n_ind
INTEGER, DIMENSION(1:3) :: coords
LOGICAL, DIMENSION(1:3) :: periods

! For communication between blocks !
INTEGER, DIMENSION(1:4) :: sizes
INTEGER, DIMENSION(1:4) :: subsizes
INTEGER, DIMENSION(1:4) :: starting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if the number of processes is consistent !

call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierror)
numcpus = NXCPU*NYCPU*NZCPU
IF(mpi_size .ne. numcpus) THEN
  IF(mpi_rank == 0) THEN
    write (*,*) 'Error, the number of MPI processes called is not equal to the total number of CPU requested'
  END IF
  STOP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check CPU number divisibility !

IF(MOD(NXTOT,NXCPU) .ne. 0) THEN
  IF(mpi_rank == 0) THEN 
    write (*,*) 'Error, number of grid along the x-direction is not divisible by the number of CPU'
  END IF
  STOP
END IF
IF(MOD(NYTOT,NYCPU) .ne. 0) THEN
  IF(mpi_rank == 0) THEN 
    write (*,*) 'Error, number of grid along the y-direction is not divisible by the number of CPU'
  END IF
  STOP
END IF
IF(MOD(NZTOT,NZCPU) .ne. 0) THEN
  IF(mpi_rank == 0) THEN 
    write (*,*) 'Error, number of grid along the z-direction is not divisible by the number of CPU'
  END IF
  STOP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create cartesian topology !

! Dimension along each CPU direction !
dims (1) = NZCPU
dims (2) = NYCPU
dims (3) = NXCPU

! Indicate if we need periodic boundary conditions !
IF(boundary_flag(1) == 0 .or. boundary_flag(2) == 0) THEN
  periods(3) = .true.
ELSE
  periods(3) = .false.
END IF
IF(boundary_flag(3) == 0 .or. boundary_flag(4) == 0) THEN
  periods(2) = .true.
ELSE
  periods(2) = .false.
END IF
IF(boundary_flag(5) == 0 .or. boundary_flag(6) == 0) THEN
  periods(1) = .true.
ELSE
  periods(1) = .false.
END IF

! Create cartesian topology !
CALL MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, new_comm, ierror)

! Find location within the cartesian topology !
! Note, the "cartesian coordinate is in the z, y, x convention !
CALL MPI_Comm_rank(new_comm, mpi_rank, ierror)
CALL MPI_Cart_coords(new_comm, mpi_rank, 3, coords, ierror)

! Find the ranks of neighbors, including edge/corner neighbors, centered around (1,1) !
! Note: Assign MPI_PROC_NULL to processes that is not used. Example topology !
!
!    |-1 |-1 |
!------------------
! -1 | 2 | 3 | -1
!------------------
! -1 | 0 | 1 | -1
!------------------
!    |-1 |-1 |
!
DO l = -1, 1
  n_ind(1) = coords(1) + l
  DO k = -1, 1
    n_ind(2) = coords(2) + k
    DO j = -1, 1
      n_ind(3) = coords(3) + j
        IF(((n_ind(1) < 0 .or. n_ind(1) >= NZCPU) .and. .not.(periods(1))) .or. &
           ((n_ind(2) < 0 .or. n_ind(2) >= NYCPU) .and. .not.(periods(2))) .or. &
           ((n_ind(3) < 0 .or. n_ind(3) >= NXCPU) .and. .not.(periods(3)))) THEN
          neighbors(l+1,k+1,j+1) = MPI_PROC_NULL 
        ELSE
          CALL MPI_Cart_rank(new_comm, n_ind, neighbors(l+1,k+1,j+1), ierror)
        END IF
    END DO
  END DO
END DO

! Assign global start and ending index !
! Note the reversal !
starts(1) = coords(3)*nx
starts(2) = coords(2)*ny
starts(3) = coords(1)*nz
stops(1) = starts(1) + nx
stops(2) = starts(2) + ny
stops(3) = starts(3) + nz

! Reduced number of variables, do not pass face-centered field to ghost zone !
reduced_nvars = ibx - imin

! Now create MPI datatype for message passing between ghost cells !
! This is for along the NZ direction !
sizes = (/no_of_eq,NX+2*NGHOST,NY+2*NGHOST,NZ+2*NGHOST/)
subsizes = (/reduced_nvars,NX+2*NGHOST,NY+2*NGHOST,NGHOST/)
starting = (/0,0,0,0/)
CALL MPI_Type_create_subarray(4, sizes, subsizes, starting, MPI_ORDER_FORTRAN, MPI_DOUBLE, face_type(1), ierror)
CALL MPI_Type_commit(face_type(1), ierror)

! Along the NY direction !
sizes = (/no_of_eq,NX+2*NGHOST,NY+2*NGHOST,NZ+2*NGHOST/)
subsizes = (/reduced_nvars,NX+2*NGHOST,NGHOST,NZ/)
starting = (/0,0,0,0/)
CALL MPI_Type_create_subarray(4, sizes, subsizes, starting, MPI_ORDER_FORTRAN, MPI_DOUBLE, face_type(2), ierror)
CALL MPI_Type_commit(face_type(2), ierror)

! Along the NX direction !
sizes = (/no_of_eq,NX+2*NGHOST,NY+2*NGHOST,NZ+2*NGHOST/)
subsizes = (/reduced_nvars,NGHOST,NY,NZ/)
starting = (/0,0,0,0/)
CALL MPI_Type_create_subarray(4, sizes, subsizes, starting, MPI_ORDER_FORTRAN, MPI_DOUBLE, face_type(3), ierror)
CALL MPI_Type_commit(face_type(3), ierror)

! Barrier !
CALL MPI_BARRIER(new_comm, ierror) 

END SUBROUTINE

!----------------------------------------------------------------------
!
! This subroutine finialize MPI parallelization
!
!----------------------------------------------------------------------

SUBROUTINE FINAL_MPI
USE DEFINITION
IMPLICIT NONE
include "mpif.h"
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Finialize !
call MPI_FINALIZE(ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!----------------------------------------------------------------------
!
! This subroutine finialize transfer message from active cell to ghost
! cell across MPI processes along the x-direction
!
!----------------------------------------------------------------------

SUBROUTINE MPI_BOUNDARYP_X
USE DEFINITION
IMPLICIT NONE
include "mpif.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Send my last NG cell to my right neighbour's ghost cell !
! Receive my left neighbour's last NG cell to my ghost cell !
CALL MPI_Sendrecv(prim(1,NX-NGHOST+1,1,1), 1, face_type(3), neighbors(1,1,2), 0, &
                  prim(1,1-NGHOST,1,1), 1, face_type(3), neighbors(1,1,0), 0, & 
                  new_comm, MPI_STATUS_IGNORE, ierror)


! Send my first NG cell to my left neighbour's ghost cell !
! Receive my right neighbour's first NG cell to my ghost cell !
CALL MPI_Sendrecv(prim(1,1,1,1), 1, face_type(3), neighbors(1,1,0), 0, &
                  prim(1,NX+1,1,1), 1, face_type(3), neighbors(1,1,2), 0, & 
                  new_comm, MPI_STATUS_IGNORE, ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!***************************************************************************************************************************!
#endif

! End of file !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!