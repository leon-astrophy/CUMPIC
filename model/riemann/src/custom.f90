!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      eps(j,k,l) = MAX(eps(j,k,l), 0.0d0)
    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do custom operator split !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print custom analysis file !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE
