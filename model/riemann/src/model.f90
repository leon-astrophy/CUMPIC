!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! integer !
INTEGER :: i, j, k, l

!---------------------------------------------------------------------------!

! Pass to wrapper !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      CALL RIEMANN_PROBLEM(j,k,l,prim(:,j,k,l),eps(j,k,l))
    END DO
  END DO
END DO

!---------------------------------------------------------------------------!

! Boundary condition !
CALL BOUNDARYP

!---------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interface for calling Riemann problem 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RIEMANN_PROBLEM(j_in,k_in,l_in,p_out,eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in,k_in,l_in

! Output !
REAL*8, INTENT(OUT) :: eps_out
REAL*8, INTENT(OUT), DIMENSION(1:no_of_eq) :: p_out

! Choose by case !
SELECT CASE (riemann_test)
!---------------------------------------------------------------------------!
CASE (brio_wu)
  CALL BRIOWU_SHOCKTUBE(j_in,k_in,l_in,p_out,eps_out)
CASE (mhd_rotor)
  CALL MHDROTOR(j_in,k_in,l_in,p_out,eps_out)

!---------------------------------------------------------------------------!
END SELECT

END SUBROUTINE


