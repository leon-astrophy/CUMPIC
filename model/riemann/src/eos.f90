!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDPRESSURE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

!--------------------------------------------------------------------------------------

! The following steps are more or less similar , so no repeat 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 1 - NGHOST, nz + NGHOST
  DO k = 1 - NGHOST, ny + NGHOST
    DO j = 1 - NGHOST, nx + NGHOST
      prim(itau,j,k,l) = prim(irho,j,k,l)*eps(j,k,l)*(ggas - 1.0D0) 
    END DO
  END DO
END DO
!$ACC END PARALLEL

!--------------------------------------------------------------------------------------

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_SOUNDSPEED (p_in, rho_in, cs_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: p_in, rho_in

! Output value !
REAL*8, INTENT (OUT) :: cs_out

! We do the DM case first !
cs_out = DSQRT(ggas*p_in/rho_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_EPSILON (rho_in, p_in, eps_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: rho_in, p_in

! Output value ! 
REAL*8, INTENT (OUT) :: eps_out

! For DM Output !
eps_out = p_in/rho_in/(ggas - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE