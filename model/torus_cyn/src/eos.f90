!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine update the pressure profile and their deriative    !
! once the density profile is being updated through rungekutta time  !
! evolution. It is being used in every time step, do not confused it !
! with subroutine GETRHOEOSRTOP                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_PRESSURE (j_in, k_in, l_in, rho_in, eps_in, p_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input !
INTEGER, INTENT (IN) :: j_in, k_in, l_in

! Input density !
REAL*8, INTENT (IN) :: rho_in, eps_in

! Output value !
REAL*8, INTENT (OUT) :: p_out

!-----------------------------------------------------!

! find pressure !
p_out = rho_in*eps_in*(ggas-1.0d0)

!------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the speed of sound !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_SOUNDSPEED (p_in, rho_in, eps_in, cs_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in

! Output value !
REAL*8, INTENT (OUT) :: cs_out

!-----------------------------------------------------!

! We do the DM case first !
cs_out = DSQRT(ggas*p_in/rho_in)

!-----------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_EPSILON (p_in, rho_in, eps_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input density !
REAL*8, INTENT (IN) :: p_in, rho_in

! Output value ! 
REAL*8, INTENT (OUT) :: eps_out

!-----------------------------------------------------!

! For DM Output !
eps_out = p_in/rho_in/(ggas - 1.0D0)

!-----------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE