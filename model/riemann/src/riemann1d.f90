!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Riemann1D.f90, storing subroutines for generating 
! Initial conditions for riemann problem test in 1D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Toro Rieamnn problem test 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TORO_SHOCKTUBE_1(j_in,k_in,l_in,p_out,eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in,k_in,l_in

! Output !
REAL*8, INTENT(OUT) :: eps_out
REAL*8, INTENT(OUT), DIMENSION(1:no_of_eq) :: p_out

! Dummy !
REAL*8 :: dummy 

! Coordinates !
REAL*8 :: x

!*********************************************************************!

CALL GET_COORD(j_in,0,0,x,dummy,dummy)
if(x <= 0.5D0) THEN
	p_out(irho) = 1.0D0
	p_out(itau) = 1.0D0
else
	p_out(irho) = 0.125D0
	p_out(itau) = 0.1D0
endif
eps_out = p_out(itau)/p_out(irho)/(ggas - 1.0d0)

!*********************************************************************!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Brio-Wu Shock tube test 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BRIOWU_SHOCKTUBE(j_in,k_in,l_in,p_out,eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in,k_in,l_in

! Output !
REAL*8, INTENT(OUT) :: eps_out
REAL*8, INTENT(OUT), DIMENSION(1:no_of_eq) :: p_out

! Dummy !
REAL*8 :: dummy 

! Coordinates !
REAL*8 :: x

!*********************************************************************!

CALL GET_COORD(j_in,0,0,x,dummy,dummy)
if(x <= 0.5D0) THEN
	p_out(irho) = 1.0D0
	p_out(itau) = 1.0D0
	p_out(iby) = 1.0D0
else
	p_out(irho) = 0.125D0
	p_out(itau) = 0.1D0
	p_out(iby) = -1.0D0
endif
p_out(ibx) = 0.75D0
eps_out = p_out(itau)/p_out(irho)/(ggas - 1.0d0)

!*********************************************************************!

END SUBROUTINE
