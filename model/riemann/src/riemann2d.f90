!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Riemann2D.f90, storing subroutines for generating 
! Initial conditions for riemann problem test in 2D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MHD Rotator test
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MHDROTOR(j_in,k_in,l_in,p_out,eps_out)
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
REAL*8 :: x, y, r

! Others !
REAL*8 :: fr, r1, r0

!*********************************************************************!

r0 = 0.1D0
r1 = 0.115D0

CALL GET_COORD(j_in,k_in,0,x,y,dummy)
r = DSQRT((x - 0.5D0)**2 + (y - 0.5D0)**2)
fr = (r1 - r)/(r1 - r0)
If(r <= r0) THEN
  p_out(irho) = 10.0D0
  p_out(ivx) = -1.0D0*(y - 0.5D0)/r0
  p_out(ivy) = 1.0D0*(x - 0.5D0)/r0
ELSEIF(r >= r1) THEN
  p_out(irho) = 1.0D0
  p_out(ivx) = 0.0D0
  p_out(ivy) = 0.0D0
ELSE
  p_out(irho) = 1.0D0 + 9.0d0*fr
  p_out(ivx) = -fr*1.0D0*(y - 0.5D0)/r
  p_out(ivy) = fr*1.0D0*(x - 0.5D0)/r
END IF
p_out(itau) = 0.5D0
p_out(ibx) = 2.5D0/DSQRT(4.0D0*pi)
eps_out = p_out(itau) / p_out(irho) / (ggas - 1.0D0)

!*********************************************************************!

END SUBROUTINE