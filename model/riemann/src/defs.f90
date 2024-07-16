!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Definition files for specific simulation models
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PARAMETER
IMPLICIT NONE
SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Adiabatic index !
REAL*8, PARAMETER :: ggas = 5.0d0/3.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define Riemann Problem Test !

! Define here !
INTEGER, PARAMETER :: toro_1 = 1
INTEGER, PARAMETER :: brio_wu = 2
INTEGER, PARAMETER :: mhd_rotor = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Which test model? 
INTEGER, PARAMETER :: riemann_test = mhd_rotor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PARAMETER
