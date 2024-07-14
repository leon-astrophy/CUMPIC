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
REAL*8, PARAMETER :: ggas = 2.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define Riemann Problem Test !

! Define here !
INTEGER, PARAMETER :: brio_wu = 1
INTEGER, PARAMETER :: mhd_rotor = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Which test model? 
INTEGER, PARAMETER :: riemann_test = brio_wu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PARAMETER
