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

! Normalization constant !
REAL*8, PARAMETER :: beta_norm = 100.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Problem setup !

! schwarzschild radius !
REAL*8, PARAMETER :: r_sh = 1.0

! angular velocity gradient !
REAL*8, PARAMETER :: q_grad = 2.0d0

! inner (equatorial) radius of the torus !
REAL*8, PARAMETER :: s_in = 3.0d0

! (equatorial) radius where the density is at maximum !
REAL*8, PARAMETER :: s_max = 4.7d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the SANE magnetic field setup 

! minimum density (fraction) to define the last contour of the vector potential !
REAL*8, PARAMETER :: rho_cut = 5.0d-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! atmospheric density !
REAL*8, PARAMETER :: rho_fac = 1.0d-4

! maximum density !
REAL*8, PARAMETER :: rho_max = 1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Atmospheric value !
REAL*8 :: rho_floor
REAL*8 :: eps_floor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Section for GPU !
#ifdef GPU
!$ACC declare create(rho_floor)
!$ACC declare create(eps_floor)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PARAMETER
