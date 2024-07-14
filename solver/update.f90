!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This subroutine prepare the data necessary for constructing
! the flux for the spatial discretization.
! It takes input/output of the U array and the 
! mode p which decides whether or not to update
! the gravitational potential
! Written by Leung Shing Chi in 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE UPDATE (p_in)
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updates to hydrodynamic variables !

! Find pressure and speed of sound !
CALL FINDPRESSURE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom updates !

! Do custom updates !
CALL CUSTOM_UPDATE (p_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE