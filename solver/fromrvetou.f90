!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine converts the primitive variables to
! conservative variables (or vice versa)
! Written by Leung Shing Chi in 2016
! If you add your own variables in the WENO scheme,
! add your conversion step here.
! This subroutine takes in the U array and conversion mode
! Mode 0: From primitive to conservative
! Mode 1: From conservative to primitive
!
! Here is a reminder in how to add new physics:
! 1. Add your own module that contains the physicsb
! 2. Remind BuildWENO to include your quantity
! 3. Add the conversion here
! 4. Write a section in how to calculate the flux in Spatial
! 5. Add a flag to give signal to the program whenever you use the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FROMRVETOU
!$ACC ROUTINE (P_to_U) SEQ
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Real variables, geometric factor !
REAL*8 :: geom_bcell_x_l
REAL*8 :: geom_bcell_x_r
REAL*8 :: geom_bcell_y_l
REAL*8 :: geom_bcell_y_r
REAL*8 :: geom_bcell_z_l
REAL*8 :: geom_bcell_z_r

!------------------------------------------------------------------------------------------!
! Get cell centered magnetic field !

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = 0.5D0*(prim(ibx,j,k,l) + prim(ibx,j-1,k,l))
      bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
      bcell(ibz,j,k,l) = 0.5D0*(prim(ibz,j,k,l) + prim(ibz,j,k,l-1))
    END DO
  END DO
END DO
!$ACC END PARALLEL

!------------------------------------------------------------------------------------------!
! Convert primitive to conservative !

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Pass the primitive variables to the primitive to conservative warpper !
      CALL P_to_U(j,k,l)

	  END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert back to primitive variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FROMUTORVE
!$ACC ROUTINE (U_to_P) SEQ
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

! Real variables, geometric factor !
REAL*8 :: geom_bcell_x_l
REAL*8 :: geom_bcell_x_r
REAL*8 :: geom_bcell_y_l
REAL*8 :: geom_bcell_y_r
REAL*8 :: geom_bcell_z_l
REAL*8 :: geom_bcell_z_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For NM sectors !

!------------------------------------------------------------------------------------------!
! Get cell centered magnetic field !

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = 0.5D0*(cons(ibx,j,k,l) + cons(ibx,j-1,k,l))
      bcell(iby,j,k,l) = 0.5D0*(cons(iby,j,k,l) + cons(iby,j,k-1,l))
      bcell(ibz,j,k,l) = 0.5D0*(cons(ibz,j,k,l) + cons(ibz,j,k,l-1))
    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do the rest conversation !

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(present)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      
      ! Pass the primitive variables to the primitive to conservative warpper !
      CALL U_to_P(j,k,l)
      
    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE



