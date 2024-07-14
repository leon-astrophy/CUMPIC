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
      CALL P_to_U(prim(:,j,k,l), bcell(:,j,k,l), eps(j,k,l), cons(:,j,k,l))

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
      CALL U_to_P(bcell(:,j,k,l), cons(:,j,k,l), prim(:,j,k,l), eps(j,k,l))
      
    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given primitive variables, construct conservative variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_to_U(prim_in, bcell_in, eps_in, cons_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL*8, INTENT (IN) :: eps_in
REAL*8, INTENT (IN), DIMENSION (ibx:ibz) :: bcell_in
REAL*8, INTENT (IN), DIMENSION (1:no_of_eq) :: prim_in
REAL*8, INTENT (OUT), DIMENSION (1:no_of_eq) :: cons_out

! Local array !
REAL*8 :: vsquare, bsquare

!--------------------------------------------------------------------------!
! Do the conversion 

! Get v*v and b*b !
vsquare = dot_product(prim_in(ivx:ivz), prim_in(ivx:ivz))
bsquare = dot_product(bcell_in(ibx:ibz), bcell_in(ibx:ibz))	  

! Assign conservative variables to core hyrodynamic variables !
cons_out(irho) = prim_in(irho)
cons_out(ivx:ivz) = prim_in(ivx:ivz)*prim_in(irho)
cons_out(itau) = prim_in(irho)*(eps_in + 0.5D0*vsquare) + 0.5D0*bsquare

! Magnetic field at face center !
cons_out(ibx:ibz) = prim_in(ibx:ibz)

! For any scalar variables 
IF(itau+1 .ne. ibx) THEN
  cons_out(itau+1:ibx-1) = prim_in(itau+1:ibx-1)*prim_in(irho)
END IF

!--------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given conservative variables, invert to get primitive variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE U_to_P(bcell_in, cons_in, prim_out, eps_out)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input
REAL*8, INTENT (IN), DIMENSION (ibx:ibz) :: bcell_in
REAL*8, INTENT (IN), DIMENSION (1:no_of_eq) :: cons_in

! Output
REAL*8, INTENT (OUT) :: eps_out
REAL*8, INTENT (OUT), DIMENSION (1:no_of_eq) :: prim_out

! Local array !
REAL*8 :: vsquare, bsquare

!--------------------------------------------------------------------------!
! Do the conversion 

! Core primitive variables !
prim_out(irho) = cons_in(irho)
prim_out(ivx:ivz) = cons_in(ivx:ivz)/cons_in(irho)

! Magnetic fields !
prim_out(ibx:ibz) = cons_in(ibx:ibz)

! Get v*v and b*b !
bsquare = dot_product(bcell_in(ibx:ibz), bcell_in(ibx:ibz))
vsquare = dot_product(prim_out(ivx:ivz), prim_out(ivx:ivz))

! For any scalar variables 
IF(itau+1 .ne. ibx) THEN
  prim_out(itau+1:ibx-1) = cons_in(itau+1:ibx-1)/cons_in(irho)
END IF

! epsilon here, CAUTION: no negativity check !
eps_out = (cons_in(itau) - 0.5D0*bsquare)/cons_in(irho) - 0.5D0 * vsquare

!--------------------------------------------------------------------------!

END SUBROUTINE



