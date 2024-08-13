!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_EQN
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_HYDRO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Populate custom arrays to GPU !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_POPULATE
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Atmospheric values ... !
!$ACC UPDATE DEVICE(rho_floor, eps_floor) 

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_GRID
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom grid !
!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_BOUNDARY
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$ACC PARALLEL DEFAULT(PRESENT)
!$ACC LOOP GANG WORKER VECTOR COLLAPSE(3) 
DO l = -2, nz + 3
  DO k = -2, ny + 3
    DO j = 1, 3
      prim(ivx,1-j,k,l) = MIN(prim(ivx,1-j,k,l), 0.0D0)
      prim(ivx,nx+j,k,l) = MAX(prim(ivx,nx+j,k,l), 0.0D0)
    END DO
  END DO               
ENDDO

!$ACC LOOP GANG WORKER VECTOR COLLAPSE(3) 
DO l = 1, 3
  DO k = -2, ny + 3
    DO j = -2, nx + 3
      prim(ivz,j,k,1-l) = MIN(prim(ivz,j,k,1-l), 0.0d0)
      prim(ivz,j,k,nz+l) = MAX(prim(ivz,j,k,nz+l), 0.0d0)
    END DO
  END DO               
ENDDO 
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!
! Custom variable floor !
!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_CHECKRHO
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Check density and internal energy density !
      IF(prim(irho,j,k,l) < rho_floor) THEN
        prim(irho,j,k,l) = rho_floor
      END IF
      IF(eps(j,k,l) < eps_floor) THEN
        eps(j,k,l) = eps_floor
      ENDIF

    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building custom equations !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_SOURCE
!$ACC ROUTINE (GET_COORD) SEQ
!$ACC ROUTINE (COORD_DX) SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! real !
REAL*8 :: radius

! Threshold for atmosphere density
REAL*8 :: dphidr, dphidz
REAL*8 :: factor, diff

! Local !
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Add black hole gravity !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(radius, dphidr, dphidz, factor, diff)
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)
      CALL COORD_DX(j,k,l,dx_loc,dy_loc,dz_loc)
      diff = prim(irho,j,k,l) - rho_floor
      factor = MAX(SIGN(1.0D0, diff), 0.0D0)
      radius = DSQRT(x_loc**2 + z_loc**2)
      dphidr = x_loc/radius/((radius - r_sh)*(radius - r_sh))
      dphidz = z_loc/radius/((radius - r_sh)*(radius - r_sh))
      sc(ivx,j,k,l) = sc(ivx,j,k,l) + (-factor*prim(irho,j,k,l)*dphidr)
      sc(ivz,j,k,l) = sc(ivz,j,k,l) + (-factor*prim(irho,j,k,l)*dphidz)
      sc(itau,j,k,l) = sc(itau,j,k,l) + (-factor*prim(irho,j,k,l)*(prim(ivx,j,k,l)*dphidr + prim(ivz,j,k,l)*dphidz))
    END DO
  END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!
! Do custom updates !
!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUSTOM_UPDATE (p_in)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Integer !
INTEGER, INTENT (IN) :: p_in

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do custom operator split !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPERATOR_SPLIT
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print custom analysis file !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OPENFILE_CUSTOM
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

END SUBROUTINE
