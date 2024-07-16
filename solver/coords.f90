!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine sets up the face position of each grid, assuming constant
! step size dx/dy/dz 
!
! Written by Leung Shing Chi in 2016
! Rewritten by H.S. Leon Chan in 2024
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETGRID
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid position of motherboard

! Then, get interface coordinate !
DO j = - NGHOST, nx + NGHOST
	xF(j) = x_start + (starts(1) + DBLE(j))*dx
END DO
DO k = - NGHOST, ny + NGHOST
	yF(k) = y_start + (starts(2) + DBLE(k))*dy
END DO
DO l = - NGHOST, nz + NGHOST
	zF(l) = z_start + (starts(3) + DBLE(l))*dz
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Or, provide your own grid function if you dislike

CALL CUSTOM_GRID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the cell-center coordinate
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GET_COORD(j_in,k_in,l_in,x_out,y_out,z_out)
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Output !
REAL*8, INTENT(OUT) :: x_out, y_out, z_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute !
x_out = 0.5D0*(xF(j_in) + xF(j_in-1))
y_out = 0.5D0*(yF(k_in) + yF(k_in-1))
z_out = 0.5D0*(zF(l_in) + zF(l_in-1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the cell-center coordinate
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE COORD_DX(j_in,k_in,l_in,dx_out,dy_out,dz_out)
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Output !
REAL*8, INTENT(OUT) :: dx_out, dy_out, dz_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute !
dx_out = (xF(j_in) - xF(j_in-1))
dy_out = (yF(k_in) - yF(k_in-1))
dz_out = (zF(l_in) - zF(l_in-1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the geometric source terms 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GEOM_SOURCE(j_in,k_in,l_in)
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the geometric factors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GEOM_FLUX(dir_in,j_in,k_in,l_in,geom_flux_p,geom_flux_c,geom_flux_m)
!$ACC ROUTINE (COORD_DX) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: dir_in, j_in, k_in, l_in

! Output !
REAL*8, INTENT (OUT) :: geom_flux_p, geom_flux_c, geom_flux_m

! Local !
REAL*8 :: dx_loc, dy_loc, dz_loc

select case(coordinate)
!-----------------------------------------------------------------------------!
case(cartesian)
	geom_flux_p = 1.0d0
	geom_flux_m = 1.0d0
  CALL COORD_DX(j_in,k_in,l_in,dx_loc,dy_loc,dz_loc)
  select case(dir_in)
	case (x_dir)
		geom_flux_c = dx_loc
	case (y_dir)
		geom_flux_c = dy_loc
	case (z_dir)
		geom_flux_c = dz_loc
	end select

!-----------------------------------------------------------------------------!
end select

END SUBROUTINE