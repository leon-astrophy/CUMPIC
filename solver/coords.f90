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
!$ACC ROUTINE (GET_COORD) SEQ
!$ACC ROUTINE (COORD_DX) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Local !
REAL*8 :: bsquare
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

!-----------------------------------------------------------------------------!
! Pre compute, maybe can move it into each case to reduce computation !
CALL GET_COORD(j_in,k_in,l_in,x_loc,y_loc,z_loc)
CALL COORD_DX(j_in,k_in,l_in,dx_loc,dy_loc,dz_loc)

select case(coordinate)
!-----------------------------------------------------------------------------!
case(cylindrical)
	bsquare = dot_product(bcell(ibx:ibz,j_in,k_in,l_in),bcell(ibx:ibz,j_in,k_in,l_in))

	sc(ivx,j_in,k_in,l_in) = sc(ivx,j_in,k_in,l_in) 
												 + (prim(itau,j_in,k_in,l_in) + prim(irho,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in) & 
						   					 + 0.5D0*bsquare - bcell(iby,j_in,k_in,l_in)*bcell(iby,j_in,k_in,l_in))/(x_loc)										

	sc(ivy,j_in,k_in,l_in) = sc(ivy,j_in,k_in,l_in) 
												 - (prim(irho,j_in,k_in,l_in)*prim(ivx,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in) & 
												 - bcell(ibx,j_in,k_in,l_in)*bcell(iby,j_in,k_in,l_in))/(x_loc)

!-----------------------------------------------------------------------------!
case(spherical)
	bsquare = dot_product(bcell(ibx:ibz,j_in,k_in,l_in),bcell(ibx:ibz,j_in,k_in,l_in))

	sc(ivx,j_in,k_in,l_in) = sc(ivx,j_in,k_in,l_in) 
												 + (2.0D0*prim(itau,j_in,k_in,l_in) + bcell(ibx,j_in,k_in,l_in)*bcell(ibx,j_in,k_in,l_in) &
												 + prim(irho,j_in,k_in,l_in)*(prim(ivy,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in) &
												 + prim(ivz,j_in,k_in,l_in)*prim(ivz,j_in,k_in,l_in)))/(x_loc)

	sc(ivy,j_in,k_in,l_in) = sc(ivy,j_in,k_in,l_in) 
												 + (prim(itau,j_in,k_in,l_in) + prim(irho,j_in,k_in,l_in)*prim(ivz,j_in,k_in,l_in) &
												 * prim(ivz,j_in,k_in,l_in) + 0.5D0*bsquare - bcell(ibz,j_in,k_in,l_in)*bcell(ibz,j_in,k_in,l_in)) &
											 	 / (x_loc*DTAN(y_loc)) & 
																			
												 - (prim(irho,j_in,k_in,l_in)*prim(ivx,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in) & 
											 	 - bcell(ibx,j_in,k_in,l_in)*bcell(iby,j_in,k_in,l_in))/(x_loc) 

	sc(ivz,j_in,k_in,l_in) = sc(ivz,j_in,k_in,l_in) 
												 - (prim(irho,j_in,k_in,l_in)*prim(ivx,j_in,k_in,l_in)*prim(ivz,j_in,k_in,l_in) & 
												 - bcell(ibx,j_in,k_in,l_in)*bcell(ibz,j_in,k_in,l_in))/(x_loc)  & 

												 - (prim(irho,j_in,k_in,l_in)*prim(ivy,j_in,k_in,l_in)*prim(ivz,j_in,k_in,l_in) & 
											   - bcell(iby,j_in,k_in,l_in)*bcell(ibz,j_in,k_in,l_in))/(x_loc*DTAN(y_loc))

!-----------------------------------------------------------------------------!
end select

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the geometric factors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GEOM_FLUX(dir_in,j_in,k_in,l_in,geom_flux_p,geom_flux_c,geom_flux_m)
!$ACC ROUTINE (GET_COORD) SEQ
!$ACC ROUTINE (COORD_DX) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: dir_in, j_in, k_in, l_in

! Output !
REAL*8, INTENT (OUT) :: geom_flux_p, geom_flux_c, geom_flux_m

! Local !
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

!-----------------------------------------------------------------------------!
! Pre compute, maybe can move it into each case to reduce computation !
CALL GET_COORD(j_in,k_in,l_in,x_loc,y_loc,z_loc)
CALL COORD_DX(j_in,k_in,l_in,dx_loc,dy_loc,dz_loc)

select case(coordinate)
!-----------------------------------------------------------------------------!
case(cartesian)
	geom_flux_p = 1.0d0
	geom_flux_m = 1.0d0
  select case(dir_in)
	case (x_dir)
		geom_flux_c = dx_loc
	case (y_dir)
		geom_flux_c = dy_loc
	case (z_dir)
		geom_flux_c = dz_loc
	end select
!-----------------------------------------------------------------------------!
case(cylindrical)
  select case(dir_in)
	case (x_dir)
		geom_flux_p = xF(j_in)
		geom_flux_m = xF(j_in-1)
		geom_flux_c = 0.5d0*(xF(j_in)**2 - xF(j_in-1)**2) 
	case (y_dir)
		geom_flux_p = 1.0D0
		geom_flux_m = 1.0D0
		geom_flux_c = x_loc*dy_loc
	case (z_dir)
		geom_flux_p = 1.0D0
		geom_flux_m = 1.0D0
		geom_flux_c = dz_loc
	end select
!-----------------------------------------------------------------------------!
case(spherical)
  select case(dir_in)
	case (x_dir)
		geom_flux_p = xF(j_in)**2
		geom_flux_m = xF(j_in-1)**2
		geom_flux_c = (xF(j_in)**3 - xF(j_in-1)**3)/3.0d0
	case (y_dir)
		geom_flux_p = DSIN(yF(k_in))
		geom_flux_m = DSIN(yF(k_in-1))
		geom_flux_c = (2.0d0*(xF(j_in)**3 - xF(j_in-1)**3)/3.0d0/(xF(j_in)**2 - xF(j_in-1)**2))*(DCOS(yF(k_in-1)) - DCOS(yF(k_in)))
	case (z_dir)
		geom_flux_p = 1.0d0
		geom_flux_m = 1.0d0
		geom_flux_c = (2.0d0*(xF(j_in)**3 - xF(j_in-1)**3)/3.0d0/(xF(j_in)**2 - xF(j_in-1)**2))*(DCOS(yF(k_in-1)) - DCOS(yF(k_in)))*dz_loc/dy_loc
	end select
!-----------------------------------------------------------------------------!
end select

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given index tuple (j,k,l), find the geometric factors in flux CT scheme
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GEOM_FLUX(j_in, k_in, l_in, g_bx_ez_m, g_bx_ez_c, g_bx_ez_p, g_bx_ey_m, g_bx_ey_c, g_bx_ey_p, &
g_by_ex_m, g_by_ex_c, g_by_ex_p, g_by_ez_m, g_by_ez_c, g_by_ez_p, g_bz_ex_m, g_bz_ex_c, g_bz_ex_p, g_bz_ey_m, &
g_bz_ey_c, g_bz_ey_p)
!$ACC ROUTINE (GET_COORD) SEQ
!$ACC ROUTINE (COORD_DX) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Output !
REAL*8, INTENT (OUT) :: g_bx_ez_m, g_bx_ez_c, g_bx_ez_p, g_bx_ey_m, g_bx_ey_c, g_bx_ey_p, &
g_by_ex_m, g_by_ex_c, g_by_ex_p, g_by_ez_m, g_by_ez_c, g_by_ez_p, g_bz_ex_m, g_bz_ex_c, g_bz_ex_p, g_bz_ey_m, &
g_bz_ey_c, g_bz_ey_p

!-----------------------------------------------------------------------------!
! Pre compute, maybe can move it into each case to reduce computation !
CALL GET_COORD(j_in,k_in,l_in,x_loc,y_loc,z_loc)
CALL COORD_DX(j_in,k_in,l_in,dx_loc,dy_loc,dz_loc)

select case(coordinate)
!-----------------------------------------------------------------------------!
case(cartesian)
	g_bx_ez_m	= 1.0d0
	g_bx_ez_c = dy_loc
	g_bx_ez_p = 1.0d0

	g_bx_ey_m	= 1.0d0
	g_bx_ey_c = dz_loc
	g_bx_ey_p = 1.0d0

	g_by_ex_m	= 1.0d0
	g_by_ex_c = dz_loc
	g_by_ex_p = 1.0d0

	g_by_ez_m	= 1.0d0
	g_by_ez_c = dx_loc
	g_by_ez_p = 1.0d0

	g_bz_ey_m	= 1.0d0
	g_bz_ey_c = dx_loc
	g_bz_ey_p = 1.0d0

	g_bz_ex_m	= 1.0d0
	g_bz_ex_c = dy_loc
	g_bz_ex_p = 1.0d0
!-----------------------------------------------------------------------------!
case(cylindrical)
	g_bx_ez_m	= 1.0d0
	g_bx_ez_c = dy_loc*xF(j_in)
	g_bx_ez_p = 1.0d0

	g_bx_ey_m	= 1.0d0
	g_bx_ey_c = dz_loc
	g_bx_ey_p = 1.0d0

	g_by_ex_m	= 1.0d0
	g_by_ex_c = dz_loc
	g_by_ex_p = 1.0d0

	g_by_ez_m	= 1.0d0
	g_by_ez_c = dx_loc
	g_by_ez_p = 1.0d0

	g_bz_ey_m	= xF(j_in-1)
	g_bz_ey_c = x_loc*dx_loc
	g_bz_ey_p = xF(j_in)

	g_bz_ex_m	= 1.0d0
	g_bz_ex_c = x_loc*dy_loc
	g_bz_ex_p = 1.0d0
!-----------------------------------------------------------------------------!
case(spherical)
	g_bx_ez_m	= DSIN(yF(k_in-1))
	g_bx_ez_c = xF(j_in)*(DCOS(yF(k_in-1)) - DCOS(yF(k_in)))
	g_bx_ez_p = DSIN(yF(k_in))

	g_bx_ey_m	= 1.0d0
	g_bx_ey_c = xF(j_in)*(DCOS(yF(k_in-1)) - DCOS(yF(k_in)))*dz_loc/dy_loc
	g_bx_ey_p = 1.0d0

	g_by_ex_m	= 1.0d0
	g_by_ex_c = x_loc*DSIN(yF(k_in))*dz_loc
	g_by_ex_p = 1.0d0

	g_by_ez_m	= xF(j_in-1)
	g_by_ez_c = x_loc*dx_loc
	g_by_ez_p = xF(j_in)

	g_bz_ey_m	= xF(j_in-1)
	g_bz_ey_c = x_loc*dx_loc
	g_bz_ey_p = xF(j_in)

	g_bz_ex_m	= 1.0d0
	g_bz_ex_c = x_loc*dy_loc
	g_bz_ex_p = 1.0d0
!-----------------------------------------------------------------------------!
end select

END SUBROUTINE