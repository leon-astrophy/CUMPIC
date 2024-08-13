!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine does one single Runge-Kutta full step
! It uses the opeator splitting and separate
! all non-gravitational source term to be done 
! after the hydro step.
! Written by Leung Shing Chi in 2016
! Updated by Leung Shing Chi in 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RUNGEKUTTA
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!---------------------------------------------------------------------------------------------!     

! Backup old arrays !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, NZ
  DO k = 0, NY
    DO j = 0, NX
			DO i = imin, imax 
				u_old (i,j,k,l) = cons (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------!
! 1st iteration

! Discretize !
CALL SPATIAL

! NM sector !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, NZ
  DO k = 0, NY
    DO j = 0, NX
			DO i = imin, imax 
				cons (i,j,k,l) = u_old (i,j,k,l) + dt * l_rk (i,j,k,l) 
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Convert from conservative to primitive
CALL FROMUTORVE

! Check density !
CALL CUSTOM_CHECKRHO

! Do conversion again !
CALL FROMRVETOU

! set boundary conditions !
CALL BOUNDARY

! Update 
CALL UPDATE (1)

!---------------------------------------------------------------------------------------------!
! 2nd iteration

! Discretize !
CALL SPATIAL

! NM sector !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, NZ
  DO k = 0, NY
    DO j = 0, NX
			DO i = imin, imax 
				cons (i,j,k,l) = rk20 * u_old(i,j,k,l) + rk21 * cons (i,j,k,l) + rk22 * dt * l_rk (i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Convert from conservative to primitive
CALL FROMUTORVE

! Check density !
CALL CUSTOM_CHECKRHO

! Do conversion again !
CALL FROMRVETOU

! set boundary conditions !
call BOUNDARY

! Update 
CALL UPDATE (2)

!---------------------------------------------------------------------------------------------!
! Prepare for next step

! Discretize !
CALL SPATIAL

! NM sector !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 0, NZ
  DO k = 0, NY
    DO j = 0, NX
			DO i = imin, imax 
				cons (i,j,k,l) = rk30 * u_old(i,j,k,l) + rk31 * cons (i,j,k,l) + rk32 * dt * l_rk (i,j,k,l)
			END DO
		END DO
	END DO
END DO 
!$ACC END PARALLEL

! Convert from conservative to primitive
CALL FROMUTORVE 

! Check density !
CALL CUSTOM_CHECKRHO

!*****************************************************************
! Section for operator splitting

CALL OPERATOR_SPLIT

!****************************************************************

! Check density again !
CALL CUSTOM_CHECKRHO

! Update again !
CALL FROMRVETOU

! set boundary conditions !
call BOUNDARY

! Update physical quantities
CALL UPDATE (3)

!---------------------------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! This subroutine calculates the maximum time step
! which satisfies the Courant condition 
! Written by Leung Shing Chi in 2016   
! If you modify the Euler equation, make sure you change this 
! part to include the new effective sound speed
! Limiters are posed based on output time and running time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE finddt
!$ACC ROUTINE (EOS_SOUNDSPEED) SEQ
!$ACC ROUTINE (GET_COORD) SEQ
!$ACC ROUTINE (coord_dx) SEQ
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
include "mpif.h"
#endif

! Dummy variables
INTEGER :: i, j, k, l

! Soundspeed !
REAL*8 :: cs_loc

! For MHD speed !
REAL*8 :: a2_mhd, b2_mhd
REAL*8 :: a4_mhd, b4_mhd
REAL*8 :: b2x_mhd, b2y_mhd, b2z_mhd
REAL*8 :: cfx_mhd, cfy_mhd, cfz_mhd

! Local maximum effective speed
REAL*8 :: lambda, lambda1, lambda2, lambda3

! Local grid size !
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

! Local minimum dt for DM, NM and 1st overlayer
REAL*8 :: dt_out, dt_temp

! Geometric factor !
REAL*8 :: geom_length_x
REAL*8 :: geom_length_y
REAL*8 :: geom_length_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initilaize !
dt_temp = 1.0d10

! Now we find the minimum time constrained by NM sector
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) &
!$ACC PRIVATE(cs_loc, a2_mhd, b2_mhd, a4_mhd, b4_mhd, b2x_mhd, b2y_mhd, b2z_mhd, & 
!$ACC cfx_mhd, cfy_mhd, cfz_mhd, lambda, lambda1, lambda2, lambda3, &
!$ACC dx_loc, dy_loc, dz_loc, geom_length_x, geom_length_y, geom_length_z) REDUCTION(MIN:dt_temp)
DO l = 1, NZ
	DO k = 1, NY
		DO j = 1, NX

			! find speed of sound cs !
			CALL EOS_SOUNDSPEED (prim(itau,j,k,l), prim(irho,j,k,l), eps(j,k,l), cs_loc)

			! Only grid with density above threshold density is counted
			a2_mhd = cs_loc*cs_loc
			a4_mhd = a2_mhd*a2_mhd
			b2x_mhd = (bcell(ibx,j,k,l)*bcell(ibx,j,k,l)/prim(irho,j,k,l))
			b2y_mhd = (bcell(iby,j,k,l)*bcell(iby,j,k,l)/prim(irho,j,k,l))
			b2z_mhd = (bcell(ibz,j,k,l)*bcell(ibz,j,k,l)/prim(irho,j,k,l))
			b2_mhd = b2x_mhd + b2y_mhd + b2z_mhd
			b4_mhd = b2_mhd*b2_mhd
			cfx_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2x_mhd)))
			cfy_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2y_mhd)))
			cfz_mhd = DSQRT(0.5D0*(a2_mhd + b2_mhd + DSQRT((a4_mhd + 2.0d0*a2_mhd*b2_mhd + b4_mhd) - 4.0D0*a2_mhd*b2z_mhd)))
			lambda1 = ABS(prim(ivx,j,k,l)) + cfx_mhd
			lambda2 = ABS(prim(ivy,j,k,l)) + cfy_mhd
			lambda3 = ABS(prim(ivz,j,k,l)) + cfz_mhd
			lambda = MAX(lambda1, lambda2, lambda3)
                        
			! Find loca grid size !
			CALL COORD_DX (j,k,l,dx_loc,dy_loc,dz_loc)
			CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)

			! Geometric factor !
			select case(coordinate)
			case(cylindrical)
				dy_loc = x_loc*dy_loc
			case(spherical)
				dy_loc = x_loc*dy_loc
				dz_loc = x_loc*DSIN(y_loc)*dz_loc
			end select
			
			! Compute time step !
			dt_out = MIN(dx_loc,dy_loc,dz_loc)*cfl/lambda
			dt_temp = MIN(dt_out, dt_temp)
	                	
		END DO
	ENDDO
ENDDO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Only the minimum one is chosen
#ifdef MPI
CALL MPI_Allreduce(dt_temp, dt, 1, MPI_DOUBLE, MPI_MIN, new_comm, ierror)
#else
dt = dt_temp
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE finddt
