!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine prepares the data for spatial discretization,
! due ask the WENO_module to do the reconstruction
! and then combines the results for one Runge-Kutta sub-step
! Prototype developed by Wong Ka Wing in 2010 (or before?)
! Extended by Leung Shing Chi in 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SPATIAL
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i,j,k,l

!---------------------------------------------------------------------------------------------!
! First, initialize !

! Initialize source term and rungekutta operator !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) 
DO l = 0, nz
	DO k = 0, ny
		DO j = 0, nx
			DO i = imin, imax
				sc(i,j,k,l) = 0.0D0
				l_rk(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Initialize electric field !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) 
DO l = 0, nz
	DO k = 0, ny
		DO j = 0, nx
			DO i = iex, iez
				ecorn(i,j,k,l) = 0.0D0
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------!
! Find source term !

! Predefined source term !
CALL GET_SOURCE

! Custom source term !
CALL CUSTOM_SOURCE

!---------------------------------------------------------------------------------------------!
! Now compute the interface flux !

! Interface flux, sweep through all direction !
CALL GET_FLUXES

! Constrained transport !
CALL flux_ct

!---------------------------------------------------------------------------------------------!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get geometric source terms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_SOURCE
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Geometric sources terms
!-----------------------------------------------------------------------------------!
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) 
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, imax
				
				! Call warpper for computing geometric sources !
				CALL GEOM_SOURCE(j, k, l, prim(:,j,k,l), bcell(:,j,k,l), sc(:,j,k,l))

			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 1. Peform reconstruction to get primitive variables at cell interfaces 
! 2. Construct conservative variables and fluxes at cell interfaces 
! 3. Solve the Riemann problem
! 4. Solve for the flux difference
! 5. Update the rungekutta operator 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_FLUXES
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, l

! Geometric factor !
REAL*8 :: geom_flux_p, geom_flux_c, geom_flux_m

!==============================================================================================================!

! First loop through the x-direction
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (x_dir, j, k, l, prim(i,j-2:j+2,k,l), primR(i,j-1,k,l), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (x_dir, j, k, l, bcell(iby,j-2:j+2,k,l), primR(iby,j-1,k,l), primL(iby,j,k,l))
			CALL INTERPOLATE (x_dir, j, k, l, bcell(ibz,j-2:j+2,k,l), primR(ibz,j-1,k,l), primL(ibz,j,k,l))
				
			! Special treatment for normal field 
			primR(ibx,j-1,k,l) = prim(ibx,j-1,k,l)
			primL(ibx,j,k,l) = prim(ibx,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Then solve the Riemann Problem !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx

			! Core hydrodynamic variables !
			CALL RIEMANN (x_dir, primL(imin:imax,j,k,l), primR(imin:imax,j,k,l), flux(imin:imax,j,k,l))

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Add the flux difference into the l-operator
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIAVTE(geom_flux_p, geom_flux_c, geom_flux_m)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(x_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j-1,k,l)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Now do the same for the electric field, using flux-CT scheme 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iey,j,k,l) = ecorn(iey,j,k,l) + 0.50D0*(flux(ibz,j,k,l) + flux(ibz,j,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l) + 0.50D0*(- flux(iby,j,k,l) - flux(iby,j,k+1,l))

    END DO
  END DO
END DO
!$ACC END PARALLEL
!==============================================================================================================!

! Then loop through the y-direction
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (y_dir, j, k, l, prim(i,j,k-2:k+2,l), primR(i,j,k-1,l), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (y_dir, j, k, l, bcell(ibx,j,k-2:k+2,l), primR(ibx,j,k-1,l), primL(ibx,j,k,l))
			CALL INTERPOLATE (y_dir, j, k, l, bcell(ibz,j,k-2:k+2,l), primR(ibz,j,k-1,l), primL(ibz,j,k,l))

			! Special treatment for normal field 
			primR(iby,j,k-1,l) = prim(iby,j,k-1,l)
			primL(iby,j,k,l) = prim(iby,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Then solve the Riemann Problem !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz + 1 
	DO k = 0, ny 
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			CALL RIEMANN (y_dir, primL(imin:imax,j,k,l), primR(imin:imax,j,k,l), flux(imin:imax,j,k,l))

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Add the flux difference into the l-operator
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIAVTE(geom_flux_p, geom_flux_c, geom_flux_m)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(y_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j,k-1,l)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Now do the same for the electric field, using flux-CT scheme 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l) + 0.50D0*(- flux(ibz,j,k,l) - flux(ibz,j,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l) + 0.50D0*(flux(ibx,j,k,l) + flux(ibx,j+1,k,l))
			
    END DO
  END DO
END DO
!$ACC END PARALLEL
!==============================================================================================================!

! Finally loop through the z-direction
IF(NZ > 1) THEN ! This switch might or might not be problematic
!--------------------------------------------------------------------------------------------------------------!
! Interpolate to get L/R state !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz + 1 
	DO k = 0, ny + 1
		DO j = 0, nx + 1
			
			! Core hydrodynamic variables !
			DO i = imin, ibx - 1
				CALL INTERPOLATE (z_dir, j, k, l, prim(i,j,k,l-2:l+2), primR(i,j,k,l-1), primL(i,j,k,l))
			END DO

			! Extra scalar !

			! Cell center magnetic fields !
			CALL INTERPOLATE (z_dir, j, k, l, bcell(ibx,j,k,l-2:l+2), primR(ibx,j,k,l-1), primL(ibx,j,k,l))
			CALL INTERPOLATE (z_dir, j, k, l, bcell(iby,j,k,l-2:l+2), primR(iby,j,k,l-1), primL(iby,j,k,l))

			! Special treatment for normal field 
			primR(ibz,j,k,l-1) = prim(ibz,j,k,l-1)
			primL(ibz,j,k,l) = prim(ibz,j,k,l)

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Then solve the Riemann Problem !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
DO l = 0, nz 
	DO k = 0, ny + 1
		DO j = 0, nx + 1

			! Core hydrodynamic variables !
			CALL RIEMANN (z_dir, primL(imin:imax,j,k,l), primR(imin:imax,j,k,l), flux(imin:imax,j,k,l))

		END DO
	END DO
END DO
!$ACC END PARALLEL

! Add the flux difference into the l-operator
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT) PRIAVTE(geom_flux_p, geom_flux_c, geom_flux_m)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1

				! Call the geometric factor !
				CALL GEOM_FLUX(z_dir,j,k,l,geom_flux_p,geom_flux_c,geom_flux_m)

				! Perform subtraction !
				l_rk(i,j,k,l) = l_rk(i,j,k,l) - (flux (i,j,k,l)*geom_flux_p - flux (i,j,k,l-1)*geom_flux_m) / geom_flux_c

			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL

! Now do the same for the electric field, using flux-CT scheme 
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l) + 0.50D0*(flux(iby,j,k,l) + flux(iby,j,k+1,l))

      ecorn(iey,j,k,l) = ecorn(iey,j,k,l) + 0.50D0*(- flux(ibx,j,k,l) - flux(ibx,j+1,k,l))

    END DO
  END DO
END DO
!$ACC END PARALLEL
!--------------------------------------------------------------------------------------------------------------!
ENDIF

! At the end, add source terms !
!--------------------------------------------------------------------------------------------------------------!
! Final step, get rungekutta operator, LHS of the hydro equation !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(4) DEFAULT(PRESENT)
DO l = 1, nz
	DO k = 1, ny
		DO j = 1, nx
			DO i = imin, ibx - 1
				l_rk(i,j,k,l) = l_rk(i,j,k,l) + sc(i,j,k,l)
			END DO
		END DO
	END DO
END DO
!$ACC END PARALLEL
!==============================================================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given a direction dir_in, a tuple (j,k,l), primitive and conservative variables
! build the flux
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_to_flux(dir_in,prim_in,bcell_in,cons_in,flux_out)
!$ACC ROUTINE SEQ
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in

! Input/Output array
REAL*8, INTENT (IN), DIMENSION (ibx:ibz) :: bcell_in
REAL*8, INTENT (IN), DIMENSION (1:no_of_eq) :: prim_in
REAL*8, INTENT (IN), DIMENSION (1:no_of_eq) :: cons_in

! Output !
REAL*8, INTENT (OUT), DIMENSION (1:no_of_eq) :: flux_out

! Integer !
INTEGER :: i, j, k, l
INTEGER :: ivn, ivt1, ivt2
INTEGER :: ibn, ibt1, ibt2

! Dummy !
REAL*8 :: vsquare, bsquare, vdotb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
select case(dir_in)
case(x_dir)
	ivn = ivx
	ivt1 = ivy
	ivt2 = ivz
	ibn = ibx
	ibt1 = iby
	ibt2 = ibz
case(y_dir)
	ivn = ivy
	ivt1 = ivz
	ivt2 = ivx
	ibn = iby
	ibt1 = ibz
	ibt2 = ibx
case(z_dir)
	ivn = ivz
	ivt1 = ivx
	ivt2 = ivy
	ibn = ibz
	ibt1 = ibx
	ibt2 = iby
end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux !

! get dot product !
bsquare = dot_product(prim_in(ibx:ibz), prim_in(ibx:ibz))
vsquare = dot_product(prim_in(ivx:ivz), prim_in(ivx:ivz))
vdotb = dot_product(prim_in(ibx:ibz), prim_in(ivx:ivz))

! First multiply by normal velocity !
flux_out (imin:ibx-1) = cons_in (imin:ibx-1) * prim_in(ivn)

! Add the pressure term to normal momentum equation !
flux_out (ivn) = flux_out (ivn) + prim_in(itau) + 0.5D0*(bsquare) - prim_in(ibn)*prim_in(ibn)

! Adjust the transverse momentum equation !
flux_out(ivt1) = flux_out(ivt1) - prim_in(ibn)*prim_in(ibt1)
flux_out(ivt2) = flux_out(ivt2) - prim_in(ibn)*prim_in(ibt2)

! Add the presusre work done term to the energy equation             
flux_out (itau) = flux_out (itau) + (prim_in(itau) + 0.5D0*(bsquare)) * prim_in(ivn) - prim_in(ibn)*vdotb

! normal and transverse magnetic fields !
flux_out(ibn) = 0.0D0
flux_out(ibt1) = (prim_in(ivn)*prim_in(ibt1) - prim_in(ivt1)*prim_in(ibn))
flux_out(ibt2) = (prim_in(ivn)*prim_in(ibt2) - prim_in(ivt2)*prim_in(ibn))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
