!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Builiding initial model 
! Newtonian torus with the Paczynsky-Wiita pseudo-Newtonean gravitational potential
! Solve the integral equations H + phi + h**2/(2-2q)/s**(2-2q) = c_const
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

!---------------------------------------------------------

! Integer !
INTEGER :: i, j, k, l

!---------------------------------------------------------

! Coordinate !
REAL*8 :: x_loc, y_loc, z_loc
REAL*8 :: dx_loc, dy_loc, dz_loc

! Real variables !
REAL*8 :: rho_a
REAL*8 :: lkep2
REAL*8 :: omega
REAL*8 :: k_poly
REAL*8 :: c_const
REAL*8 :: enthalpy
REAL*8 :: rand_num
REAL*8 :: rho_local

! Normalization !
REAL*8 :: pgas
REAL*8 :: pmag
REAL*8 :: pgas_max
REAL*8 :: pmag_max
REAL*8 :: beta_loc
REAL*8 :: norm
REAL*8 :: maxdb
REAL*8 :: div_b

!---------------------------------------------------------

! Vector potential !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: A_phi, A_corner

!---------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Allocate ! 
ALLOCATE(A_phi(-2:nx+3,-2:ny+3,-2:nz+3))
ALLOCATE(A_corner(-2:nx+3,-2:ny+3,-2:nz+3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Square of kepler angular momentum !
lkep2 = s_max**(2.0d0*q_grad - 1.0d0)/(s_max - r_sh)**2

! Find integration constant !
c_const = schwarzschild(s_in,0.0d0,r_sh) - lkep2*s_in**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad)

! Find enthalpy !
enthalpy = - schwarzschild(s_max,0.0d0,r_sh) + lkep2*s_max**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad) + c_const

! Find polytropic constant, defined at the position of maximum density along equator, and we set the max density to be 1 !
k_poly = enthalpy*(ggas - 1.0d0)/ggas/rho_max**(ggas - 1.0d0)

! Exit condition !
IF(k_poly < 0.0d0) THEN
  STOP 'The polytropic constant is negative'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! solve for the density profile !
rho_a = rho_fac*rho_max
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx

      ! Get coordinate !
      CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)

      ! Find enthalpy !
      enthalpy = - schwarzschild(x_loc,z_loc,r_sh) + lkep2*x_loc**(2.0d0 - 2.0d0*q_grad)/(2.0d0 - 2.0d0*q_grad) + c_const

      ! Find density !
      IF(enthalpy < 0.0d0) THEN
        enthalpy = 0.0d0
      END IF

      ! torus density and pressure
      rho_local = (enthalpy*(ggas - 1.0d0)/ggas/k_poly)**(1.0d0/(ggas - 1.0d0))

      ! Select based on conditions !
      IF (rho_local > rho_a .AND. x_loc >= s_in) THEN
        prim(irho,j,k,l) = rho_local
        omega = DSQRT(lkep2)*x_loc**(-q_grad)
        prim(ivy,j,k,l) = omega*x_loc
        prim(itau,j,k,l) = k_poly*prim(irho,j,k,l)**(ggas)
        eps(j,k,l) = k_poly*prim(irho,j,k,l)**(ggas - 1.0D0)/(ggas - 1.0D0)
      ELSE
        prim(irho,j,k,l) = rho_a
        prim(ivy,j,k,l) = 0.0d0   
        prim(itau,j,k,l) = k_poly*prim(irho,j,k,l)**(ggas)
        eps(j,k,l) = k_poly*prim(irho,j,k,l)**(ggas - 1.0D0)/(ggas - 1.0D0)        
      END IF

      ! vector potential !
      IF(prim(irho,j,k,l) >= rho_cut) THEN
        A_phi(j,k,l) = (prim(irho,j,k,l) - rho_cut)/rho_max
      END IF

    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get cell-corner vector potential !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      A_corner(j,k,l) = 0.25D0*(A_phi(j,k,l) + A_phi(j+1,k,l) + A_phi(j,k,l+1) + A_phi(j+1,k,l+1))
    END DO
  END DO
END DO

! First, initialize magnetic fields !
prim(ibx:ibz,:,:,:) = 0.0d0
bcell(ibx:ibz,:,:,:) = 0.0d0

! Get face-magnetic fields by cross product !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)
      CALL COORD_DX(j,k,l,dx_loc,dy_loc,dz_loc)
      prim(ibx,j,k,l) = - (A_corner(j,k,l) - A_corner(j,k,l-1))/(dz_loc)
      prim(ibz,j,k,l) = (xF(j)*A_corner(j,k,l) - xF(j-1)*A_corner(j-1,k,l))/(x_loc*dx_loc)
    END DO
  END DO
END DO

! Get cell-centered magnetic fields by averaging !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = 0.5D0*(prim(ibx,j,k,l) + prim(ibx,j-1,k,l))
      bcell(iby,j,k,l) = 0.5D0*(prim(iby,j,k,l) + prim(iby,j,k-1,l))
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find maxmium pressure and b^2 !
pgas_max = -1e30
pmag_max = -1e30
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      pgas = prim(itau,j,k,l)
      pmag = 0.5*(bcell(ibx,j,k,l)**2 + bcell(iby,j,k,l)**2 + bcell(ibz,j,k,l)**2)
      pgas_max = MAX(pgas_max, pgas)
      pmag_max = MAX(pmag_max, pmag)
    END DO
  END DO
END DO

! Communicate across all processor 
#ifdef MPI
CALL MPI_Allreduce(pgas_max, pgas_max, 1, MPI_DOUBLE, MPI_MAX, new_comm, ierror)
CALL MPI_Allreduce(pmag_max, pmag_max, 1, MPI_DOUBLE, MPI_MAX, new_comm, ierror)
#endif

! Assign beta !
beta_loc = pgas_max/pmag_max
norm = sqrt(beta_loc/beta_norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scale the magnetic field !

! Get face-magnetic fields by cross product !
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx
      prim(ibx,j,k,l) = prim(ibx,j,k,l)*norm
      prim(ibz,j,k,l) = prim(ibz,j,k,l)*norm
    END DO
  END DO
END DO
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      bcell(ibx,j,k,l) = bcell(ibx,j,k,l)*norm
      bcell(iby,j,k,l) = bcell(iby,j,k,l)*norm
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check divergence-B = 0 constraint !
maxdb = 0.0d0       
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j,k,l,x_loc,y_loc,z_loc)
      CALL COORD_DX(j,k,l,dx_loc,dy_loc,dz_loc)
      div_b = (xF(j)*prim(ibx,j,k,l) - xF(j-1)*prim(ibx,j-1,k,l))*dy_loc*dz_loc &
            + (prim(iby,j,k,l) - prim(iby,j,k-1,l))*dx_loc*dz_loc &
            + (prim(ibz,j,k,l) - prim(ibz,j,k,l-1))*(x_loc*dx_loc)*dy_loc
      maxdb = MAX(maxdb, div_b)
    END DO
  END DO
END DO

! Communicate across all processor 
#ifdef MPI
CALL MPI_Allreduce(maxdb, maxdb, 1, MPI_DOUBLE, MPI_MAX, new_comm, ierror)
#endif

! Print out !
WRITE (*,*)
WRITE (*,*) 'Maximum initial divergence B', maxdb
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign floor !
rho_floor = rho_fac*rho_max
eps_floor = k_poly*rho_floor**(ggas - 1.0)/(ggas - 1.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! perturb the torus !
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      IF(prim(irho,j,k,l) > rho_floor) THEN
        CALL RANDOM_NUMBER(rand_num)
        prim(irho,j,k,l) = prim(irho,j,k,l) + prim(irho,j,k,l)*(rand_num - 0.5d0)/(0.5d0)*1.0d-4
      END IF
    END DO
  END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate !
DEALLOCATE(A_phi)
DEALLOCATE(A_corner)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! choose your pseudo-newtonian potentials !

contains

  REAL*8 function schwarzschild(r_in,z_in,r_g)
  implicit none
  REAL*8 :: r_in, z_in, r_g
  schwarzschild = -1.0d0/(DSQRT(r_in**2 + z_in**2) - r_g)
  end function

END SUBROUTINE