!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
!
! mhd.f90, written by H.S. Leon Chan is 2024
! Uses Flux-CT (Balsara 1999 b, Toth 2000) to ensure divergence-less 
! 2nd order accurate scheme
! Updated to the Upwind-CT Scheme (Gardiner and Stone 2015)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************************
!
! Constrained transport on the mangetic field
!
!***********************************************************************
SUBROUTINE flux_ct
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l

! Real !
REAL*8 :: sum1, sum2 

! Real, geometric factor !
REAL*8 :: geom_ecorn_x_m,  geom_ecorn_x_c, geom_ecorn_x_p
REAL*8 :: geom_ecorn_y_m,  geom_ecorn_y_c, geom_ecorn_y_p
REAL*8 :: geom_ecorn_z_m,  geom_ecorn_z_c, geom_ecorn_z_p

!---------------------------------------------------------------------------------------------------------!

! Find cell-centered electric fields !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz + 1
  DO k = 0, ny + 1
    DO j = 0, nx + 1
      ecell (iex,j,k,l) = prim(ivz,j,k,l)*bcell(iby,j,k,l) - prim(ivy,j,k,l)*bcell(ibz,j,k,l)
      ecell (iey,j,k,l) = prim(ivx,j,k,l)*bcell(ibz,j,k,l) - prim(ivz,j,k,l)*bcell(ibx,j,k,l)
      ecell (iez,j,k,l) = prim(ivy,j,k,l)*bcell(ibx,j,k,l) - prim(ivx,j,k,l)*bcell(iby,j,k,l)
    END DO
  END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------------------!

! upwind constrained transport !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Add emf !
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l) - 0.25D0*(ecell(iex,j,k,l) + ecell(iex,j,k,l+1) + ecell(iex,j,k+1,l) + ecell(iex,j,k+1,l+1))

      ecorn(iey,j,k,l) = ecorn(iey,j,k,l) - 0.25D0*(ecell(iey,j,k,l) + ecell(iey,j,k,l+1) + ecell(iey,j+1,k,l) + ecell(iey,j+1,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l) - 0.25D0*(ecell(iez,j,k,l) + ecell(iez,j,k+1,l) + ecell(iez,j+1,k,l) + ecell(iez,j+1,k+1,l))

    END DO
  END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------------------!

! Update rungekutta operator !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) 
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! dbx/dt !
      l_rk(ibx,j,k,l) = - ((ecorn(iez,j,k,l) - ecorn(iez,j,k-1,l))/(dy) - (ecorn(iey,j,k,l) - ecorn(iey,j,k,l-1))/(dz))

      ! dby/dt !
      l_rk(iby,j,k,l) = - ((ecorn(iex,j,k,l) - ecorn(iex,j,k,l-1))/(dz) - (ecorn(iez,j,k,l) - ecorn(iez,j-1,k,l))/(dx))

      ! dbz/dt !
      l_rk(ibz,j,k,l) = - ((ecorn(iey,j,k,l) - ecorn(iey,j-1,k,l))/(dx) - (ecorn(iex,j,k,l) - ecorn(iex,j,k-1,l))/(dy))

    END DO
  END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------------------!

END SUBROUTINE

