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
! Find the divergence of magnetic field 
!
!***********************************************************************
SUBROUTINE find_divb
!$ACC ROUTINE (GEOM_AREA) SEQ
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
include "mpif.h"
#endif

! Integer !
INTEGER :: j, k, l

! Real !
REAL*8 :: maxdb, div_b

! For the divergence condition !
REAL*8 :: axp
REAL*8 :: axm
REAL*8 :: ayp
REAL*8 :: aym
REAL*8 :: azp
REAL*8 :: azm

!---------------------------------------------------------------------------------------------------------!

! Check divergence-B = 0 constraint !
maxdb = 0.0d0       
DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GEOM_AREA(j,k,l,axp,axm,ayp,aym,azp,azm)
      div_b = (axp*prim(ibx,j,k,l) - axm*prim(ibx,j-1,k,l)) &
            + (ayp*prim(iby,j,k,l) - aym*prim(iby,j,k-1,l)) &
            + (azp*prim(ibz,j,k,l) - azm*prim(ibz,j,k,l-1))
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
WRITE (*,*) 'Maximum divergence B', maxdb
WRITE (*,*)

!---------------------------------------------------------------------------------------------------------!

END SUBROUTINE
  
!***********************************************************************
!
! Constrained transport on the mangetic field
!
!***********************************************************************
SUBROUTINE flux_ct
!$ACC ROUTINE (GEOM_CT) SEQ
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: j, k, l

! Real !
REAL*8 :: sum1, sum2 

! Real, geometric factor !
REAL*8 :: g_bx_ez_m, g_bx_ez_c, g_bx_ez_p
REAL*8 :: g_bx_ey_m, g_bx_ey_c, g_bx_ey_p
REAL*8 :: g_by_ex_m, g_by_ex_c, g_by_ex_p
REAL*8 :: g_by_ez_m, g_by_ez_c, g_by_ez_p
REAL*8 :: g_bz_ex_m, g_bz_ex_c, g_bz_ex_p
REAL*8 :: g_bz_ey_m, g_bz_ey_c, g_bz_ey_p

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
      ecorn(iex,j,k,l) = ecorn(iex,j,k,l)*2.0d0 - 0.25D0*(ecell(iex,j,k,l) + ecell(iex,j,k,l+1) + ecell(iex,j,k+1,l) + ecell(iex,j,k+1,l+1))

      ecorn(iey,j,k,l) = ecorn(iey,j,k,l)*2.0d0 - 0.25D0*(ecell(iey,j,k,l) + ecell(iey,j,k,l+1) + ecell(iey,j+1,k,l) + ecell(iey,j+1,k,l+1))

      ecorn(iez,j,k,l) = ecorn(iez,j,k,l)*2.0d0 - 0.25D0*(ecell(iez,j,k,l) + ecell(iez,j,k+1,l) + ecell(iez,j+1,k,l) + ecell(iez,j+1,k+1,l))

    END DO
  END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------------------!

! Update rungekutta operator !
!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE(g_bx_ez_m, g_bx_ez_c, g_bx_ez_p, &
!$ACC g_bx_ey_m, g_bx_ey_c, g_bx_ey_p, g_by_ex_m, g_by_ex_c, g_by_ex_p, g_by_ez_m, g_by_ez_c, g_by_ez_p, &
!$ACC g_bz_ex_m, g_bz_ex_c, g_bz_ex_p, g_bz_ey_m, g_bz_ey_c, g_bz_ey_p)
DO l = 0, nz
  DO k = 0, ny
    DO j = 0, nx

      ! Get geometric factor !
      CALL GEOM_CT(j,k,l,g_bx_ez_m, g_bx_ez_c, g_bx_ez_p, g_bx_ey_m, g_bx_ey_c, g_bx_ey_p, g_by_ex_m, g_by_ex_c, g_by_ex_p, &
      g_by_ez_m, g_by_ez_c, g_by_ez_p, g_bz_ex_m, g_bz_ex_c, g_bz_ex_p, g_bz_ey_m, g_bz_ey_c, g_bz_ey_p)

      ! dbx/dt !
      l_rk(ibx,j,k,l) = - ((ecorn(iez,j,k,l)*g_bx_ez_p - ecorn(iez,j,k-1,l)*g_bx_ez_m)/g_bx_ez_c - (ecorn(iey,j,k,l)*g_bx_ey_p - ecorn(iey,j,k,l-1)*g_bx_ey_m)/g_bx_ey_c)

      ! dby/dt !
      l_rk(iby,j,k,l) = - ((ecorn(iex,j,k,l)*g_by_ex_p - ecorn(iex,j,k,l-1)*g_by_ex_m)/g_by_ex_c - (ecorn(iez,j,k,l)*g_by_ez_p - ecorn(iez,j-1,k,l)*g_by_ez_m)/g_by_ez_c)

      ! dbz/dt !
      l_rk(ibz,j,k,l) = - ((ecorn(iey,j,k,l)*g_bz_ey_p - ecorn(iey,j-1,k,l)*g_bz_ey_m)/g_bz_ey_c - (ecorn(iex,j,k,l)*g_bz_ex_p - ecorn(iex,j,k-1,l)*g_bz_ex_m)/g_bz_ex_c)

    END DO
  END DO
END DO
!$ACC END PARALLEL

!---------------------------------------------------------------------------------------------------------!

END SUBROUTINE

