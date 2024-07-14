!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Riemann1D.f90, storing subroutines for generating 
! Initial conditions for riemann problem test in 1D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------!

!r0 = 0.1D0
!r1 = 0.115D0
!DO k = 0, ny
!  DO j = 0, nx
!    CALL GET_COORD(j,k,0,x,y,dummy)
!    r = DSQRT((x - 0.5D0)**2 + (y - 0.5D0)**2)
!    fr = (r1 - r)/(r1 - r0)
!    If(r <= r0) THEN
!      prim(irho,j,k,:) = 10.0D0
!      prim(ivx,j,k,:) = -1.0D0*(y - 0.5D0)/r0
!      prim(ivy,j,k,:) = 1.0D0*(x - 0.5D0)/r0
!    ELSEIF(r >= r1) THEN
!      prim(irho,j,k,:) = 1.0D0
!      prim(ivx,j,k,:) = 0.0D0
!      prim(ivy,j,k,:) = 0.0D0
!    ELSE
!      prim(irho,j,k,:) = 1.0D0 + 9.0d0*fr
!      prim(ivx,j,k,:) = -fr*1.0D0*(y - 0.5D0)/r
!      prim(ivy,j,k,:) = fr*1.0D0*(x - 0.5D0)/r
!    END IF
!    prim(itau,j,k,:) = 0.5D0
!    prim(ibx,j,k,:) = 2.5D0/DSQRT(4.0D0*pi)
!    eps(j,k,:) = prim(itau,j,k,:) / prim(irho,j,k,:) / (ggas - 1.0D0)
!  ENDDO
!ENDDO