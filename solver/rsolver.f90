!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! rsolver.f90, contains all Riemann solvers 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************
!
! Interface that decide reconstruction method
!
!*******************************************************************
SUBROUTINE RIEMANN (dir_in, j_in, k_in, l_in)
!$ACC ROUTINE (HLL_SOLVER) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Index tuple !
INTEGER, INTENT (IN) :: dir_in

! Location index !
INTEGER, INTENT (IN) :: j_in, k_in, l_in

! Select based on case 
!-----------------------------------------------------------------------
select case(SOLVER)

! TVD MM !
!-----------------------------------------------------------------------
case(HLL)
  CALL HLL_SOLVER (dir_in, j_in, k_in, l_in)

!-----------------------------------------------------------------------
end select

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The HLL flux, see for example, TORO's book on numerical hydrodynamics 
!			    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLL_SOLVER (dir_in, j_in, k_in, l_in)
!$ACC ROUTINE (EOS_EPSILON) SEQ
!$ACC ROUTINE (P_to_U_face) SEQ
!$ACC ROUTINE (P_to_flux) SEQ
!$ACC ROUTINE (EOS_SOUNDSPEED) SEQ 
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE
  
! integer !
INTEGER, INTENT(IN) :: dir_in

! Location index !
INTEGER, INTENT (IN) :: j_in, k_in, l_in

! Signal speeds !
REAL*8 :: sL, sR
REAL*8 :: epsL, epsR
REAL*8 :: csL, csR
REAL*8 :: cfsL, cfsR
REAL*8 :: ubar, cbar

! Integer !
INTEGER :: i
INTEGER :: ivn
INTEGER :: ibn, ibt1, ibt2

!***********************************************************************************!

! Assign !
select case(dir_in)
!**************************************!
case(x_dir)
  ibn = ibx 
  ibt1 = iby
  ibt2 = ibz
  ivn = ivx
case(y_dir)
  ibn = iby
  ibt1 = ibz
  ibt2 = ibx
  ivn = ivy
case(z_dir)
  ibn = ibz
  ibt1 = ibx
  ibt2 = iby
  ivn = ivz
!**************************************!
end select
  
!***********************************************************************************!

! Find epsilon !
CALL EOS_EPSILON (primL(itau, j_in, k_in, l_in), primL(irho, j_in, k_in, l_in), epsL)
CALL EOS_EPSILON (primR(itau, j_in, k_in, l_in), primR(irho, j_in, k_in, l_in), epsR)

! Find the sound speed !
CALL EOS_SOUNDSPEED (primL(itau, j_in, k_in, l_in), primL(irho, j_in, k_in, l_in), epsL, csL)
CALL EOS_SOUNDSPEED (primR(itau, j_in, k_in, l_in), primR(irho, j_in, k_in, l_in), epsR, csR)

! First, compute conservative variables !
CALL P_to_U_face(j_in, k_in, l_in, epsL, epsR)

! Then, compute local fluxes !
CALL P_to_flux(dir_in, j_in, k_in, l_in)

! Signal speed !
cfsL = compute_signalspeed(csL, primL(ibn, j_in, k_in, l_in), primL(ibt1, j_in, k_in, l_in), primL(ibt2, j_in, k_in, l_in), primL(irho, j_in, k_in, l_in))
cfsR = compute_signalspeed(csR, primR(ibn, j_in, k_in, l_in), primR(ibt1, j_in, k_in, l_in), primR(ibt2, j_in, k_in, l_in), primR(irho, j_in, k_in, l_in))
ubar = compute_roe(primL(ivn, j_in, k_in, l_in),primR(ivn, j_in, k_in, l_in),primL(irho, j_in, k_in, l_in),primR(irho, j_in, k_in, l_in))
cbar = compute_roe(cfsL,cfsR,primL(irho, j_in, k_in, l_in),primR(irho, j_in, k_in, l_in))
sL = min(primL(ivn, j_in, k_in, l_in) - cfsL, ubar - cbar)
sR = max(primR(ivn, j_in, k_in, l_in) + cfsR, ubar + cbar)

! Find the flux !
IF(sL >= 0.0D0) THEN
  flux(imin:imax, j_in, k_in, l_in) = fluxL(imin:imax, j_in, k_in, l_in)
ELSEIF(sL <= 0.0D0 .AND. sR >= 0.0D0) THEN
  DO i = imin, imax
    flux(i, j_in, k_in, l_in) = compute_fluxhll(fluxL(i,j_in,k_in,l_in), fluxR(i,j_in,k_in,l_in), consL(i,j_in,k_in,l_in), consR(i,j_in,k_in,l_in), sL, sR) 
  END DO
ELSEIF(sR <= 0.0D0) THEN
  flux(imin:imax, j_in, k_in, l_in) = fluxR(imin:imax, j_in, k_in, l_in)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE