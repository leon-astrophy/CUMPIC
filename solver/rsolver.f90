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
SUBROUTINE RIEMANN (dir_in, pl_in, pr_in, f_out)
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Index tuple !
INTEGER, INTENT (IN) :: dir_in

! Input variables !
REAL*8, INTENT (IN), DIMENSION(1:no_of_eq) :: pl_in
REAL*8, INTENT (IN), DIMENSION(1:no_of_eq) :: pr_in

! Output fluxes !
REAL*8, INTENT (OUT), DIMENSION(1:no_of_eq) :: f_out

! Select based on case 
!-----------------------------------------------------------------------
select case(SOLVER)

! TVD MM !
!-----------------------------------------------------------------------
case(HLL)
  CALL HLL_SOLVER (dir_in, pl_in, pr_in, f_out)

!-----------------------------------------------------------------------
end select

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The HLL flux, see for example, TORO's book on numerical hydrodynamics 
!			    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HLL_SOLVER (dir_in, pl_in, pr_in, f_out)
!$ACC ROUTINE SEQ
USE DEFINITION 
IMPLICIT NONE
  
! integer !
INTEGER, INTENT(IN) :: dir_in

! Input variables !
REAL*8, INTENT (IN), DIMENSION(1:no_of_eq) :: pl_in
REAL*8, INTENT (IN), DIMENSION(1:no_of_eq) :: pr_in

! Output fluxes !
REAL*8, INTENT (OUT), DIMENSION(1:no_of_eq) :: f_out

! local conservative variables  
REAL*8, DIMENSION(1:no_of_eq) :: ul
REAL*8, DIMENSION(1:no_of_eq) :: ur

! local fluxes
REAL*8, DIMENSION(1:no_of_eq) :: fl
REAL*8, DIMENSION(1:no_of_eq) :: fr

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
end select
  
!***********************************************************************************!

! Find epsilon !
CALL EOS_EPSILON (pl_in(irho), pl_in(itau), epsL)
CALL EOS_EPSILON (pr_in(irho), pr_in(itau), epsR)

! First, compute conservative variables !
CALL P_to_U(pl_in, pl_in(ibx:ibz), epsL, ul)
CALL P_to_U(pr_in, pr_in(ibx:ibz), epsR, ur)

! Then, compute local fluxes !
CALL P_to_flux(dir_in,pl_in,pl_in(ibx:ibz),ul,fl)
CALL P_to_flux(dir_in,pr_in,pr_in(ibx:ibz),ur,fr)

! Find the sound speed !
CALL EOS_SOUNDSPEED (pl_in(itau), pl_in(irho), csL)
CALL EOS_SOUNDSPEED (pr_in(itau), pr_in(irho), csR)

! Signal speed !
cfsL = compute_signalspeed(csL, pl_in(ibn), pl_in(ibt1), pl_in(ibt2), pl_in(irho))
cfsR = compute_signalspeed(csR, pr_in(ibn), pr_in(ibt1), pr_in(ibt2), pr_in(irho))
ubar = compute_roe(pl_in(ivn),pr_in(ivn),pl_in(irho),pr_in(irho))
cbar = compute_roe(cfsL,cfsR,pl_in(irho),pr_in(irho))
sL = min(pl_in(ivn) - cfsL, ubar - cbar)
sR = max(pr_in(ivn) + cfsR, ubar + cbar)

! Find the flux !
IF(sL >= 0.0D0) THEN
  f_out(imin:imax) = fl(imin:imax)
ELSEIF(sL <= 0.0D0 .AND. sR >= 0.0D0) THEN
  DO i = imin, imax
    f_out(i) = compute_fluxhll(fl(i),fr(i),uL(i),uR(i),sL,sR)
  END DO
ELSEIF(sR <= 0.0D0) THEN
  f_out(imin:imax) = fr(imin:imax)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE