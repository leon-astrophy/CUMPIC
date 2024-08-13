!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Reconstruction.f90 - contains wrapper for calling per-cell
! reconstruction subroutine
! Inputs (p_-2, p_-1, p_c, p_+1, p_+2)
! Outputs (p_l, p_r)
! 
! | -2 | -1 | c | +1 | +2 |
!             L | R
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************************************
!
! Interface that decide reconstruction method
!
!*******************************************************************
SUBROUTINE INTERPOLATE (dir_in, j_in, k_in, l_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE (TVD_MM) SEQ
!$ACC ROUTINE (PPMC) SEQ  
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Index tuple !
INTEGER, INTENT (IN) :: dir_in, j_in, k_in, l_in

! Input conservative variables !
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Select based on case 
!-----------------------------------------------------------------------
select case(RECON)

! TVD MM !
!-----------------------------------------------------------------------
case(TVDMM)
  CALL TVD_MM (dir_in, j_in, k_in, l_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)

! PPM !
!-----------------------------------------------------------------------
case(PPMC)
  CALL PPM (dir_in, j_in, k_in, l_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)

!-----------------------------------------------------------------------
end select

END SUBROUTINE

!*******************************************************************
!
! TVD with minmod limiter, reference: Mignone 2014
!
!*******************************************************************
SUBROUTINE TVD_MM (dir_in, j_in, k_in, l_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE (coord_dx) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION
IMPLICIT NONE

! Index tuple !
INTEGER, INTENT (IN) :: dir_in, j_in, k_in, l_in

! Input conservative variables !
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Real !
REAL*8 :: dummy
REAL*8 :: dh_m1, dh_c, dh_p1

REAL*8 :: dm1, dc, dp1
REAL*8 :: cf, cb, dq2
REAL*8 :: dqb, dqf
REAL*8 :: slope, v_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

select case(dir_in)
! Choose by direction !
case(x_dir)

  ! Get geometric factor !
  CALL COORD_DX(j_in-1,k_in,l_in,dh_m1,dummy,dummy)
  CALL COORD_DX(j_in,k_in,l_in,dh_c,dummy,dummy)
  CALL COORD_DX(j_in+1,k_in,l_in,dh_p1,dummy,dummy)

case(y_dir)

  ! Get geometric factor !
  CALL COORD_DX(j_in,k_in-1,l_in,dummy,dh_m1,dummy)
  CALL COORD_DX(j_in,k_in,l_in,dummy,dh_c,dummy)
  CALL COORD_DX(j_in,k_in+1,l_in,dummy,dh_p1,dummy)

case(z_dir)

  ! Get geometric factor !
  CALL COORD_DX(j_in,k_in,l_in-1,dummy,dummy,dh_m1)
  CALL COORD_DX(j_in,k_in,l_in,dummy,dummy,dh_c)
  CALL COORD_DX(j_in,k_in,l_in+1,dummy,dummy,dh_p1)

end select

! Assign !
cf = (dh_c + dh_p1)/dh_c
cb = (dh_c + dh_m1)/dh_c

! Assign !
dqf = 2.0d0*(vp1 - vc)/cf
dqb = 2.0d0*(vc - vm1)/cb

! Assign slope !
slope = 0.5d0*(SIGN(1.0d0, dqf) + SIGN(1.0d0, dqb))*MIN(ABS(dqf), ABS(dqb))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Output !
vp_out = vc + 0.5d0*slope
vm_out = vc - 0.5d0*slope
  
!----------------------------------------------------------------------!

END SUBROUTINE

!*******************************************************************
!
! Interpolate cell average values to interface values using PPM 
! Assumed non-uniform gridding. Using the original PPM algorithm
! See Colella 1984. No steepening or flattening is performed.
!
!*******************************************************************
SUBROUTINE PPM (dir_in, j_in, k_in, l_in, vm2, vm1, vc, vp1, vp2, vm_out, vp_out)
!$ACC ROUTINE (coord_dx) SEQ
!$ACC ROUTINE SEQ 
USE DEFINITION
IMPLICIT NONE
  
! Index tuple !
INTEGER, INTENT (IN) :: dir_in, j_in, k_in, l_in

! Input primitive variables !
REAL*8, INTENT (IN) :: vm2, vm1, vc, vp1, vp2

! The reconstructed states at cell boundary !
REAL*8, INTENT (OUT) :: vm_out, vp_out

! Temporal variables !
REAL*8 :: vl, vr
REAL*8 :: deltam1, deltac, deltap1
REAL*8 :: dmp1, dmc, dmm1
REAL*8 :: condition, dummy

! grid size !
REAL*8 :: dm2, dm1, dc, dp1, dp2

! Reconstruction coefficients !
REAL*8 :: a0m, a1m, a2m, a3m
REAL*8 :: a0p, a1p, a2p, a3p
REAL*8 :: b0m, b1m
REAL*8 :: b0c, b1c
REAL*8 :: b0p, b1p

! Real !
REAL*8 :: a1R, aR, deltaXR, z1R, zR
REAL*8 :: a1L, aL, deltaXL, z1L, zL
REAL*8 :: cp1, cp2, cp3
REAL*8 :: cc1, cc2, cc3
REAL*8 :: cm1, cm2, cm3
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Choose by direction !
select case(dir_in)
  case(x_dir)
  
    ! Get geometric factor !
    CALL COORD_DX(j_in-2,k_in,l_in,dm2,dummy,dummy)
    CALL COORD_DX(j_in-1,k_in,l_in,dm1,dummy,dummy)
    CALL COORD_DX(j_in,k_in,l_in,dc,dummy,dummy)
    CALL COORD_DX(j_in+1,k_in,l_in,dp1,dummy,dummy)
    CALL COORD_DX(j_in+2,k_in,l_in,dp2,dummy,dummy)
  
  case(y_dir)
  
    ! Get geometric factor !
    CALL COORD_DX(j_in,k_in-2,l_in,dummy,dm2,dummy)
    CALL COORD_DX(j_in,k_in-1,l_in,dummy,dm1,dummy)
    CALL COORD_DX(j_in,k_in,l_in,dummy,dc,dummy)
    CALL COORD_DX(j_in,k_in+1,l_in,dummy,dp1,dummy)
    CALL COORD_DX(j_in,k_in+2,l_in,dummy,dp2,dummy)
  
  case(z_dir)
  
    ! Get geometric factor !
    CALL COORD_DX(j_in,k_in,l_in-2,dummy,dummy,dm2)
    CALL COORD_DX(j_in,k_in,l_in-1,dummy,dummy,dm1)
    CALL COORD_DX(j_in,k_in,l_in,dummy,dummy,dc)
    CALL COORD_DX(j_in,k_in,l_in+1,dummy,dummy,dp1)
    CALL COORD_DX(j_in,k_in,l_in+2,dummy,dummy,dp2)
  
  end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! Get coefficient for interpolants!
a1R = dc/(dc + dp1)
a1L = dm1/(dm1 + dc)
aR = 2.0D0*dp1*dc/(dp1 + dc)
aL = 2.0D0*dc*dm1/(dc + dm1)
deltaXR = dm1 + dc + dp1 + dp2
deltaXL = dm2 + dm1 + dc + dp1
z1R = (dm1 + dc)/(2.0D0*dc + dp1)
z1L = (dm2 + dm1)/(2.0D0*dm1 + dc)
zR = (dp2 + dp1)/(2.0D0*dp1 + dc)
zL = (dp1 + dc)/(2.0D0*dc + dm1)

! For slope estimations !
cc1 = dc/(dm1 + dc + dp1)
cp1 = dp1/(dc + dp1 + dp2)
cm1 = dm1/(dm2 + dm1 + dc)
cc2 = (2.0D0*dm1 + dc)/(dp1 + dc)
cp2 = (2.0D0*dc + dp1)/(dp2 + dp1)
cm2 = (2.0D0*dm2 + dm1)/(dc + dm1)
cc3 = (dc + 2.0D0*dp1)/(dm1 + dc)
cp3 = (dp1 + 2.0D0*dp2)/(dc + dp1)
cm3 = (dm1 + 2.0D0*dc)/(dm2 + dm1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign weight !
a0m = a1L
a1m = aL*(z1L - zL)/deltaXL
a2m = (-dm1*z1L)/deltaXL
a3m = dc*zL/deltaXL
a0p = a1R
a1p = aR*(z1R - zR)/deltaXR
a2p = (-dc*z1R)/deltaXR
a3p = dp1*zR/deltaXR
b0m = cm1*cm2
b1m = cm1*cm3
b0c = cc1*cc2
b1c = cc1*cc3
b0p = cp1*cp2
b1p = cp1*cp3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get slopes !
deltap1 = b0p*(vp2 - vp1) + b1p*(vp1 - vc)
deltac = b0c*(vp1 - vc) + b1c*(vc - vm1)
deltam1 = b0m*(vc - vm1) + b1m*(vm1 - vm2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Van Leer limiter !

condition = (vp1 - vc)*(vc - vm1)
IF(condition > 0.0D0) THEN
  dmc = min(abs(deltac), 2.0D0*abs(vc - vm1), 2.0D0*abs(vc - vp1))*sign(1.0D0,deltac)
ELSE
  dmc = 0.0D0
END IF

condition = (vp2 - vp1)*(vp1 - vc)
IF(condition > 0.0D0) THEN
  dmp1 = min(abs(deltap1), 2.0D0*abs(vp1 - vc), 2.0D0*abs(vp1 - vp2))*sign(1.0D0,deltap1)
ELSE
  dmp1 = 0.0D0
END IF

condition = (vc - vm1)*(vm1 - vm2)
IF(condition > 0.0D0) THEN
  dmm1 = min(abs(deltam1), 2.0D0*abs(vm1 - vm2), 2.0D0*abs(vm1 - vc))*sign(1.0D0,deltam1)
ELSE
  dmm1 = 0.0D0
END IF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get the interpolant !
vp_out = vc + a0p*(vp1 - vc) + a1p*(vp1 - vc) + a2p*dmp1 + a3p*dmc
vm_out = vm1 + a0m*(vc - vm1) + a1m*(vc - vm1) + a2m*dmc + a3m*dmm1

! backup !
vl = vm_out
vr = vp_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check for extremum conditions !
IF((vr - vc)*(vc - vl) <= 0.0D0) THEN
  vp_out = vc
  vm_out = vc
ELSE
  condition = (vr - vl)*(vl - 3.0D0*vc + 2.0D0*vr)
  IF(condition < 0.0D0) THEN
    vm_out = 3.0d0*vc - 2.0D0*vr
  END IF
  condition = (vp_out - vm_out)*(3.0D0*vc - 2.0d0*vl - vr)
  IF(condition < 0.0D0) THEN
    vp_out = 3.0d0*vc - 2.0D0*vl
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE