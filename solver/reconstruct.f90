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