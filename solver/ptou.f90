!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given primitive variables, construct conservative variables 
! This is the interface function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_to_U(j_in, k_in, l_in)
!$acc routine (PtoU_MHD) seq
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**************************************************************!
! Choose according to the fluid model !
CALL PtoU_MHD(j_in, k_in, l_in)

!********************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given primitive variables, construct conservative variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PtoU_MHD(j_in, k_in, l_in)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Local array !
REAL*8 :: vsquare, bsquare

!--------------------------------------------------------------------------!
! Do the conversion 

! Get v*v and b*b !
vsquare = dot_product(prim(ivx:ivz,j_in,k_in,l_in), prim(ivx:ivz,j_in,k_in,l_in))
bsquare = dot_product(bcell(ibx:ibz,j_in,k_in,l_in), bcell(ibx:ibz,j_in,k_in,l_in))	  

! Assign conservative variables to core hyrodynamic variables !
cons(irho,j_in,k_in,l_in) = prim(irho,j_in,k_in,l_in)
cons(ivx:ivz,j_in,k_in,l_in) = prim(ivx:ivz,j_in,k_in,l_in)*prim(irho,j_in,k_in,l_in)
cons(itau,j_in,k_in,l_in) = prim(irho,j_in,k_in,l_in)*(eps(j_in,k_in,l_in) + 0.5D0*vsquare) + 0.5D0*bsquare

! Magnetic field at face center !
cons(ibx:ibz,j_in,k_in,l_in) = prim(ibx:ibz,j_in,k_in,l_in)

! For any scalar variables 
IF(itau + 1 .ne. ibx) THEN
  cons(itau+1:ibx-1,j_in,k_in,l_in) = prim(itau+1:ibx-1,j_in,k_in,l_in)*prim(irho,j_in,k_in,l_in)
END IF

!--------------------------------------------------------------------------!

END SUBROUTINE