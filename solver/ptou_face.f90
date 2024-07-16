!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given primLitive variables, consLtruct consLervative variables 
! This is the interface function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_to_U_face(j_in, k_in, l_in, epsL_in, epsR_in)
!$acc routine (PtoU_MHD_face) seq
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input !
REAL*8, INTENT(IN) :: epsL_in, epsR_in
INTEGER, INTENT(IN) :: j_in, k_in, l_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**************************************************************!
! Choose according to the fluid model !
CALL PtoU_MHD_face(j_in, k_in, l_in, epsL_in, epsR_in)

!********************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given primLitive variables, consLtruct consLervative variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PtoU_MHD_face(j_in, k_in, l_in, epsL_in, epsR_in)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input integer !
REAL*8, INTENT(IN) :: epsL_in, epsR_in
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Local array !
REAL*8 :: vsquare, bsquare

!--------------------------------------------------------------------------!
! Do the conversion, left face !

! Get v*v and b*b !
vsquare = dot_product(primL(ivx:ivz,j_in,k_in,l_in), primL(ivx:ivz,j_in,k_in,l_in))
bsquare = dot_product(primL(ibx:ibz,j_in,k_in,l_in), primL(ibx:ibz,j_in,k_in,l_in))	  

! Assign consLervative variables to core hyrodynamic variables !
consL(irho,j_in,k_in,l_in) = primL(irho,j_in,k_in,l_in)
consL(ivx:ivz,j_in,k_in,l_in) = primL(ivx:ivz,j_in,k_in,l_in)*primL(irho,j_in,k_in,l_in)
consL(itau,j_in,k_in,l_in) = primL(irho,j_in,k_in,l_in)*(epsL_in + 0.5D0*vsquare) + 0.5D0*bsquare

! Magnetic field at face center !
consL(ibx:ibz,j_in,k_in,l_in) = primL(ibx:ibz,j_in,k_in,l_in)

! For any scalar variables 
IF(itau + 1 .ne. ibx) THEN
  consL(itau+1:ibx-1,j_in,k_in,l_in) = primL(itau+1:ibx-1,j_in,k_in,l_in)*primL(irho,j_in,k_in,l_in)
END IF

!--------------------------------------------------------------------------!
! Do the conversion, right face !

! Get v*v and b*b !
vsquare = dot_product(primR(ivx:ivz,j_in,k_in,l_in), primR(ivx:ivz,j_in,k_in,l_in))
bsquare = dot_product(primR(ibx:ibz,j_in,k_in,l_in), primR(ibx:ibz,j_in,k_in,l_in))	  

! Assign consRervative variables to core hyrodynamic variables !
consR(irho,j_in,k_in,l_in) = primR(irho,j_in,k_in,l_in)
consR(ivx:ivz,j_in,k_in,l_in) = primR(ivx:ivz,j_in,k_in,l_in)*primR(irho,j_in,k_in,l_in)
consR(itau,j_in,k_in,l_in) = primR(irho,j_in,k_in,l_in)*(epsR_in + 0.5D0*vsquare) + 0.5D0*bsquare

! Magnetic field at face center !
consR(ibx:ibz,j_in,k_in,l_in) = primR(ibx:ibz,j_in,k_in,l_in)

! For any scalar variables 
IF(itau + 1 .ne. ibx) THEN
  consR(itau+1:ibx-1,j_in,k_in,l_in) = primR(itau+1:ibx-1,j_in,k_in,l_in)*primR(irho,j_in,k_in,l_in)
END IF

!--------------------------------------------------------------------------!

END SUBROUTINE