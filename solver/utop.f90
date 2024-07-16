!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given conservative variables, invert to get primitive variables 
! This is the interface function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE U_to_P(j_in, k_in, l_in)
!$acc routine (UtoP_MHD) seq 
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************************************************!
! Choose according to fluid model !
CALL UtoP_MHD(j_in, k_in, l_in)

!********************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given conservative variables, invert to get primitive variables 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UtoP_MHD(j_in, k_in, l_in)
!$acc routine seq
USE DEFINITION
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: j_in, k_in, l_in

! Local array !
REAL*8 :: vsquare, bsquare

!--------------------------------------------------------------------------!
! Do the conversion 

! Core primitive variables !
prim(irho,j_in,k_in,l_in) = cons(irho,j_in,k_in,l_in)
prim(ivx:ivz,j_in,k_in,l_in) = cons(ivx:ivz,j_in,k_in,l_in)/cons(irho,j_in,k_in,l_in)

! Magnetic fields !
prim(ibx:ibz,j_in,k_in,l_in) = cons(ibx:ibz,j_in,k_in,l_in)

! Get v*v and b*b !
bsquare = dot_product(bcell(ibx:ibz,j_in,k_in,l_in), bcell(ibx:ibz,j_in,k_in,l_in))
vsquare = dot_product(prim(ivx:ivz,j_in,k_in,l_in), prim(ivx:ivz,j_in,k_in,l_in))

! For any scalar variables 
IF(itau + 1 .ne. ibx) THEN
  prim(itau+1:ibx-1,j_in,k_in,l_in) = cons(itau+1:ibx-1,j_in,k_in,l_in)/cons(irho,j_in,k_in,l_in)
END IF

! epsilon here, CAUTION: no negativity check !
eps(j_in,k_in,l_in) = (cons(itau,j_in,k_in,l_in) - 0.5D0*bsquare)/cons(irho,j_in,k_in,l_in) - 0.5D0 * vsquare

!--------------------------------------------------------------------------!

END SUBROUTINE