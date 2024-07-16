!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given a direction dir_in, a tuple (j,k,l), primitive and conservative variables
! build the flux. This is the interface function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_to_flux(dir_in, j_in, k_in, l_in)
!$ACC ROUTINE (Ptof_MHD) SEQ
!$ACC ROUTINE SEQ
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in, j_in, k_in, l_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*****************************************************************!
! Choose according to fluid model !
CALL Ptof_MHD(dir_in, j_in, k_in, l_in)

!*****************************************************************!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Given a direction dir_in, a tuple (j,k,l), primitive and conservative variables
! build the flux
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Ptof_MHD(dir_in, j_in, k_in, l_in)
!$ACC ROUTINE SEQ
USE DEFINITION 
IMPLICIT NONE

! Integer !
INTEGER, INTENT(IN) :: dir_in, j_in, k_in, l_in

! Integer !
INTEGER :: i, j, k, l
INTEGER :: ivn, ivt1, ivt2
INTEGER :: ibn, ibt1, ibt2

! Dummy !
REAL*8 :: vsquare, bsquare, vdotb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Assign !
select case(dir_in)
!**************************************!
case(x_dir)
  ivn = ivx
  ivt1 = ivy
  ivt2 = ivz
  ibn = ibx
  ibt1 = iby
  ibt2 = ibz
case(y_dir)
  ivn = ivy
  ivt1 = ivz
  ivt2 = ivx
  ibn = iby
  ibt1 = ibz
  ibt2 = ibx
case(z_dir)
  ivn = ivz
  ivt1 = ivx
  ivt2 = ivy
  ibn = ibz
  ibt1 = ibx
  ibt2 = iby
!**************************************!
end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux, left face !

! get dot product !
bsquare = dot_product(primL(ibx:ibz,j_in,k_in,l_in), primL(ibx:ibz,j_in,k_in,l_in))
vsquare = dot_product(primL(ivx:ivz,j_in,k_in,l_in), primL(ivx:ivz,j_in,k_in,l_in))
vdotb = dot_product(primL(ibx:ibz,j_in,k_in,l_in), primL(ivx:ivz,j_in,k_in,l_in))

! First multiply by normal velocity !
fluxL (imin:ibx-1,j_in,k_in,l_in) = consL (imin:ibx-1,j_in,k_in,l_in) * primL(ivn,j_in,k_in,l_in)

! Add the pressure term to normal momentum equation !
fluxL (ivn,j_in,k_in,l_in) = fluxL (ivn,j_in,k_in,l_in) + primL(itau,j_in,k_in,l_in) + 0.5D0*bsquare - primL(ibn,j_in,k_in,l_in)*primL(ibn,j_in,k_in,l_in)

! Adjust the transverse momentum equation !
fluxL(ivt1,j_in,k_in,l_in) = fluxL(ivt1,j_in,k_in,l_in) - primL(ibn,j_in,k_in,l_in)*primL(ibt1,j_in,k_in,l_in)
fluxL(ivt2,j_in,k_in,l_in) = fluxL(ivt2,j_in,k_in,l_in) - primL(ibn,j_in,k_in,l_in)*primL(ibt2,j_in,k_in,l_in)

! Add the presusre work done term to the energy equation             
fluxL (itau,j_in,k_in,l_in) = fluxL (itau,j_in,k_in,l_in) + (primL(itau,j_in,k_in,l_in) + 0.5D0*bsquare) * primL(ivn,j_in,k_in,l_in) - primL(ibn,j_in,k_in,l_in)*vdotb

! normal and transverse magnetic fields !
fluxL(ibn,j_in,k_in,l_in) = 0.0D0
fluxL(ibt1,j_in,k_in,l_in) = (primL(ivn,j_in,k_in,l_in)*primL(ibt1,j_in,k_in,l_in) - primL(ivt1,j_in,k_in,l_in)*primL(ibn,j_in,k_in,l_in))
fluxL(ibt2,j_in,k_in,l_in) = (primL(ivn,j_in,k_in,l_in)*primL(ibt2,j_in,k_in,l_in) - primL(ivt2,j_in,k_in,l_in)*primL(ibn,j_in,k_in,l_in))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get flux, right face !

! get dot product !
bsquare = dot_product(primR(ibx:ibz,j_in,k_in,l_in), primR(ibx:ibz,j_in,k_in,l_in))
vsquare = dot_product(primR(ivx:ivz,j_in,k_in,l_in), primR(ivx:ivz,j_in,k_in,l_in))
vdotb = dot_product(primR(ibx:ibz,j_in,k_in,l_in), primR(ivx:ivz,j_in,k_in,l_in))

! First multiply by normal velocity !
fluxR (imin:ibx-1,j_in,k_in,l_in) = consR (imin:ibx-1,j_in,k_in,l_in) * primR(ivn,j_in,k_in,l_in)

! Add the pressure term to normal momentum equation !
fluxR (ivn,j_in,k_in,l_in) = fluxR (ivn,j_in,k_in,l_in) + primR(itau,j_in,k_in,l_in) + 0.5D0*bsquare - primR(ibn,j_in,k_in,l_in)*primR(ibn,j_in,k_in,l_in)

! Adjust the transverse momentum equation !
fluxR(ivt1,j_in,k_in,l_in) = fluxR(ivt1,j_in,k_in,l_in) - primR(ibn,j_in,k_in,l_in)*primR(ibt1,j_in,k_in,l_in)
fluxR(ivt2,j_in,k_in,l_in) = fluxR(ivt2,j_in,k_in,l_in) - primR(ibn,j_in,k_in,l_in)*primR(ibt2,j_in,k_in,l_in)

! Add the presusre work done term to the energy equation             
fluxR (itau,j_in,k_in,l_in) = fluxR (itau,j_in,k_in,l_in) + (primR(itau,j_in,k_in,l_in) + 0.5D0*bsquare) * primR(ivn,j_in,k_in,l_in) - primR(ibn,j_in,k_in,l_in)*vdotb

! normal and transverse magnetic fields !
fluxR(ibn,j_in,k_in,l_in) = 0.0D0
fluxR(ibt1,j_in,k_in,l_in) = (primR(ivn,j_in,k_in,l_in)*primR(ibt1,j_in,k_in,l_in) - primR(ivt1,j_in,k_in,l_in)*primR(ibn,j_in,k_in,l_in))
fluxR(ibt2,j_in,k_in,l_in) = (primR(ivn,j_in,k_in,l_in)*primR(ibt2,j_in,k_in,l_in) - primR(ivt2,j_in,k_in,l_in)*primR(ibn,j_in,k_in,l_in))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE
  