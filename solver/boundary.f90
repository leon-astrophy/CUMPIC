!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign boundary conditions in a lump-sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY 
USE DEFINITION
IMPLICIT NONE

! Call boundary condition !
call BOUNDARYP
call BOUNDARY1D (eps, even, even, even, even, even, even)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine creates the suitable boundary for WENO (DM and NM seperately)
! Written by Leung Shing Chi in 2016
! The subroutine takes ARRAY as input/output and SIGN
! for doing odd/even parity extension
! Notice that this subroutines worked for a reduced
! size array, (1:length_step_r_part, 1:length_step_z_part)
! For full array extension, check Boundary1D_FULL.f90
! For hybrid boundaries, such as the quadrant star 
! Specific modifications are needed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARY1D (array, signx_in, signx_out, signy_in, signy_out, signz_in, signz_out)
USE DEFINITION
IMPLICIT NONE

! Input/Output array
REAL*8, INTENT (INOUT), DIMENSION (1-NGHOST:NX+3,1-NGHOST:NY+3,1-NGHOST:NZ+3) :: array

! Input parity
INTEGER, INTENT (IN) :: signx_in, signx_out
INTEGER, INTENT (IN) :: signy_in, signy_out
INTEGER, INTENT (IN) :: signz_in, signz_out

! Dummy variables
INTEGER :: i, j, k, l

! Parity factor
INTEGER :: fac_xin, fac_yin, fac_zin
INTEGER :: fac_xout, fac_yout, fac_zout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up the parity factor according to the input sign
IF(signx_in == 0) THEN
  fac_xin = 1
ELSEIF(signx_in == 1) THEN
  fac_xin = -1
END IF
IF(signx_out == 0) THEN
  fac_xout = 1
ELSEIF(signx_in == 1) THEN
  fac_xout = -1
END IF

! y-direction !
IF(signy_in == 0) THEN
  fac_yin = 1
ELSEIF(signy_in == 1) THEN
  fac_yin = -1
END IF
IF(signy_out == 0) THEN
  fac_yout = 1
ELSEIF(signy_in == 1) THEN
  fac_yout = -1
END IF

! z-direciton !
IF(signz_in == 0) THEN
  fac_zin = 1
ELSEIF(signz_in == 1) THEN
  fac_zin = -1
END IF
IF(signz_out == 0) THEN
  fac_zout = 1
ELSEIF(signz_in == 1) THEN
  fac_zout = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary

!------------------------------------------------------------------------------------
! Do the inner boundary

!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(1) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(1-j,k,l) = array(NX+1-j,k,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(1-j,k,l) = array(1,k,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(1) >= 2) THEN     
 
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)       
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(1-j,k,l) = fac_xin * array(j,k,l)
      END DO
    END DO
  ENDDO

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(2) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(NX+j,k,l) = array(j,k,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(2) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(NX+j,k,l) = array(NX,k,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(2) >= 2) THEN 

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NZ
    DO k =  1, NY
      DO j = 1, NGHOST 
        array(NX+j,k,l) = fac_xout * array(NX+1-j,k,l)
      END DO
    END DO
  ENDDO

ENDIF
!****************************************************************************
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary

!------------------------------------------------------------------------------------
! Do the inner boundary

!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(3) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,1-k,l) = array(j,NY+1-k,l)                     
      END DO
    END DO
  ENDDO
ELSEIF(boundary_flag(3) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,1-k,l) = array(j,1,l)
      END DO
    END DO
  ENDDO
ELSEIF(boundary_flag(3) >= 2) THEN   

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)                  
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,1-k,l) = fac_yin * array(j,k,l)
      END DO
    END DO
  ENDDO

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(3) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,NY+k,l) = array(j,k,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(4) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,NY+k,l) = array(j,NY,l)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(4) >= 2) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,NY+k,l) = fac_yout * array(j,NY+1-k,l)
      END DO
    END DO
  ENDDO

ENDIF
!****************************************************************************
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary

!------------------------------------------------------------------------------------
! Do the inner boundary

!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(5) == 0) THEN   

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)      
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,1-l) = array(j,k,NZ+1-l)                     
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(5) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)   
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,1-l) = array(j,k,1)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(5) >= 2) THEN   
 
  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)         
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,1-l) = fac_zin * array(j,k,l)
      END DO
    END DO
  ENDDO

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(6) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)    
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,NZ+l) = array(j,k,l)
      END DO
    END DO
  ENDDO 

ELSEIF(boundary_flag(6) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,NZ+l) = array(j,k,NZ)
      END DO
    END DO
  ENDDO

ELSEIF(boundary_flag(6) >= 2) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 1, NGHOST
    DO k =  1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        array(j,k,NZ+l) = fac_zout * array(j,k,NZ+1-l)
      END DO
    END DO
  ENDDO

END IF
!****************************************************************************
!$ACC END PARALLEL

!****************************************************************************
! MPI passing boundary values !
#ifdef MPI
CALL MPI_BOUNDARY_X(array)
CALL MPI_BOUNDARY_Y(array)
CALL MPI_BOUNDARY_Z(array)
#endif
!****************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE boundary1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To copy values to boundary ghost cell, for normal matter primitive variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNDARYP
USE DEFINITION
IMPLICIT NONE

! Dummy variables
INTEGER :: i, j, k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! x-boundary 

!------------------------------------------------------------------------------------

! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(1) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)     
  DO l = 0, NZ
    DO k = 0, NY 
      DO j = 1, NGHOST
        prim(imin:ibx-1,1-j,k,l) = prim(imin:ibx-1,NX+1-j,k,l)    
        bcell(ibx:ibz,1-j,k,l) = bcell(ibx:ibz,NX+1-j,k,l)
        prim(iby,1-j,k,l) = prim(iby,NX+1-j,k,l)
        prim(ibz,1-j,k,l) = prim(ibz,NX+1-j,k,l)
      END DO
    END DO               
  ENDDO

ELSEIF(boundary_flag(1) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 0, NY 
      DO j = 1, NGHOST
        prim(imin:ibx-1,1-j,k,l) = prim(imin:ibx-1,1,k,l)
        bcell(ibx:ibz,1-j,k,l) = bcell(ibx:ibz,1,k,l)
        prim(iby,1-j,k,l) = prim(iby,1,k,l)
        prim(ibz,1-j,k,l) = prim(ibz,1,k,l)
      END DO
    END DO               
  ENDDO

ELSEIF(boundary_flag(1) >= 2) THEN    

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 0, NY  
      DO j = 1, NGHOST
        prim(imin:ibx-1,1-j,k,l) = bfac_xin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,1-j,k,l) = bfac_xin(ibx:ibz) * bcell(ibx:ibz,j,k,l)
        prim(iby,1-j,k,l) = bfac_xin(iby)*prim(iby,j,k,l)
        prim(ibz,1-j,k,l) = bfac_xin(ibz)*prim(ibz,j,k,l)
      END DO
    END DO               
  ENDDO

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(2) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 0, NY 
      DO j = 1, NGHOST
        prim(imin:ibx-1,NX+j,k,l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,NX+j,k,l) = bcell(ibx:ibz,j,k,l)
        prim(iby,NX+j,k,l) = prim(iby,j,k,l)
        prim(ibz,NX+j,k,l) = prim(ibz,j,k,l)
      END DO
    END DO               
  ENDDO
    
ELSEIF(boundary_flag(2) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 0, NY 
      DO j = 1, NGHOST
        prim(imin:ibx-1,NX+j,k,l) = prim(imin:ibx-1,NX,k,l)
        bcell(ibx:ibz,NX+j,k,l) = bcell(ibx:ibz,NX,k,l)
        prim(iby,NX+j,k,l) = prim(iby,NX,k,l)
        prim(ibz,NX+j,k,l) = prim(ibz,NX,k,l)
      END DO
    END DO               
  ENDDO

ELSEIF(boundary_flag(2) >= 2) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 0, NY 
      DO j = 1, NGHOST
        prim(imin:ibx-1,NX+j,k,l) = bfac_xout(imin:ibx-1) * prim(imin:ibx-1,NX+1-j,k,l)
        bcell(ibx:ibz,NX+j,k,l) = bfac_xout(ibx:ibz) * bcell(ibx:ibz,NX+1-j,k,l)
        prim(iby,NX+j,k,l) = bfac_xout(iby)*prim(iby,NX+1-j,k,l)
        prim(ibz,NX+j,k,l) = bfac_xout(ibz)*prim(ibz,NX+1-j,k,l)
      END DO
    END DO               
  ENDDO

ENDIF
!****************************************************************************
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! y-boundary 

!------------------------------------------------------------------------------------
! Do the inner boundary

!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(3) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,1-k,l) = prim(imin:ibx-1,j,NY+1-k,l) 
        bcell(ibx:ibz,j,1-k,l) = bcell(ibx:ibz,j,NY+1-k,l)    
        prim(ibx,j,1-k,l) = prim(ibx,j,NY+1-k,l) 
        prim(ibz,j,1-k,l) = prim(ibz,j,NY+1-k,l) 
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(3) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,1-k,l) = prim(imin:ibx-1,j,1,l)
        bcell(ibx:ibz,j,1-k,l) = bcell(ibx:ibz,j,1,l)      
        prim(ibx,j,1-k,l) = prim(ibx,j,1,l)
        prim(ibz,j,1-k,l) = prim(ibz,j,1,l)  
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(3) >= 2) THEN    

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,1-k,l) = bfac_yin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,1-k,l) = bfac_yin(ibx:ibz) * bcell(ibx:ibz,j,k,l)  
        prim(ibx,j,1-k,l) = bfac_yin(ibx)*prim(ibx,j,k,l)  
        prim(ibz,j,1-k,l) = bfac_yin(ibz)*prim(ibz,j,k,l)   
      END DO
    END DO               
  ENDDO 

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(4) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,NY+k,l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,NY+k,l) = bcell(ibx:ibz,j,k,l) 
        prim(ibx,j,NY+k,l) = prim(ibx,j,k,l) 
        prim(ibz,j,NY+k,l) = prim(ibz,j,k,l)       
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(4) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,NY+k,l) = prim(imin:ibx-1,j,NY,l)
        bcell(ibx:ibz,j,NY+k,l) = bcell(ibx:ibz,j,NY,l)    
        prim(ibx,j,NY+k,l) = prim(ibx,j,NY,l)    
        prim(ibz,j,NY+k,l) = prim(ibz,j,NY,l)    
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(4) >= 2) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 0, NZ
    DO k = 1, NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,NY+k,l) = bfac_yout(imin:ibx-1) * prim(imin:ibx-1,j,NY+1-k,l)
        bcell(ibx:ibz,j,NY+k,l) = bfac_yout(ibx:ibz) * bcell(ibx:ibz,j,NY+1-k,l)
        prim(ibx,j,NY+k,l) = bfac_yout(ibx)*prim(ibx,j,NY+1-k,l)
        prim(ibz,j,NY+k,l) = bfac_yout(ibz)*prim(ibz,j,NY+1-k,l) 
      END DO
    END DO               
  ENDDO 

ENDIF
!****************************************************************************
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z-boundary 

!------------------------------------------------------------------------------------

! Do the inner boundary
!$ACC PARALLEL DEFAULT(PRESENT)
!****************************************************************************
IF(boundary_flag(5) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,1-l) = prim(imin:ibx-1,j,k,NZ+1-l)           
        bcell(ibx:ibz,j,k,1-l) = bcell(ibx:ibz,j,k,NZ+1-l) 
        prim(ibx,j,k,1-l) = prim(ibx,j,k,NZ+1-l)
        prim(iby,j,k,1-l) = prim(iby,j,k,NZ+1-l)                           
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(5) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,1-l) = prim(imin:ibx-1,j,k,1)
        bcell(ibx:ibz,j,k,1-l) = bcell(ibx:ibz,j,k,1)       
        prim(ibx,j,k,1-l) = prim(ibx,j,k,1)
        prim(iby,j,k,1-l) = prim(iby,j,k,1)        
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(5) >= 2) THEN  

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,1-l) = bfac_zin(imin:ibx-1) * prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,k,1-l) = bfac_zin(ibx:ibz) * bcell(ibx:ibz,j,k,l)    
        prim(ibx,j,k,1-l) = bfac_zin(ibx)*prim(ibx,j,k,l)
        prim(iby,j,k,1-l) = bfac_zin(iby)*prim(iby,j,k,l)     
      END DO
    END DO               
  ENDDO 

ENDIF

!****************************************************************************
! Do the outer boundary
IF(boundary_flag(6) == 0) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,NZ+l) = prim(imin:ibx-1,j,k,l)
        bcell(ibx:ibz,j,k,NZ+l) = bcell(ibx:ibz,j,k,l)       
        prim(ibx,j,k,NZ+l) = prim(ibx,j,k,l)
        prim(iby,j,k,NZ+l) = prim(iby,j,k,l)     
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(6) == 1) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,NZ+l) = prim(imin:ibx-1,j,k,NZ)   
        bcell(ibx:ibz,j,k,NZ+l) = bcell(ibx:ibz,j,k,NZ)       
        prim(ibx,j,k,NZ+l) = prim(ibx,j,k,NZ)
        prim(iby,j,k,NZ+l) = prim(iby,j,k,NZ)   
      END DO
    END DO               
  ENDDO 

ELSEIF(boundary_flag(6) >= 2) THEN

  !$ACC LOOP GANG WORKER VECTOR COLLAPSE(3)
  DO l = 1, NGHOST
    DO k = 1 - NGHOST, NY + NGHOST
      DO j = 1 - NGHOST, NX + NGHOST
        prim(imin:ibx-1,j,k,NZ+l) = bfac_zout(imin:ibx-1) * prim(imin:ibx-1,j,k,NZ+1-l)   
        bcell(ibx:ibz,j,k,NZ+l) = bfac_zout(ibx:ibz) * bcell(ibx:ibz,j,k,NZ+1-l)    
        prim(ibx,j,k,NZ+l) = bfac_zout(ibx)*prim(ibx,j,k,NZ+1-l)
        prim(iby,j,k,NZ+l) = bfac_zout(iby)*prim(iby,j,k,NZ+1-l)     
      END DO
    END DO               
  ENDDO 

ENDIF
!****************************************************************************
!$ACC END PARALLEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Or, add your custom boundary conditions !
CALL CUSTOM_BOUNDARY

!****************************************************************************
! MPI passing boundary values !
#ifdef MPI
CALL MPI_BOUNDARYP_X
CALL MPI_BOUNDARYP_Y
CALL MPI_BOUNDARYP_Z
#endif
!****************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE