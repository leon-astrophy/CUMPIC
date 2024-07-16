!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Written by Leung Shing Chi in 2016
! The subroutine setup the primitive and conservative equation indexes
! Notice that if you want your to add your own quantities, you need to specify them here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SETUP_EQNS
USE DEFINITION
IMPLICIT NONE

#ifdef MPI
include "mpif.h"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get rank of MPI !
#ifdef MPI
	call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
#else
	mpi_rank = 0
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize no_of_eq
no_of_eq = 0

! Initialize imax imin !
imin = 0
imax = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for buliding equation indexes !

IF(mpi_rank == 0) THEN
        WRITE(*,*) 
	WRITE(*,*) 'Now we arrange no_of_eq accordingly'
	WRITE(*,*) 'For each variables, no_of_eq increases by 1'
	WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we do the NM sector

! Set up minimum equation for NM !
imin = no_of_eq + 1

! NM density
no_of_eq = no_of_eq + 1
imax = no_of_eq
irho = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make irho = ', no_of_eq
END IF

! NM vel-x
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivx = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make ivx = ', no_of_eq
END IF

! NM vel-y
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivy = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make ivy = ', no_of_eq
END IF

! NM vel-z
no_of_eq = no_of_eq + 1
imax = no_of_eq
ivz = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make ivz = ', no_of_eq
END IF

! NM epsilon
no_of_eq = no_of_eq + 1
imax = no_of_eq
itau = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make itau = ', no_of_eq
END IF

! Custom equations !
CALL CUSTOM_EQN

! Magnetic fields, Bx !
no_of_eq = no_of_eq + 1
imax = no_of_eq
ibx = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make ibx = ', no_of_eq
END IF

! Magnetic fields, By !
no_of_eq = no_of_eq + 1
imax = no_of_eq
iby = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make iby = ', no_of_eq
END IF

! Magnetic fields, Bz !
no_of_eq = no_of_eq + 1
imax = no_of_eq
ibz = no_of_eq
IF(mpi_rank == 0) THEN 
	write(*,*) 'Make ibz = ', no_of_eq
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs according to boundary conditions !

! First allocate !
ALLOCATE(bfac_xin(1:no_of_eq))
ALLOCATE(bfac_yin(1:no_of_eq))
ALLOCATE(bfac_zin(1:no_of_eq))
ALLOCATE(bfac_xout(1:no_of_eq)) 
ALLOCATE(bfac_yout(1:no_of_eq)) 
ALLOCATE(bfac_zout(1:no_of_eq)) 

! Assume all variables are even !
bfac_xin(:) = 1
bfac_yin(:) = 1
bfac_zin(:) = 1
bfac_xout(:) = 1
bfac_yout(:) = 1
bfac_zout(:) = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, x inner boundary !
IF(boundary_flag(1) == 2) THEN
	bfac_xin(ivx) = -1
	bfac_xin(ibx) = -1
ELSEIF(boundary_flag(1) == 3) THEN
	bfac_xin(ivx) = -1
	bfac_xin(ibx) = -1
	bfac_xin(ivz) = -1
	bfac_xin(ibz) = -1
ELSEIF(boundary_flag(1) == 4) THEN
	IF(mpi_rank == 0) THEN 
		WRITE (*,*) 'Equatorial symmetry is not allowed for the x-boundary'
	END IF
	STOP
END IF

! Flip signs, x outer boundary !
IF(boundary_flag(2) == 2) THEN
	bfac_xout(ivx) = -1
	bfac_xout(ibx) = -1
ELSEIF(boundary_flag(2) == 3) THEN
	bfac_xout(ivx) = -1
	bfac_xout(ibx) = -1
	bfac_xout(ivz) = -1
	bfac_xout(ibz) = -1
ELSEIF(boundary_flag(2) == 4) THEN
	IF(mpi_rank == 0) THEN 
		WRITE (*,*) 'Equatorial symmetry is not allowed for the x-boundary'
	END IF
	STOP
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, y inner boundary !
IF(boundary_flag(3) == 2) THEN
	bfac_yin(ivy) = -1
	bfac_yin(iby) = -1
ELSEIF(boundary_flag(3) == 3) THEN
	bfac_yin(ivy) = -1
	bfac_yin(iby) = -1
	bfac_yin(ivz) = -1
	bfac_yin(ibz) = -1
ELSEIF(boundary_flag(3) == 4) THEN
	bfac_yin(ivy) = -1
	bfac_yin(ibx) = -1
	bfac_yin(ibz) = -1
END IF

! Flip signs, y outer boundary !
IF(boundary_flag(4) == 2) THEN
	bfac_yout(ivy) = -1
	bfac_yout(iby) = -1
ELSEIF(boundary_flag(4) == 3) THEN
	bfac_yout(ivy) = -1
	bfac_yout(iby) = -1
	bfac_yout(ivz) = -1
	bfac_yout(ibz) = -1
ELSEIF(boundary_flag(4) == 4) THEN
	bfac_yout(ivy) = -1
	bfac_yout(ibx) = -1
	bfac_yout(ibz) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Flip signs, z inner boundary !
IF(boundary_flag(5) == 2) THEN
	bfac_zin(ivz) = -1
	bfac_zin(ibz) = -1
ELSEIF(boundary_flag(5) == 3) THEN
	IF(mpi_rank == 0) THEN 
		WRITE (*,*) 'Axial symmetry is not allowed for the z-boundary'
	END IF
	STOP
ELSEIF(boundary_flag(5) == 4) THEN
	bfac_zin(ivz) = -1
	bfac_zin(ibx) = -1
	bfac_zin(iby) = -1
END IF

! Flip signs, z outer boundary !
IF(boundary_flag(6) == 2) THEN
	bfac_zout(ivz) = -1
	bfac_zout(ibz) = -1
ELSEIF(boundary_flag(6) == 3) THEN
	IF(mpi_rank == 0) THEN 
		WRITE (*,*) 'Axial symmetry is not allowed for the z-boundary'
	END IF
	STOP
ELSEIF(boundary_flag(6) == 4) THEN
	bfac_zout(ivz) = -1
	bfac_zout(ibx) = -1
	bfac_zout(iby) = -1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for printing out !

IF(mpi_rank == 0) THEN 
	WRITE(*,*)
	WRITE(*,*) 'Finished setting up equation indexes'
	WRITE(*,*) 'There are', no_of_eq, 'number of equations'
	WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE SETUP_EQNS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILD_HYDRO
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Core arrays for solving hyperbolic PDE

! Primitive variables 
ALLOCATE(prim(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! Conservative variables 
ALLOCATE(cons(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! Epsilon !
ALLOCATE(eps(1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! L/R primitive variables at cell interfaces !
ALLOCATE(primL(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))
ALLOCATE(primR(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! L/R conservative variables at cell interfaces !
ALLOCATE(consL(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))
ALLOCATE(consR(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! L/R flux variables at cell interfaces !
ALLOCATE(fluxL(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))
ALLOCATE(fluxR(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! Rungekutta time evolution !
ALLOCATE(l_rk(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))
ALLOCATE(u_old(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! Flux and source terms !
ALLOCATE(sc(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))
ALLOCATE(flux(1:no_of_eq,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! cell centered magnetic fields !
ALLOCATE(bcell(ibx:ibz,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! cell centered electric fields !
ALLOCATE(ecell(iex:iez,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

! cell cornered electrin fields
ALLOCATE(ecorn(iex:iez,1-NGHOST:NX+NGHOST,1-NGHOST:NY+NGHOST,1-NGHOST:NZ+NGHOST))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! face coordinate of x !
ALLOCATE(xf(-NGHOST:NX+NGHOST))

! face coordinate of y !
ALLOCATE(yf(-NGHOST:NY+NGHOST))

! face coordinate of z !
ALLOCATE(zf(-NGHOST:NZ+NGHOST))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build you custom variables !
CALL CUSTOM_HYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE BUILD_HYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
