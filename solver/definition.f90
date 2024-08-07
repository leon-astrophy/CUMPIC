!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! definition.f90
! contains key simulation parameteres, define global variables, and declare arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DEFINITION
USE HDF5
IMPLICIT NONE
SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUMPIC3D (Last Modified: Leon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Foreword:
! Developed by Leung Shing Chi based on the WENO prototype developed by Wong Ka Wing in 2010
! c.f. Leung et al., MNRAS 454, 1238 (2015).
! Further developed by Leon H.S. Chan to extend to 3D, included MHD, parallelized by MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control work distribution across CPUs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of CPUs along the first dimension
INTEGER, PARAMETER :: NXCPU = 1

! Number of CPUs along the second dimension
INTEGER, PARAMETER :: NYCPU = 1

! Number of CPUs along the third dimension
INTEGER, PARAMETER :: NZCPU = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Constant and universal values 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Boundary flag notation for boundary conditions
INTEGER, PARAMETER :: even = 0, odd = 1

! Which direction of reconstruction
INTEGER, PARAMETER :: x_dir = 1, y_dir = 2, z_dir = 3

! pi constant 
REAL*8, PARAMETER :: pi = 4.D0*DATAN(1.D0)

! define small number to avoid coordinate singularity !
REAL*8, PARAMETER :: small_num = TINY(1.0D0)

! Number of ghost shell !
INTEGER, PARAMETER :: NGHOST = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Core part of the simulation box
! This is where the system variables are controled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Define coordinate here !
INTEGER, PARAMETER :: cartesian = 1
INTEGER, PARAMETER :: cylindrical = 2
INTEGER, PARAMETER :: spherical = 3

! Use which coordiante? !
INTEGER, PARAMETER :: coordinate = cylindrical

! Flag for boundary condition
! The boundary flag is defined by four scalar
! 1st one for x-inner boundary
! 2nd one for x-outer boundary
! 3rd one for y-inner boundary
! 4th one for y-outer boundary
! 5rd one for z-inner boundary
! 6th one for z-outer boundary
! 0 = periodic
! 1 = outgoing (1st derivative = 0)
! 2 = reflecting boundary (depend on scalar/vector)
! 3 = axis-symmetric
! 4 = equatorial-symmetric
INTEGER, PARAMETER :: boundary_flag(6) = (/1,1,0,0,1,1/)

! Starting position of the grid !
REAL*8, PARAMETER :: x_start = 1.5d0
REAL*8, PARAMETER :: y_start = 0.0d0
REAL*8, PARAMETER :: z_start = -5.0d0

! Ending position of the grid !
REAL*8, PARAMETER :: x_end = 11.5d0
REAL*8, PARAMETER :: y_end = 2.0d0*pi
REAL*8, PARAMETER :: z_end = 5.0d0

! The total number of grid in the x, y, z direction
INTEGER, PARAMETER :: nxtot = 128
INTEGER, PARAMETER :: nytot = 128
INTEGER, PARAMETER :: nztot = 128

! Grid sizes for uniform grid 
REAL*8, PARAMETER :: dx = (x_end - x_start)/DBLE(nxtot)	
REAL*8, PARAMETER :: dy = (y_end - y_start)/DBLE(nytot)	
REAL*8, PARAMETER :: dz = (z_end - z_start)/DBLE(nztot)	

! Number of grid per MPI processes !
INTEGER, PARAMETER :: nx = nxtot/NXCPU
INTEGER, PARAMETER :: ny = nytot/NYCPU
INTEGER, PARAMETER :: nz = nztot/NZCPU

! Cournat-Friedrich-Levy constant
! Defined as dt = cfl * dx / MAX(vel + cs)
REAL*8, PARAMETER :: cfl = 0.30D0			

! Maximum time to be simulated in the model
REAL*8, PARAMETER :: total_time = 1000.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for Output setting
! This sets how frequent each type of profile is output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Physical time interval for all log file
REAL*8 :: output_logtime_last = 0.0D0
REAL*8, PARAMETER :: output_logtime = 1.0D1                            

! Physical time interval for each hydro profile
REAL*8 :: output_profiletime_last = 0.0D0
REAL*8, PARAMETER :: output_profiletime = 0.05d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for riemann solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the LF (two-shocks) Riemann solver !
INTEGER, PARAMETER :: LF = 1

! Use the HLL Riemann solver !
INTEGER, PARAMETER :: HLL = 2

! Use the HLLC Riemann solver !
INTEGER, PARAMETER :: HLLC = 3

! Use the HLLD Riemann solver !
INTEGER, PARAMETER :: HLLD = 4

! Define reconstruction method #
INTEGER, PARAMETER :: SOLVER = HLL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for primitive reconstructions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the TVD (Mignone 2014) reconstruction scheme with Min-Mod limiter !
INTEGER, PARAMETER :: TVDMM = 1

! Use the TVD (Mignone 2014) reconstruction scheme with MC limiter !
INTEGER, PARAMETER :: TVDMC = 2

! Use the TVD (Mignone 2014) reconstruction scheme with Van-Leer limiter !
INTEGER, PARAMETER :: TVDVL = 3

! Use the PPM (Colella 1984) reconstruction scheme !
INTEGER, PARAMETER :: PPMC = 4

! Use the WENO (Shu 1997) reconstruction scheme !
INTEGER, PARAMETER :: WENO = 5

! Define reconstruction method #
INTEGER, PARAMETER :: RECON = TVDMM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for rk3 time evolution 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Storing coefficients for the RK method !
REAL*8, PARAMETER :: rk20 = 3.0d0/4.0d0
REAL*8, PARAMETER :: rk21 = 1.0d0/4.0d0
REAL*8, PARAMETER :: rk22 = 1.0d0/4.0d0
REAL*8, PARAMETER :: rk30 = 1.0d0/3.0d0
REAL*8, PARAMETER :: rk31 = 2.0d0/3.0d0
REAL*8, PARAMETER :: rk32 = 2.0d0/3.0d0

!----------------------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Equations identifiers !

! Min/max number of equations
INTEGER :: imin
INTEGER :: imax

! Total number of equations !
INTEGER :: no_of_eq

! density 
INTEGER :: irho

! x-velocity
INTEGER :: ivx

! y-velocity
INTEGER :: ivy

! z-velocity
INTEGER :: ivz

! total energy density
INTEGER :: itau

! magnetic fields !
INTEGER :: ibx, iby, ibz

!=====================================================================================================!
! Boundary flag for variables !

INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_xin
INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_yin
INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_zin
INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_xout
INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_yout
INTEGER, ALLOCATABLE, DIMENSION (:) :: bfac_zout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gloabl indexes for array loop !

! start and ending index !
INTEGER :: starts (3)
INTEGER :: stops (3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are the hydrodynamics 

! primitive and conservative variables !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: prim
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: cons

! Primitive variables at the L/R of cell boundaries 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: primL
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: primR

! Conservative variables at the L/R of cell boundaries 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: consL
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: consR

! Fluxes at the L/R of cell boundaries 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: fluxL
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: fluxR

! Epsilon !
REAL*8, ALLOCATABLE, DIMENSION (:,:,:) :: eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for time-evolution !

! For RK-Time evolution 
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: l_rk
REAL*8, ALLOCATABLE, DIMENSION (:,:,:,:) :: u_old

! source terms !
REAL*8, allocatable, DIMENSION (:,:,:,:) :: sc
REAL*8, allocatable, DIMENSION (:,:,:,:) :: flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for MHD !

! cell centered magnetic fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: bcell

! cell centered electric fields !
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecell

! cell cornered electrin fields
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecorn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computational Grid !

! face coordiante of x !
REAL*8, ALLOCATABLE, DIMENSION(:) :: xf

! face coordiante of y !
REAL*8, ALLOCATABLE, DIMENSION(:) :: yf

! face coordiante of z !
REAL*8, ALLOCATABLE, DIMENSION(:) :: zf

!=====================================================================================================!
! MHD identifiers !

! Indices for cell cornered and centered electric fields !
INTEGER, PARAMETER :: iex = 1
INTEGER, PARAMETER :: iey = 2
INTEGER, PARAMETER :: iez = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Real/integer variable scalars !

! Time step !
REAL*8 :: dt

! RK step number !
INTEGER :: n_step

! File counter !
INTEGER :: n_iter

! Global simulation time !
REAL*8 :: global_time

!****************************************************************************************************!
! MPI stuff !

! Integer for MPI stuff !
INTEGER :: mpi_rank, mpi_size, ierror, numcpus

! For cartesian topology !
INTEGER :: new_comm

! For communication !
INTEGER, DIMENSION(1:3) :: face_type
INTEGER, DIMENSION(1:3) :: bcell_type
INTEGER, DIMENSION(1:3) :: bface_type
INTEGER, DIMENSION(1:3) :: scalar_type
INTEGER, DIMENSION(0:2,0:2,0:2) :: neighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HDF5 stuff !

! for HDF5 !
integer :: error
integer(HID_T) :: file_id, dset_id, plist_id, space_id, mem_id

!****************************************************************************************************!

! Section for GPU !
#ifdef GPU
!$ACC declare create(no_of_eq)
!$ACC declare create(imin, imax) 
!$ACC declare create(irho, itau)
!$ACC declare create(ivx, ivy, ivz)
!$ACC declare create(ibx, iby, ibz)
!$ACC declare create(xF, yF, zF)
!$ACC declare create(bcell, prim, sc)
!$ACC declare create(primL, primR, flux)
!$ACC declare create(consL, consR, eps)
!$ACC declare create(prim, cons, bcell)
!$ACC declare create(fluxL, fluxR)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing cartesian dot_product
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! dot product between two vectors !
	REAL*8 function dot_product(vec_a, vec_b)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8, DIMENSION (1:3) :: vec_a, vec_b
	dot_product = vec_a(1)*vec_b(1) + vec_a(2)*vec_b(2) + vec_a(3)*vec_b(3)
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for computing the HLL flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	REAL*8 function compute_hll(ul,ur,fl,fr,sl,sr)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: ul, ur, fl, fr, sl, sr
	compute_hll = (sr*ur - sl*ul - (fr - fl))/(sr - sl)
	end function
	
	REAL*8 function compute_fluxhll(yl,yr,xl,xr,sl,sr)	
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: yl, yr, xl, xr, sl, sr
	compute_fluxhll = (sr*yl-sl*yr+sl*sr*(xr-xl))/(sr-sl)
	end function
	
	REAL*8 function compute_roe(xl,xr,rhol,rhor)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: xl,xr,rhol,rhor
	compute_roe = (DSQRT(rhol)*xl + DSQRT(rhor)*xr)/(DSQRT(rhol) + DSQRT(rhor))
	end function
	
	REAL*8 function compute_signalspeed(cs, bn, bt1, bt2, rho)
	!$ACC ROUTINE SEQ
	implicit none
	REAL*8 :: cs, bn, bt1, bt2, rho
	REAL*8 :: a2_mhd, b2_mhd, a4_mhd, b4_mhd
	REAL*8 :: b2n_mhd, b2t1_mhd, b2t2_mhd
	a2_mhd = cs*cs
	a4_mhd = a2_mhd*a2_mhd
	b2n_mhd = (bn*bn/rho)
	b2t1_mhd = (bt1*bt1/rho)
	b2t2_mhd = (bt2*bt2/rho)
	b2_mhd = b2n_mhd + b2t1_mhd + b2t2_mhd
	b4_mhd = b2_mhd*b2_mhd
	compute_signalspeed = DSQRT(0.5D0*(a2_mhd + b2_mhd + & 
												DSQRT(MAX((a4_mhd + b4_mhd + 2.0d0*b2_mhd*a2_mhd) & 
												- 4.0D0*a2_mhd*b2n_mhd, 0.0d0))))
	end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
END MODULE DEFINITION
