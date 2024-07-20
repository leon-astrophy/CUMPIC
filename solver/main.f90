

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! CUMPIC3D - three-dimensional, newtonian hydrodynamic code 
! Fully parallized by MPI, support multi-GPU computing
! Written by Leung Shing Chi in 2016, major update by H.S. Leon Chan
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM CUMC3D
USE DEFINITION
IMPLICIT NONE

integer :: n

! Include MPI only if MPI is used !
#ifdef MPI
include "mpif.h" 
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI initialization

! Initialize MPI if MPI is used !
#ifdef MPI
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
#else
  mpi_rank = 0
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print only for the master rank !
IF(mpi_rank == 0) THEN 
  WRITE(*,*)
  WRITE(*,*) '---------------------'
  WRITE(*,*) '-Welcome to CUMPIC3D-'
  WRITE(*,*) '---------------------'
  WRITE(*,*)
  WRITE(*,*) '-----------------------------------------'
  WRITE(*,*) '- First written by Ka Wing Wong in 2010 -'
  WRITE(*,*) '- Updated by Shing Chi Leung in 2016    -'
  WRITE(*,*) '- Rewritten by H.S. Leon Chan in 2022   -'
  WRITE(*,*) '-----------------------------------------'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Setup the number of equations !
IF(mpi_rank == 0) THEN 
  WRITE(*,*) 'Build hydro equations'
END IF
CALL SETUP_EQNS
IF(mpi_rank == 0) THEN 
  WRITE(*,*) 'End building hydro equations'
  WRITE(*,*)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize MPI !
#ifdef MPI  
  CALL SETUP_MPI
#else
  starts(1) = 0
  starts(2) = 0
  starts(3) = 0
  stops(1) = starts(1) + nx
  stops(2) = starts(2) + ny
  stops(3) = starts(3) + nz
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for GPU !

#ifdef GPU
CALL ASSIGN_DEVICE
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial settings !

! Initialise !
n_step = 0
n_iter = 0
global_time = 0.0D0

! setup initial conditions !
CALL initial_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for GPU !

#ifdef GPU
CALL POPULATE_DEVICE
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! do initial updates !
CALL initial_update

! print primitive profiles !
CALL print_hydroprofile
n_iter = n_iter + 1

! Custom quantities !
CALL OPENFILE_CUSTOM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main rungekutta loop !

! Print out !
IF(mpi_rank == 0) THEN 
  WRITE (*,*) 'We now start the time evolution'
  WRITE (*,*)
  WRITE (*,*) '------------------------------------'
  WRITE (*,*) 'iteration, timestep, simulation time'
  WRITE (*,*) '------------------------------------'
END If

! Loop !
DO while (global_time < total_time)

  ! Find time step !
  CALL finddt

  ! Adjust dt for the end of simulation !
  IF(global_time + dt > total_time) THEN
    dt = total_time - global_time
  END IF

  ! Rungekutta step ! 
  CALL Rungekutta

  ! Update !
  n_step = n_step + 1
  global_time = global_time + dt

  ! Print out !
  IF(mpi_rank == 0) THEN 
    WRITE (*,*) n_step, dt, global_time
  END If

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!s	
  ! Section for data I/O !
  if(ABS(global_time - output_profiletime_last) >= output_profiletime) THEN

    ! print out !
	  output_profiletime_last = global_time
    call print_hydroprofile 
    n_iter = n_iter + 1

  END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output again

call print_hydroprofile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finalize MPI !

#ifdef MPI  
  CALL FINAL_MPI
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print out message 

! Print only for the master rank !
IF(mpi_rank == 0) THEN 
  WRITE(*,*) '----------------'
  WRITE(*,*) '-Simulation End-'
  WRITE(*,*) '----------------'
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM
