!BOP
!
! !ROUTINE: integration_pdaf  --- integration routine for L96 with assimilation
!
! !INTERFACE:
SUBROUTINE integration_pdaf(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the Lorenz96 model.
!
! For simplicity of the implementation with PDAF,
! the time stepping is separated into a single routine.
! This allows to simply put the assimilation routine
! assimilation\_pdaf() in between the main program and
! the integration routine. If the time stepping is not
! available as a separate routine, a different 
! implementation style is required.
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! 2025-03 - Lars Nerger - Adaption to using PDAF_assimilate with PDAF3
! Later revisions - see svn log
!
! !USES:
  USE PDAF, &
       ONLY: PDAF_get_state, PDAF_get_fcst_info
  
  USE timer, &            ! Timing
       ONLY: timeit      
  USE mod_model, &        ! Model variables
       ONLY: x, dt, dim_state
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: model_error
  USE output_netcdf, &    ! NetCDF output
       ONLY: write_netcdf, close_netcdf

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(inout) :: time   ! Model time
  INTEGER, INTENT(inout) :: nsteps ! Number of time steps to be performed
!EOP

! local variables
  INTEGER :: step               ! Time step counter
  REAL, ALLOCATABLE :: x1(:), x2(:), x3(:), x4(:) ! Temporary arrays for RK4
  REAL, ALLOCATABLE :: x_tmp(:)
! Variables for PDAF
  INTEGER :: doexit      ! Whether to exit forecasting (1=true)
  INTEGER :: status      ! Status flag for filter routines
  REAL :: timenow        ! Current model time

! ! External subroutines 
! Interface between model and PDAF, and prepoststep
  EXTERNAL :: distribute_state_pdaf, & ! Distribute a state vector to model fields
       next_observation_pdaf, &        ! Provide time step of next observation
       prepoststep_pdaf                ! User supplied pre/poststep routine


! **********************
! *** Initialization ***
! **********************

  ! Allocate temporary arrays for RK4
  ALLOCATE(x1(dim_state))
  ALLOCATE(x2(dim_state))
  ALLOCATE(x3(dim_state))
  ALLOCATE(x4(dim_state))
  ALLOCATE(x_tmp(dim_state))


! *************************************************
! *** Initialize forecast for data assimilation ***
! *************************************************

!++ Addition for PDAF     
  ! *** PDAF: Get state and forecast information (nsteps,time) at initial time ***
  CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_pdaf, status)

  ! PDAF: External loop around model time stepper loop
  pdaf_modelloop: DO  
!++ End addition for PDAF     


! *********************************
! *** Perform model integration ***
! *********************************

     CALL timeit(5, 'new')

     ! *** time stepping loop ***
     integrate: DO step = 1, nsteps

        ! Intermediate steps
        CALL lorenz96_dxdt(dim_state, x, x1)
        x1 = dt * x1

        x_tmp = x + x1/2.0
        CALL lorenz96_dxdt(dim_state, x_tmp, x2)
        x2 = dt * x2

        x_tmp = x + x2/2.0
        CALL lorenz96_dxdt(dim_state, x_tmp, x3)
        x3 = dt * x3

        x_tmp = x + x3
        CALL lorenz96_dxdt(dim_state, x_tmp, x4)
        x4 = dt * x4

        ! New value of x
        x = x + x1/6.0 + x2/3.0 + x3/3.0 + x4/6.0

        ! Increment time
        time = time + dt

!++ Addition for PDAF     
        ! *** PDAF: Add model error ***      
        IF (model_error) CALL add_model_noise(dt, dim_state, x)

        ! *** PDAF: check forecast progress and perform analysis ***
        CALL assimilate_pdaf()
!++ End addition for PDAF     

     END DO integrate

!++ Addition for PDAF     
     ! *** PDAF: Get state and forecast information (nsteps,time)  ***
     CALL PDAF_get_fcst_info(nsteps, timenow, doexit)

     ! *** Check exit flag ***
     IF (doexit==1) EXIT pdaf_modelloop
!++ End addition for PDAF     

  END DO pdaf_modelloop

  DEALLOCATE(x1, x2, x3, x4, x_tmp)

  CALL timeit(5, 'old')

END SUBROUTINE integration_pdaf
