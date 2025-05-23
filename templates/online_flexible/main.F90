!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! This variant is for the flexible parallelization
!! variant of PDAF using PDAF3_assimilate routines. 
!! It shows the structure of the required outper loop
!! which enables to integrate an ensemble of model states.
!! 
!! In the online implementation with a real model
!! The outer loop and control structure would be
!! inserted in the actual model code.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code using PDAF3_assimilate
!!
PROGRAM MAIN

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: init_parallel, finalize_parallel, &
       n_modeltasks, mype_world
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: time
  USE mod_model, &             ! Module provided by model code
       ONLY: dt
  USE PDAF , &                 ! Interface definitions to PDAF core routines
       ONLY: PDAF_init_forecast, PDAF_get_fcst_info

  IMPLICIT NONE

! local variables
  INTEGER :: istep       ! Counter
  INTEGER :: nsteps      ! Number of time steps to be performed in current forecast
  INTEGER :: doexit      ! Whether to exit forecasting (1=true)
  INTEGER :: status_pdaf ! PDAF status flag      
  REAL :: timenow        ! Current model time

! ! External subroutines 
  EXTERNAL :: distribute_state_pdaf, &  ! Distribute a state vector to model fields
       prepoststep_pdaf, &              ! User supplied pre/poststep routine
       next_observation_pdaf            ! Provide time step of next observation

! *** Initialize MPI ***

  ! If the model itself is parallelized this step is done by the model

  CALL init_parallel() ! initializes MPI

  ! FOR TESTING: 
  n_modeltasks = 1

  IF (mype_world==0) THEN
     WRITE (*,*) '**********************************************************************'
     WRITE (*,*) '*   THIS IS A TEST PROGRAM TO CHECK THE TEMPLATE CODE CONSISTENCY    *'
     WRITE (*,*) '*                   Run this program with:                           *'
     WRITE (*,*) '*                ./PDAF_online -dim_ens NENS                         *'
     WRITE (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
     WRITE (*,*) '**********************************************************************'
  END IF

  
! *** Initialize MPI communicators for PDAF (model, filter and coupling) ***

  ! This step is always inserted directly after the MPI initialization

  CALL init_parallel_pdaf(0, 1)

  
  ! MODEL: Here the model would perform its initialization


! *** Initialize PDAF ***

  ! This step is always inserted after the model initialization
  ! is complete and just before the time stepping starts

  CALL init_pdaf(nsteps, timenow, doexit)


! *** PDAF: Get state and forecast information (nsteps,time) at initial time ***

!  CALL PDAF_init_forecast(nsteps, timenow, doexit, next_observation_pdaf, &
!       distribute_state_pdaf, prepoststep_pdaf, status_pdaf)

! *** Ensemble forecasting and analysis steps ***

  ! PDAF: External loop around model time stepper loop
  pdaf_modelloop: DO  

     ! *** Forecast ensemble state ***
 
     ! Initialize current model time
     time = timenow

     ! *** run time stepper ***  

     ! MODEL: Here the model code would do the time stepping
     DO istep = 1, nsteps
        WRITE (*,'(3x, a, f6.2)') 'main.F90: Do stepping, time', time

        ! The model would increment the time information
        time = time + dt  

        ! *** Let PDAF check forecast progress and perform analysis ***
        CALL assimilate_pdaf()

     ENDDO

     ! *** PDAF: Get forecast information ***
     CALL PDAF_get_fcst_info(nsteps, timenow, doexit)

     ! *** Check exit flag ***
     IF (doexit==1) EXIT pdaf_modelloop

  END DO pdaf_modelloop



! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  CALL finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  CALL finalize_parallel()

END PROGRAM MAIN
