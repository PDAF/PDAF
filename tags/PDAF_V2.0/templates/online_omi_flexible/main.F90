!$Id: main_offline.F90 871 2021-11-22 16:44:34Z lnerger $
!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! This variant is for the flexible parallelization
!! variant of PDAF. It shows the structure of the
!! required outper loop which enables to integrate
!! an ensemble of model states.
!! 
!! In the online implementation with a real model
!! this driver program is replaced by the actual
!! model code.
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code
!!
PROGRAM MAIN

  USE mpi                        ! MPI
  USE mod_parallel_pdaf, &       ! Parallelization
       ONLY: init_parallel, finalize_parallel, &
       n_modeltasks, mype_world
  USE mod_assimilation, &        ! Assimilation variables
       ONLY: time
  USE mod_model, &               ! Module provided by model code
       ONLY: dt
  USE pdaf_interfaces_module, &  ! Interface definitions to PDAF core routines
       ONLY: PDAF_get_state

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

  CALL init_pdaf()


! *** Ensemble forecasting and analysis steps ***

  ! PDAF: External loop around model time stepper loop
  pdaf_modelloop: DO  

     ! *** PDAF: Get state and forecast information (nsteps,time)  ***
     CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
          distribute_state_pdaf, prepoststep_pdaf, status_pdaf)

     ! Check whether forecast has to be performed
     checkforecast: IF (doexit /= 1 .AND. status_pdaf == 0) THEN

        ! *** Forecast ensemble state ***
      
        IF (nsteps > 0) THEN

           ! Initialize current model time
           time = timenow

           ! *** call time stepper ***  

            ! MODEL: Here the model code would do the time stepping
           DO istep = 1, nsteps
              WRITE (*,'(3x, a, f6.2)') 'main.F90: Do stepping, time', time

              ! The model would increment the time information
              time = time + dt  
           ENDDO

        END IF

        ! *** Perform analysis ***

        ! This step is inserted after the inner time stepping loop

        CALL assimilate_pdaf()


     ELSE checkforecast

        ! *** No more work, exit modeling loop
        EXIT pdaf_modelloop

     END IF checkforecast

  END DO pdaf_modelloop



! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  CALL finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  CALL finalize_parallel()

END PROGRAM MAIN
