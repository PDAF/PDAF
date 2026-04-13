!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! This variant is for the flexible parallelization
!! variant of PDAF using PDAF3_assimilate routines. 
!! It shows the structure of the required outer loop
!! which enables to integrate an ensemble of model states.
!! 
!! In the online implementation with a real model
!! the outer loop and control structure would be
!! inserted in the actual model code.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code using PDAF3_assimilate
!!
PROGRAM MAIN

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: n_modeltasks, mype_world
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


! *** Initialize MPI ***

  ! If the model itself is parallelized this step is done by the model
  ! The initialization of ensemble-parallelization is added to this routine

  CALL init_parallel() ! initializes MPI

  ! FOR TESTING: 
  n_modeltasks = 1

  IF (mype_world==0) THEN
     WRITE (*,*) '**********************************************************************'
     WRITE (*,*) '*   THIS IS A TEST PROGRAM TO CHECK THE TEMPLATE CODE CONSISTENCY    *'
     WRITE (*,*) '*                   Run this program with:                           *'
     WRITE (*,*) '*                ./PDAF_online -dim_ens NENS                         *'
     WRITE (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
     WRITE (*,*) '*           Alternatively run this program with:                     *'
     WRITE (*,*) '*       mpirun -np NP ./PDAF_online -dim_ens NENS -n_tasks NTSK      *'
     WRITE (*,*) '* with number of processes NP and number of model tasks NTSK         *'
     WRITE (*,*) '* (required are NTSK<=NP, NTSK<=NENS)                                *' 
     WRITE (*,*) '**********************************************************************'
  END IF

  
  ! MODEL: Here the model would perform its initialization


! *** Initialize PDAF ***

  ! This step is always inserted after the model initialization
  ! is complete and just before the time stepping starts

  CALL init_pdaf()


! *** Ensemble forecasting and analysis steps ***

  ! PDAF: External loop around model time stepper loop
  pdaf_modelloop: DO  

     ! *** Forecast ensemble state ***

     ! *** PDAF: Get forecast information ***
     CALL PDAF_get_fcst_info(nsteps, timenow, doexit)

     ! *** Check exit flag ***
     IF (doexit==1) EXIT pdaf_modelloop
 
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

  END DO pdaf_modelloop



! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  CALL finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  CALL finalize_parallel()

END PROGRAM MAIN

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!! The two following subroutines are helpers to initialize
!! and finalize MPI. With a real parallellized model, 
!! this functionality would be in the model code.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE init_parallel()

  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: MPI_COMM_WORLD, npes_world, mype_world, &
       COMM_model, npes_model, mype_model
       
  IMPLICIT NONE

  INTEGER :: MPIerr
  
  CALL MPI_INIT(MPIerr);
  CALL MPI_Comm_Size(MPI_COMM_WORLD,npes_world,MPIerr)
  CALL MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,MPIerr)

  ! Initialize model communicator, its size and the process rank
  ! Here the same as for MPI_COMM_WORLD
  Comm_model = MPI_COMM_WORLD
  npes_model = npes_world
  mype_model = mype_world

  ! Initialize parallelization for PDAF
  ! This step is always inserted directly after the MPI initialization

  CALL init_parallel_pdaf(1, COMM_model, mype_model, npes_model)
   
END SUBROUTINE init_parallel
!-------------------------------------------------------------------------------
!> Finalize MPI
!!
!! Routine to finalize MPI
!!
SUBROUTINE finalize_parallel()

  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: MPI_COMM_WORLD

  IMPLICIT NONE
    
  INTEGER :: MPIerr

  CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  CALL  MPI_Finalize(MPIerr)

END SUBROUTINE finalize_parallel
