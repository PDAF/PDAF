!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! The program shows the setup for the fully-parallel
!! implementation variant of PDAF with a parallelized
!! model.
!! 
!! In the online implementation with a real model
!! this driver program is replaced by the actual
!! model code.
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code
!!
PROGRAM MAIN

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: n_modeltasks, npes_world, mype_world
  USE mod_model, &             ! Module provided by model code
       ONLY: step_final        

  IMPLICIT NONE

! local variables
  INTEGER :: istep       ! Counter


! *** Initialize MPI ***

  ! If the model itself is parallelized this step is done by the model
  ! The initialization of ensemble-parallelization is added to this routine

  CALL init_parallel() ! initializes MPI

  ! FOR TESTING: 
  n_modeltasks = npes_world

  IF (mype_world==0) THEN
     WRITE (*,*) '**********************************************************************'
     WRITE (*,*) '*   THIS IS A TEST PROGRAM TO CHECK THE TEMPLATE CODE CONSISTENCY    *'
     WRITE (*,*) '*                   Run this program with:                           *'
     WRITE (*,*) '*          mpirun -np NENS ./PDAF_online -dim_ens NENS               *'
     WRITE (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
     WRITE (*,*) '**********************************************************************'
  END IF

  
  ! MODEL: Here the model would perform its initialization


! *** Initialize PDAF ***

  ! This step is always inserted after the model initialization
  ! is complete and just before the time stepping starts

  CALL init_pdaf()


! *** Ensemble forecasting and analysis steps ***

  ! MODEL: In the real model this is the time stepping loop of the model
  timesteps: DO istep = 1, step_final

     ! MODEL: Here the model code would compute the time stepping


     ! *** Perform analysis ***

     ! This step is inserted in the time stepping loop
     ! usually just before the 'end do'

     CALL assimilate_pdaf()

  ENDDO timesteps


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
