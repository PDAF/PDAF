!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This variant is for online coupling with a model that is not parallelized,
!! Thus there are no MPI operations in the model code.
!!
!! This routine stores the parallelization information provided by the
!! model. Afterwards, it calls the parallelization routine for PDAF, which 
!! initializes the communicators for handling the ensemble and the
!! the analysis of the data assimilation. This overwrites the communicator
!! provided by the model by the ensemble configuration.
!!
!! The parallelization variables returned from the PDAF initialization
!! routines are stored in variables from the module parallel_pdaf_mod
!! so that they can be used in the different user-provided routines.
!!
!! The routine is generic, but has to be compiled with the user code
!! because it uses the module parallel_pdaf_mod.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * 2026-02 - Lars Nerger - Revision for using PDAF3_init_forecast 
!! * Later revisions - see repository log
!!
SUBROUTINE init_parallel_pdaf(screen)

  USE PDAF, &                     ! Command line parser
       ONLY: PDAF_parse, PDAF3_init_parallel
  USE mod_parallel_pdaf, &        ! PDAF parallelization variables
       ONLY: n_modeltasks, task_id, &
       mype_world, npes_world, MPI_COMM_WORLD, &
       mype_model, npes_model, COMM_model, &
       mype_filter, npes_filter, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)    :: screen           !< Whether screen information is shown

! *** Local variables ***
  INTEGER :: dim_ens                         ! Ensemble size
  CHARACTER(len=32) :: handle                ! Handle for command line parser
  INTEGER :: MPIerr                          ! MPI error flag


  ! Parse ensemble size
  handle = 'dim_ens'
  CALL PDAF_parse(handle, dim_ens)

  ! Set number of model tasks for fully-parallel mode
  n_modeltasks = dim_ens

  ! Initialize ensemble parallelization
  CALL PDAF3_init_parallel(screen, 0, 1, dim_ens, n_modeltasks, &
     COMM_model, mype_model, npes_model, &
     COMM_filter, mype_filter, npes_filter, &
     task_id)

  ! *** Initialize PE information on COMM_world ***
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)

END SUBROUTINE init_parallel_pdaf
