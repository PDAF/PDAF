!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This variant is for online coupling with a model that is parallelized,
!! Thus MPI_init was called before calling this routine.
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
SUBROUTINE init_parallel_pdaf(screen, model_comm, model_comm_rank, model_comm_size)

  USE PDAF, &                     ! Command line parser
       ONLY: PDAF_parse, PDAF3_init_parallel
  USE mod_parallel_pdaf, &        ! PDAF parallelization variables
       ONLY: n_modeltasks, task_id, mype_world, npes_world, &
       mype_model, npes_model, COMM_model, &
       mype_filter, npes_filter, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)    :: screen           !< Whether screen information is shown

  ! Model parallelization variables (one can keep these generic names)
  INTEGER, INTENT(inout) :: model_comm       !< Model MPI communicator for model tasks
  INTEGER, INTENT(inout) :: model_comm_size  !< Number of processes in model_comm
  INTEGER, INTENT(inout) :: model_comm_rank  !< Process rank in model_comm

! *** Local variables ***
  INTEGER :: dim_ens                         ! Ensemble size
  CHARACTER(len=32) :: handle                ! Handle for command line parser


  ! Parse ensemble size
  handle = 'dim_ens'
  CALL PDAF_parse(handle, dim_ens)

  ! Set number of model tasks for fully-parallel mode
  n_modeltasks = 1

  ! Parse number of model tasks for flexible-parallel mode
  handle = 'n_tasks'
  CALL PDAF_parse(handle, n_modeltasks)


  ! For a parallelized model, MPI is initialized at this point
  ! Thus, dtore the parallelization variables provided by the model
  ! At this point, they describe all processes doing model integrations
  mype_world = model_comm_rank
  npes_world  = model_comm_size

  ! Initialize ensemble parallelization
  CALL PDAF3_init_parallel(screen, 0, 1, dim_ens, n_modeltasks, &
     model_comm, model_comm_rank, model_comm_size, &
     COMM_filter, mype_filter, npes_filter, &
     task_id)

  ! Initialize parallelization variables for parallel_pdaf_mod
  COMM_model = model_comm
  mype_model = model_comm_rank
  npes_model = model_comm_size

END SUBROUTINE init_parallel_pdaf
