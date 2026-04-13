!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This variant is for online coupling with a model that is parallelized,
!! Thus MPI_init was called before calling this routine. The parallelization
!! The parallelization configuration used here is that there are dim_ens
!! model tasks and a separate task that exclusively computes the analysis
!! step, with the same number of processes as a model task.
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
!! * 2014-04 - Lars Nerger - Variant for separate processes for model and filter
!! * 2026-02 - Lars Nerger - Revision for using PDAF3_init_forecast 
!! * Other revisions - see repository log
!!
SUBROUTINE init_parallel_pdaf(screen)

  USE mpi                         ! MPI
  USE PDAF, &                     ! PDAF routines
       ONLY: PDAF_parse, PDAF3_set_parallel, PDAF3_init_parallel
  USE mod_parallel_pdaf, &        ! PDAF parallelization variables
       ONLY: mype_filter, npes_filter, COMM_filter, &
       n_modeltasks, task_id
  USE mod_parallel_model, &       ! Model parallelization variables
       ONLY: mype_world, npes_world, mype_2dmodel, npes_2dmodel, &
       COMM_2dmodel, MPIerr, modelproc

  IMPLICIT NONE    
  
! *** Arguments ***
  INTEGER, INTENT(in) :: screen       !< Whether screen information is shown

! *** local variables ***
  INTEGER :: i, j                     ! Counters
  INTEGER :: dim_ens                  ! Ensemble size
  INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
  INTEGER :: pe_index                 ! Index of PE
  INTEGER :: my_color, color_couple   ! Variables for communicator-splitting 
  LOGICAL :: iniflag                  ! Flag whether MPI is initialized
  INTEGER :: flag                     ! Status flag
  INTEGER, ALLOCATABLE :: local_npes_model(:) ! Array for npes per model task
  INTEGER :: COMM_couple              ! MPI communicator for coupling filter and model
  CHARACTER(len=32) :: handle         ! handle for command line parser
  INTEGER :: t_id                     ! Variable for storing task id for communicator splitting


  ! Parse ensemble size
  handle = 'dim_ens'
  CALL PDAF_parse(handle, dim_ens)

  ! Set number of model tasks for fully-parallel mode
  n_modeltasks = dim_ens

  ! Initialize ensemble parallelization
  call PDAF3_init_parallel(1, 1, 1, dim_ens, n_modeltasks, &
     comm_2dmodel, mype_2dmodel, npes_2dmodel, &
     COMM_filter, mype_filter, npes_filter, &
     task_id)

  ! Set flag whether a process computes the model integrations
  IF (comm_2dmodel/= MPI_COMM_NULL) THEN
     modelproc = .TRUE.
  ELSE
     modelproc = .FALSE.
  END IF

END SUBROUTINE init_parallel_pdaf
