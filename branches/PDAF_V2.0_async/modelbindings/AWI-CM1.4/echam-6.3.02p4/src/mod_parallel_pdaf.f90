!$Id$
!>  Module for MPI parallelization of the ensemble 
!!
!! This modules provides variables for the MPI parallelization
!! to be shared between model-related routines. The are variables
!! that are used in the model, even without PDAF and additional
!! variables that are only used, if data assimialtion with PDAF
!! is performed.
!! In addition methods to initialize and finalize MPI are provided.
!! The initialization routine is only for the model itself, the 
!! more complex initialization of communicators for xecution with
!! PDAF is performed in init_parallel_pdaf.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
MODULE mod_parallel_pdaf

  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

  ! Basic variables for model state integrations
  INTEGER :: COMM_model  ! MPI communicator for model tasks
  INTEGER :: mype_model  ! Number of PEs in COMM_model
  INTEGER :: npes_model  ! PE rank in COMM_model
  INTEGER :: mype_world  ! Number of PEs in MPI_COMM_WORLD
  INTEGER :: npes_world  ! PE rank in MPI_COMM_WORLD
  INTEGER :: mype_submodel  ! Number of PEs in model compartment task

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks = 1         ! Number of parallel model tasks
  INTEGER :: n_filterpes  = 1         ! Number of PEs for filter analysis
  INTEGER :: COMM_filter ! MPI communicator for filter PEs 
  INTEGER :: mype_filter, npes_filter ! # PEs and PE rank in COMM_filter
  INTEGER :: COMM_couple ! MPI communicator for coupling filter and model
  LOGICAL :: modelpe     ! Whether we are on a PE in a COMM_model
  LOGICAL :: filterpe    ! Whether we are on a PE in a COMM_filter
  LOGICAL :: pairs       ! Whether we use pairs of fesom.x and echam.x
  INTEGER :: task_id     ! Index of my model task (1,...,n_modeltasks)
  INTEGER :: MPIerr      ! Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  INTEGER, ALLOCATABLE :: local_npes_model(:) ! # PEs per ensemble
  LOGICAL :: writepe     ! Whether the process does file writing

  ! Specific variables for atmosphere compartment
  INTEGER :: COMM_filter_echam ! MPI communicator for filter PEs - only ECHAM part
  INTEGER :: mype_filter_echam, npes_filter_echam ! # PEs and PE rank in COMM_filter - only ECHAM part

!-------------------------------------------------------------------------------

CONTAINS

!> Routine to abort MPI program
!!
  SUBROUTINE abort_parallel()

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_pdaf
