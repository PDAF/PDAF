!$Id$
!> Module for ensemble parallelization
!!
!! This modules provides variables for the MPI parallelization
!! to be shared between model-related routines. The are variables
!! that are used in the model, even without PDAF and additional
!! variables that are only used, if data assimialtion with PDAF
!! is performed.
!!
!! Revision history:
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_parallel_pdaf

  USE mpi

  IMPLICIT NONE
  SAVE 

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks = 1         !< Number of parallel model tasks
  INTEGER :: n_filterpes  = 1         !< Number of PEs for filter analysis
  INTEGER :: COMM_filter              !< MPI communicator for filter PEs 
  INTEGER :: mype_filter              !< Number of processes in COMM_filter
  INTEGER :: npes_filter              !< Process rank in COMM_filter
  INTEGER :: COMM_couple              !< MPI communicator for coupling filter and model
  LOGICAL :: modelpe                  !< Whether we are on a PE in a COMM_model
  LOGICAL :: filterpe                 !< Whether we are on a PE in a COMM_filter
  INTEGER :: task_id                  !< Index of my model task (1,...,n_modeltasks)
  INTEGER :: MPIerr                   !< Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       !< Status array for MPI
  INTEGER, ALLOCATABLE :: local_npes_model(:) !< Number of processes per ensemble

END MODULE mod_parallel_pdaf
