!$Id$
!> Module for ensemble parallelization
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. The are variables
!! that are used in the model, even without PDAF and additional
!! variables that are only used, if data assimilation with PDAF
!! is performed.
!!
!! In addition methods to initialize and finalize MPI are provided.
!! The initialization routine is only for the model itself, the 
!! more complex initialization of communicators for execution with
!! PDAF is peformed in init_parallel_pdaf.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_parallel_pdaf

  USE mpi

  IMPLICIT NONE
  SAVE 

  ! Basic variables for model state integrations
  INTEGER :: COMM_model  !< MPI communicator for model tasks
  INTEGER :: mype_model  !< Number of PEs in COMM_model
  INTEGER :: npes_model  !< PE rank in COMM_model

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks = 1         !< Number of parallel model tasks
  INTEGER :: n_filterpes  = 1         !< Number of PEs for filter analysis
  INTEGER :: npes_world               !< Number of processes in MPI_COMM_WORLD
  INTEGER :: mype_world               !< Process rank in MPI_COMM_WORLD
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

CONTAINS
!-------------------------------------------------------------------------------
!> Initialize MPI
!!
!! Routine to initialize MPI, the number of PEs
!! (npes_world) and the rank of a PE (mype_world).
!! The model is executed within the scope of the
!! communicator Comm_model. It is also initialized
!! here together with its size (npes_model) and 
!! the rank of a PE (mype_model) within Comm_model.
!!
  SUBROUTINE init_parallel()

    IMPLICIT NONE

    INTEGER :: i
  
    CALL MPI_INIT(i);
    CALL MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for MPI_COMM_WORLD
    Comm_model = MPI_COMM_WORLD
    npes_model = npes_world
    mype_model = mype_world
   
  END SUBROUTINE init_parallel
!-------------------------------------------------------------------------------
!> Finalize MPI
!!
!! Routine to finalize MPI
!!
  SUBROUTINE finalize_parallel()

    IMPLICIT NONE
    
    CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    CALL  MPI_Finalize(MPIerr)

  END SUBROUTINE finalize_parallel
!-------------------------------------------------------------------------------
!> Abort MPI
!!
!! Routine for abort MPI program.
!!
  SUBROUTINE abort_parallel()

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_pdaf
