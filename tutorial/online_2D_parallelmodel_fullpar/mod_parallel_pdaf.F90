!> Module for ensemble parallelization
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. The are variables
!! that are used in the model, even without PDAF and additional
!! variables that are only used, if data assimilation with PDAF
!! is performed.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_parallel_pdaf

  USE mpi

  IMPLICIT NONE
  SAVE 

  ! Parallelization variables that can be used in the user code

  ! Variables for each model task
  INTEGER :: COMM_model=0               !< MPI communicator for model tasks
  INTEGER :: mype_model=0               !< PE rank in COMM_model
  INTEGER :: npes_model=1               !< Number of PEs in COMM_model

  ! Variables describing all processes involved in model integrations
  INTEGER :: mype_world=0               !< Process rank in MPI_COMM_WORLD
  INTEGER :: npes_world=1               !< Number of processes in MPI_COMM_WORLD

  ! Variables describing the processes involved in the analysis step
  integer :: COMM_filter=0              !< MPI communicator processes in analysis step
  integer :: mype_filter=1              !< Process rank in COMM_da
  integer :: npes_filter=0              !< Number of processes in COMM_da

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks=1             !< Number of parallel model tasks
  INTEGER :: task_id=1                  !< Index of my model task (1,...,n_modeltasks)

  INTEGER :: MPIerr                     !< Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE) !< Status array for MPI

END MODULE mod_parallel_pdaf
