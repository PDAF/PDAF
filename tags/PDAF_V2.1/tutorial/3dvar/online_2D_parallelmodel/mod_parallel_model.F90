!$Id$
!> Module for model parallelization
!!
!! This module provides variables for the MPI parallelization
!! of the tutorial model to be shared between model-related routines. 
!!
!! In addition, methods to initialize and finalize MPI are provided.
!!
!! Revision history:
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_parallel_model

  USE mpi

  IMPLICIT NONE
  SAVE 

  ! Basic variables for model state integrations
  INTEGER :: COMM_model  !< MPI communicator for model tasks
  INTEGER :: npes_model  !< Number of PEs in COMM_model
  INTEGER :: mype_model  !< PE rank in COMM_model
  INTEGER :: npes_world  !< Number of PEs in MPI_COMM_WORLD
  INTEGER :: mype_world  !< PE rank in MPI_COMM_WORLD
  INTEGER :: MPIerr      !< Error flag for MPI
  
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

END MODULE mod_parallel_model
