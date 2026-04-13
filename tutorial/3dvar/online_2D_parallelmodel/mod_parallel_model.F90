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
  INTEGER :: COMM_2dmodel  !< MPI communicator for model tasks
  INTEGER :: npes_2dmodel  !< Number of PEs in COMM_2dmodel
  INTEGER :: mype_2dmodel  !< PE rank in COMM_2dmodel
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
!! communicator Comm_2dmodel. It is also initialized
!! here together with its size (npes_2dmodel) and 
!! the rank of a PE (mype_2dmodel) within Comm_2dmodel.
!!
  SUBROUTINE init_parallel()

    IMPLICIT NONE

    INTEGER :: i
  
    CALL MPI_INIT(i);
    CALL MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for MPI_COMM_WORLD
    Comm_2dmodel = MPI_COMM_WORLD
    npes_2dmodel = npes_world
    mype_2dmodel = mype_world

#ifdef USE_PDAF
    ! Revise parallelization for ensemble assimilation
    CALL init_parallel_pdaf(1, Comm_2Dmodel, mype_2Dmodel, npes_2Dmodel)
#endif
   
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
