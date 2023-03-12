!$Id: mod_parallel_model.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_parallel_model

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! of the tutorial model to be shared between model-related routines. 
!
! In addition, methods to initialize and finalize MPI are provided.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  ! Basic variables for model state integrations
  INTEGER :: COMM_model  ! MPI communicator for model tasks
  INTEGER :: mype_model  ! Number of PEs in COMM_model
  INTEGER :: npes_model  ! PE rank in COMM_model
  INTEGER :: mype_world  ! Number of PEs in MPI_COMM_WORLD
  INTEGER :: npes_world  ! PE rank in MPI_COMM_WORLD
  INTEGER :: MPIerr      ! Error flag for MPI
!EOP
  
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_parallel - Initialize MPI
!
! !INTERFACE:
  SUBROUTINE init_parallel()

! !DESCRIPTION:
! Routine to initialize MPI, the number of PEs
! (npes\_world) and the rank of a PE (mype\_world).
! The model is executed within the scope of the
! communicator Comm_model. It is also initialized
! here together with its size (npes\_model) and 
! the rank of a PE (mype\_model) within Comm_model.
!EOP

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
!BOP
!
! !ROUTINE: finalize_parallel - Finalize MPI
!
! !INTERFACE:
  SUBROUTINE finalize_parallel()

! !DESCRIPTION:
! Routine to finalize MPI
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    CALL  MPI_Finalize(MPIerr)

  END SUBROUTINE finalize_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: abort_parallel - Abort MPI
!
! !INTERFACE:
  SUBROUTINE abort_parallel()

! !DESCRIPTION:
! Routine to abort MPI program
!EOP

    IMPLICIT NONE
    
    CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_model
