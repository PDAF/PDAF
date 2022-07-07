C     $Id$
C!BOP
C
C     !DESCRIPTION:
C     This modules provides variables for the MPI parallelization
C     of the tutorial model to be shared between model-related routines.
C
C     In addition, methods to initialize and finalize MPI are provided.
C
C     !REVISION HISTORY:
C     2004-10 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !MODULE:
      MODULE mod_parallel_model

C     !USES:

#ifdef ALLOW_USE_MPI
      USE mpi
#else
These lines are here on purpose to STOP the compilation,
because online-PDAF needs MPI
#endif

      IMPLICIT NONE
      SAVE

C     !PUBLIC DATA MEMBERS:
C     Basic variables for model state integrations
      INTEGER :: COMM_model  ! MPI communicator for model tasks
      INTEGER :: mype_model  ! Number of PEs in COMM_model
      INTEGER :: npes_model  ! PE rank in COMM_model
      INTEGER :: mype_world  ! Number of PEs in MPI_COMM_WORLD
      INTEGER :: npes_world  ! PE rank in MPI_COMM_WORLD
      INTEGER :: MPIerr      ! Error flag for MPI

C!EOP

      CONTAINS
C--------------------------------------------------------------
C!BOP
C
C     !ROUTINE: init_parallel - Initialize MPI
C
C     !INTERFACE:
      SUBROUTINE init_parallel()

C     !DESCRIPTION:
C     Routine to initialize MPI, the number of PEs
C     (npes\_world) and the rank of a PE (mype\_world).
C     The model is executed within the scope of the
C     communicator Comm_model. It is also initialized
C     here together with its size (npes\_model) and
C     the rank of a PE (mype\_model) within Comm_model.
C!EOP

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
C-------------------------------------------------------------------------------
C!BOP
C
C     !ROUTINE: finalize_parallel - Finalize MPI
C
C     !INTERFACE:
      SUBROUTINE finalize_parallel()

C     !DESCRIPTION:
C     Routine to finalize MPI
C!EOP

      IMPLICIT NONE

      CALL  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
      CALL  MPI_Finalize(MPIerr)

      END SUBROUTINE finalize_parallel
C-------------------------------------------------------------------------------
C!BOP
C
C     !ROUTINE: abort_parallel - Abort MPI
C
C     !INTERFACE:
      SUBROUTINE abort_parallel()

C     !DESCRIPTION:
C     Routine to abort MPI program
C!EOP

      IMPLICIT NONE

      CALL  MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

      END SUBROUTINE abort_parallel

      END MODULE mod_parallel_model
