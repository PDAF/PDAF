C     $Id$
C!BOP
C
C     !DESCRIPTION:
C     This modules provides variables for the MPI parallelization
C     to be shared between model-related routines. The are variables
C     that are used in the model, even without PDAF and additional
C     variables that are only used, if data assimialtion with PDAF
C     is performed.
C     In addition methods to initialize and finalize MPI are provided.
C     The initialization routine is only for the model itself, the
C     more complex initialization of communicators for xecution with
C     PDAF is peformed in init\_parallel\_pdaf.
C
C     !REVISION HISTORY:
C     2004-10 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !MODULE:
      MODULE mod_parallel_pdaf
      
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
C     Additional variables for use with PDAF
      INTEGER :: n_modeltasks = 1         ! Number of parallel model tasks
      INTEGER :: n_filterpes  = 1         ! Number of PEs for filter analysis
      INTEGER :: COMM_filter ! MPI communicator for filter PEs
      INTEGER :: mype_filter, npes_filter ! # PEs and PE rank in COMM_filter
      INTEGER :: COMM_couple ! MPI communicator for coupling filter and model
      LOGICAL :: modelpe     ! Whether we are on a PE in a COMM_model
      LOGICAL :: filterpe    ! Whether we are on a PE in a COMM_filter
      INTEGER :: task_id     ! Index of my model task (1,...,n_modeltasks)
      INTEGER :: MPIerr      ! Error flag for MPI
      INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
      INTEGER, ALLOCATABLE :: local_npes_model(:) ! # PEs per ensemble
C!EOP

      END MODULE mod_parallel_pdaf
