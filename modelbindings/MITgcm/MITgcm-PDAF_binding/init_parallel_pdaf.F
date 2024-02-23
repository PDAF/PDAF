C     $Id$
C
C     !ROUTINE: init_parallel_pdaf --- Initialize communicators for PDAF
C     !DESCRIPTION:
C     Parallelization routine for a model with
C     attached PDAF. The subroutine is called in
C     the main program subsequently to the
C     initialization of MPI. It initializes
C     MPI communicators for the model tasks, filter
C     tasks and the coupling between model and
C     filter tasks. In addition some other variables
C     for the parallelization are initialized.
C     The communicators and variables are handed
C     over to PDAF in the call to
C     PDAF\_init\_parallel.
C
C     3 Communicators are generated:\\
C     - COMM\_filter: Communicator in which the
C     filter itself operates\\
C     - COMM\_model: Communicators for parallel
C     model forecasts\\
C     - COMM\_couple: Communicator for coupling
C     between models and filter\\
C     Other variables that have to be initialized are:\\
C     - filterpe - Logical: Does the PE execute the
C     filter?\\
C     - my\_ensemble - Integer: The index of the
C     model task of PE\\
C     - local\_npes\_model - Integer array holding
C     numbers of PEs per model task
C
C     For COMM\_filter and COMM\_model also
C     the size of the communicators (npes\_filter and
C     npes\_model) and the rank of each PE
C     (mype\_filter, mype\_model) are initialized.
C     These variables can be used in the model part
C     of the program, but are not handed over to PDAF.
C
C     This variant is for a domain decomposed
C     model.
C
C     This is a template that is expected to work
C     with many domain-decomposed models. However,
C     it might be necessary to adapt the routine
C     for a particular model. Inportant is that the
C     communicator COMM_model equals the communicator
C     used in the model. If one plans to run the
C     ensemble forecast in parallel COMM_model cannot
C     be MPI_COMM_WORLDC Thus, if the model uses
C     MPI_COMM_WORLD it has to be replaced by an
C     alternative communicator named, e.g., COMM_model.
C
C     !REVISION HISTORY:
C     2004-11 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE init_parallel_pdaf(dim_ens, screen,
     &      MPI_COMM_MODEL, MPI_COMM_WORLD, mpi_task_id)
C     !USES:
      USE mod_parallel_model,
     &      ONLY: mype_model, npes_model, COMM_model,mype_world
      USE mod_parallel_pdaf,
     &      ONLY: mype_filter, npes_filter, COMM_filter,
     &      filterpe, n_modeltasks, local_npes_model, task_id,
     &      COMM_couple, MPIerr
      USE parser,
     &      ONLY: parse

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(inout) :: dim_ens ! Ensemble size
C     Often dim_ens=0 when calling this routine, because the real ensemble size
C     is initialized later in the program. For dim_ens=0 no consistency check
C     for ensemble size with number of model tasks is performed.
      INTEGER, INTENT(in)    :: screen ! Whether screen information is shown
      INTEGER, INTENT(out)   :: MPI_COMM_MODEL ! model communicator
      INTEGER, INTENT(in)    :: MPI_COMM_WORLD ! model communicator
      INTEGER, INTENT(out)   :: mpi_task_id

C     !CALLING SEQUENCE:
C     Called by: main program
C     Calls: MPI_Comm_size
C     Calls: MPI_Comm_rank
C     Calls: MPI_Comm_split
C     Calls: MPI_Barrier


C     local variables
      INTEGER :: i, j                      ! Counters
      INTEGER :: COMM_ensemble             ! Communicator of all PEs doing model tasks
      INTEGER :: mype_ens, npes_ens        ! rank and size in COMM_ensemble
      INTEGER :: mype_couple, npes_couple  ! Rank and size in COMM_couple
      INTEGER :: pe_index                  ! Index of PE
      INTEGER :: my_color, color_couple    ! Variables for communicator-splitting
      LOGICAL :: iniflag                   ! Flag whether MPI is initialized
      CHARACTER(len=32) :: handle          ! handle for command line parser
      INTEGER :: npes_world                ! Rank and size on MPI_COMM_WORLD
      INTEGER :: dim_ens_parsed            ! Parsed ensemble size

C     *** Initialize MPI if not yet initialized ***
      CALL MPI_Initialized(iniflag, MPIerr)
      IF (.not.iniflag) THEN
         CALL MPI_Init(MPIerr)
      END IF

C     *** Initialize PE information on COMM_world ***
      CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)

C     *** Parse number of model tasks ***
      handle = 'dim_ens'
      CALL parse(handle, dim_ens_parsed)
      n_modeltasks = dim_ens_parsed

C     *** Initialize communicators for ensemble evaluations ***
      IF (mype_world == 0)
     &     WRITE (*, '(/1x, a, a)') 'Initialize communicators for',
     &     ' assimilation with PDAF'


C     *** Check consistency of number of parallel ensemble tasks ***
      consist1: IF (n_modeltasks > npes_world) THEN
C     *** # parallel tasks is set larger than available PEs ***
         n_modeltasks = npes_world
         IF (mype_world == 0) WRITE (*, '(3x, a)')
     &        '!!! Resetting number of parallel ensemble tasks',
     &        ' to total number of PEs '
      END IF consist1
      IF (dim_ens > 0) THEN
C     Check consistency with ensemble size
        consist2: IF (n_modeltasks > dim_ens) THEN
C     # parallel ensemble tasks is set larger than ensemble size
        n_modeltasks = dim_ens
        IF (mype_world == 0) WRITE (*, '(5x, a)')
     &      '!!! Resetting number of parallel ensemble tasks',
     &      ' to number of ensemble states '
        END IF consist2
      END IF


C     ***              COMM_ENSEMBLE                ***
C     *** Generate communicator for ensemble runs   ***
C     *** only used to generate model communicators ***
      COMM_ensemble = MPI_COMM_WORLD

      npes_ens = npes_world
      mype_ens = mype_world


C     *** Store # PEs per ensemble                 ***
C     *** used for info on PE 0 and for generation ***
C     *** of model communicators on other Pes      ***
      ALLOCATE(local_npes_model(n_modeltasks))

      local_npes_model = FLOOR(REAL(npes_world) /
     &      REAL(n_modeltasks))
      DO i = 1, (npes_world - n_modeltasks *
     &      local_npes_model(1))
        local_npes_model(i) = local_npes_model(i) + 1
      END DO


C     ***              COMM_MODEL               ***
C     *** Generate communicators for model runs ***
C     *** (Split COMM_ENSEMBLE)                 ***
      pe_index = 0
      doens1: DO i = 1, n_modeltasks
         DO j = 1, local_npes_model(i)
            IF (mype_ens == pe_index) THEN
               task_id = i
               EXIT doens1
            END IF
            pe_index = pe_index + 1
         END DO
      END DO doens1

      CALL MPI_Comm_split(COMM_ensemble, task_id, mype_ens,
     &     COMM_model, MPIerr)

C     *** Re-initialize PE informations   ***
C     *** according to model communicator ***
      CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
      CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

      IF (screen > 1) THEN
         WRITE (*,*) 'MODEL: mype(w)= ', mype_world,
     &        ';model task: ', task_id, '; mype(m)= ',
     &        mype_model,'; npes(m)= ', npes_model
      END IF


C     Init flag FILTERPE (all PEs of model task 1)
      IF (task_id == 1) THEN
         filterpe = .TRUE.
      ELSE
         filterpe = .FALSE.
      END IF

C     ***         COMM_FILTER                 ***
C     *** Generate communicator for filter    ***
C     *** For simplicity equal to COMM_couple ***
      my_color = task_id

      CALL MPI_Comm_split(MPI_COMM_WORLD, my_color, mype_world,
     &     COMM_filter, MPIerr)

C     *** Initialize PE informations         ***
C     *** according to coupling communicator ***
      CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
      CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)


C     ***              COMM_COUPLE                 ***
C     *** Generate communicators for communication ***
C     *** between model and filter PEs             ***
C     *** (Split COMM_ENSEMBLE)                    ***

      color_couple = mype_filter + 1

      CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple,
     &     mype_world, COMM_couple, MPIerr)

C     *** Initialize PE informations         ***
C     *** according to coupling communicator ***
      CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
      CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

      IF (screen > 0) THEN
         IF (mype_world == 0) THEN
            WRITE (*, '(/12x, a)') 'Configuration of parallelization:'
            WRITE (*, '(2x,a6,a9,a10,a14,a13,/2x,a5,a9,5a7,/2x,53a)')
     &           'world', 'filter',
     &           'model', 'couple', 'filterPE', 'rank', 'rank',
     &           'task', 'rank', 'task', 'rank', 'T/F',
     &           ('-', i=1,53)
         END IF
         CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
         IF (task_id == 1) THEN
            WRITE (*, '(2x,i4,4x,i4,4x,i3,4x,i3,4x,i3,4x,i3,5x,l3)')
     &           mype_world, mype_filter, task_id,
     &           mype_model, color_couple, mype_couple, filterpe
         ENDIF
         IF (task_id > 1) THEN
            WRITE (*,'(2x,i4,12x,i3,4x,i3,4x,i3,4x,i3,5x,l3)')
     &            mype_world, task_id, mype_model,
     &           color_couple, mype_couple, filterpe
         END IF
         CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

         IF (mype_world == 0) WRITE (*, '(/a)') ' '

      END IF


C     ***********************************************************
C     *** Initialize model equivalents to COMM_model, task_id ***
C     *** to be passed back to MITgcm                         ***
C     ***********************************************************

C     Set model communicator of MITgcm to COMM_model
      MPI_COMM_MODEL = COMM_model

C     Set task ID of MITgcm to task_id initialize here
      mpi_task_id = task_id

      END SUBROUTINE init_parallel_pdaf
