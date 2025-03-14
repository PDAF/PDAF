!$Id: init_parallel_pdaf.F90 1091 2011-08-16 13:55:36Z lnerger $
!BOP
!
! !ROUTINE: init_parallel_pdaf --- Initialize communicators for PDAF
!
! !INTERFACE:
SUBROUTINE init_parallel_pdaf(dim_ens, screen)

! !DESCRIPTION:
! Parallelization routine for a model with 
! attached PDAF. The subroutine is called in 
! the main program subsequently to the 
! initialization of MPI. It initializes
! MPI communicators for the model tasks, filter 
! tasks and the coupling between model and
! filter tasks. In addition some other variables 
! for the parallelization are initialized.
! The communicators and variables are handed
! over to PDAF in the call to 
! PDAF\_filter\_init.
!
! 3 Communicators are generated:\\
! - COMM\_filter: Communicator in which the
!   filter itself operates\\
! - COMM\_model: Communicators for parallel
!   model forecasts\\
! - COMM\_couple: Communicator for coupling
!   between models and filter\\
! Other variables that have to be initialized are:\\
! - filterpe - Logical: Does the PE execute the 
! filter?\\
! - my\_ensemble - Integer: The index of the PE's 
! model task\\
! - local\_npes\_model - Integer array holding 
! numbers of PEs per model task
!
! For COMM\_filter and COMM\_model also
! the size of the communicators (npes\_filter and 
! npes\_model) and the rank of each PE 
! (mype\_filter, mype\_model) are initialized. 
! These variables can be used in the model part 
! of the program, but are not handed over to PDAF.
!
! This variant is for a domain decomposed 
! model and the case that the filter runs on
! a single process that is separate from the 
! processes that run the ensemble integrations
! (practically a model task 0).
!
! NOTE: 
! This is a template that is expected to work 
! with many domain-decomposed models. However, 
! it might be necessary to adapt the routine 
! for a particular model. Inportant is that the
! communicator COMM_model equals the communicator 
! used in the model. If one plans to run a parallel 
! ensemble forecast (that is using multiple model
! tasks), COMM_model cannot be MPI_COMM_WORLD! Thus,
! if the model uses MPI_COMM_WORLD it has to be
! replaced by an alternative communicator named,
! e.g., COMM_model.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! 2014-04 - Lars Nerger - Variant for separate processes for model and filter
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE mod_parallel_model, &
       ONLY: mype_world, npes_world, mype_model, npes_model, &
       COMM_model, MPIerr, modelpe
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, filterpe, n_modeltasks, &
       local_npes_model, task_id, COMM_couple
  USE parser, &
       ONLY: parse

  IMPLICIT NONE    
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: dim_ens ! Ensemble size or number of EOFs (only SEEK)
  ! Often dim_ens=0 when calling this routine, because the real ensemble size
  ! is initialized later in the program. For dim_ens=0 no consistency check
  ! for ensemble size with number of model tasks is performed.
  INTEGER, INTENT(in)    :: screen ! Whether screen information is shown

! !CALLING SEQUENCE:
! Called by: main program
! Calls: MPI_Comm_size
! Calls: MPI_Comm_rank
! Calls: MPI_Comm_split
! Calls: MPI_Barrier
!EOP

  ! local variables
  INTEGER :: i, j               ! Counters
  INTEGER :: COMM_ensemble      ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens ! rank and size in COMM_ensemble
  INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
  INTEGER :: pe_index           ! Index of PE
  INTEGER :: my_color, color_couple ! Variables for communicator-splitting 
  LOGICAL :: iniflag            ! Flag whether MPI is initialized
  CHARACTER(len=32) :: handle   ! handle for command line parser
  INTEGER :: t_id


  ! *** Initialize MPI if not yet initialized ***
  CALL MPI_Initialized(iniflag, MPIerr)
  IF (.not.iniflag) THEN
     CALL MPI_Init(MPIerr)
  END IF

  ! *** Initialize PE information on COMM_world ***
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)


  ! *** Parse number of model tasks ***
  ! *** The module variable is N_MODELTASKS. Since it has to be equal
  ! *** to the ensemble size we parse dim_ens from the command line.
  handle = 'dim_ens'
  CALL parse(handle, n_modeltasks)


  ! *** Initialize communicators for ensemble evaluations ***
  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'Initialize communicators for assimilation with PDAF'


  ! *** Check consistency of number of parallel ensemble tasks ***
  consist1: IF (n_modeltasks > npes_world) THEN
     ! *** # parallel tasks is set larger than available PEs ***
     n_modeltasks = npes_world
     IF (mype_world == 0) WRITE (*, '(3x, a)') &
          '!!! Resetting number of parallel ensemble tasks to total number of PEs!'
  END IF consist1
  IF (dim_ens > 0) THEN
     ! Check consistency with ensemble size
     consist2: IF (n_modeltasks > dim_ens) THEN
        ! # parallel ensemble tasks is set larger than ensemble size
        n_modeltasks = dim_ens
        IF (mype_world == 0) WRITE (*, '(5x, a)') &
             '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
     END IF consist2
  END IF


  ! ***              COMM_ENSEMBLE                ***
  ! *** Generate communicator for ensemble runs   ***
  ! *** only used to generate model communicators ***
  COMM_ensemble = MPI_COMM_WORLD

  npes_ens = npes_world
  mype_ens = mype_world


  ! *** Store # PEs per ensemble                 ***
  ! *** used for info on PE 0 and for generation ***
  ! *** of model communicators on other Pes      ***
  ALLOCATE(local_npes_model(n_modeltasks+1))

  local_npes_model = FLOOR(REAL(npes_world-1) / REAL(n_modeltasks))
  local_npes_model(1) = 1
  DO i = 1, (npes_world - (n_modeltasks) * local_npes_model(2) - 1)
     local_npes_model(i+1) = local_npes_model(i+1) + 1
  END DO


  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
  pe_index = 0
  doens1: DO i = 1, n_modeltasks+1
     DO j = 1, local_npes_model(i)
        IF (mype_ens == pe_index) THEN
           task_id = i
           EXIT doens1
        END IF
        pe_index = pe_index + 1
     END DO
  END DO doens1
  IF (task_id==1) THEN
     t_id = MPI_UNDEFINED
  ELSE
     t_id = task_id
  END IF

  CALL MPI_Comm_split(COMM_ensemble, t_id, mype_ens, &
       COMM_model, MPIerr)

  ! Init flag FILTERPE (all PEs of model task 1)
  IF (task_id == 1) THEN
     filterpe = .TRUE.
     modelpe = .FALSE.
  ELSE
     filterpe = .FALSE.
     modelpe = .TRUE.
  END IF

  ! *** Reset task IDs ***
  task_id = task_id - 1

  IF (modelpe) THEN
     ! *** Re-initialize PE informations   ***
     ! *** according to model communicator ***
     CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
     CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

     IF (screen > 1) THEN
        WRITE (*,*) 'MODEL: mype(w)= ', mype_world, '; model task: ', task_id, &
             '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
     END IF
  END IF

  ! ***         COMM_FILTER                 ***
  ! *** Generate communicator for filter    ***
  ! *** For simplicity equal to COMM_couple ***

  IF (filterpe) THEN
     my_color = task_id
  ELSE
     my_color = MPI_UNDEFINED
  ENDIF

  CALL MPI_Comm_split(MPI_COMM_WORLD, my_color, mype_world, &
       COMM_filter, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  IF (filterpe) THEN
     CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
     CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

     IF (screen > 1) THEN
        WRITE (*,*) 'FILTER: mype(w)= ', mype_world, '; mype(f)= ', mype_filter, '; npes(f)= ', npes_filter
     END IF
  ENDIF

  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter PEs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  IF (modelpe .AND. mype_model==0) THEN
     color_couple = mype_model + 1
  ELSEIF (modelpe .AND. mype_model/=0) THEN
     color_couple = MPI_UNDEFINED
  ELSE
     color_couple = mype_filter + 1
  ENDIF

  CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple, mype_world, &
       COMM_couple, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  IF (COMM_couple /= MPI_COMM_NULL) THEN
     CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
     CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)
  ENDIF

  IF (screen > 0) THEN
     IF (mype_world == 0) THEN
        WRITE (*, '(/18x, a)') 'PE configuration:'
        WRITE (*, '(2x, a6, a9, a10, a14, 4x, a12, /2x, a5, a9, a7, a7, a7, a7, a8, /2x, a)') &
             'world', 'filter', 'model', 'couple', 'model/filter', &
             'rank', 'rank', 'task', 'rank', 'task', 'rank', 'M/F', &
             '----------------------------------------------------------'
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
     IF (filterpe) THEN
        IF (COMM_COUPLE/=MPI_COMM_NULL) THEN
           WRITE (*, '(2x, i4, 4x, i4, 4x, 3x, 4x, 3x, 4x, i3, 4x, i3, 5x, 3x, a1)') &
                mype_world, mype_filter, color_couple, mype_couple, 'F'
        ELSE
           WRITE (*, '(2x, i4, 4x, i4, 4x, 3x, 4x, 3x, 4x, 3x, 4x, 3x, 5x, 3x, a1)') &
                mype_world, mype_filter, 'F'
        END IF
     ENDIF
     IF (modelpe) THEN
        IF (COMM_COUPLE/=MPI_COMM_NULL) THEN
           WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, 3x, a1)') &
                mype_world, task_id, mype_model, color_couple, mype_couple, 'M'
        ELSE
           WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, 3x, 4x, 3x, 5x, 3x, a1)') &
                mype_world, task_id, mype_model, 'M'
        END IF
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

     IF (mype_world == 0) WRITE (*, '(/a)') ''

  END IF


! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

  ! If the names of the variables for COMM_model, npes_model, and 
  ! mype_model are different in the numerical model, the 
  ! model-internal variables should be initialized at this point.


END SUBROUTINE init_parallel_pdaf
