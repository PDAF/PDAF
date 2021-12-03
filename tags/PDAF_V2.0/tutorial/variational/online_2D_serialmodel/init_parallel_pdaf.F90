!$Id$
!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with 
!! attached PDAF. The subroutine is called in 
!! the main program subsequently to the 
!! initialization of MPI. It initializes
!! MPI communicators for the model tasks, filter 
!! tasks and the coupling between model and
!! filter tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed
!! over to PDAF in the call to 
!! PDAF_filter_init.
!!
!! 3 Communicators are generated:
!! * _COMM_filter_: Communicator in which the
!!   filter itself operates
!! * _COMM_model_: Communicators for parallel
!!   model forecasts
!! * _COMM_couple_: Communicator for coupling
!!   between models and filter
!!
!! Other variables that have to be initialized are:
!! * _filterpe_ - Logical: Does the PE execute the 
!! filter?
!! * _my_ensemble_ - Integer: The index of the PE's 
!! model task
!! * _local_npes_model_ - Integer array holding 
!! numbers of PEs per model task
!!
!! For COMM_filter and COMM_model also
!! the size of the communicators (npes_filter and 
!! npes_model) and the rank of each PE 
!! (mype_filter, mype_model) are initialized. 
!! These variables can be used in the model part 
!! of the program, but are not handed over to PDAF.
!!
!! This variant is for a model without parallelization
!!
!! This is a template that is expected to work 
!! with many models without parallelization. However, 
!! it might be necessary to adapt the routine 
!! for a particular model. Inportant is that the
!! communicator COMM_model equals the communicator 
!! used in the model. If one plans to run the
!! ensemble forecast in parallel COMM_model cannot 
!! be MPI_COMM_WORLD! Thus, if the model uses 
!! MPI_COMM_WORLD it has to be replaced by an 
!! alternative communicator named, e.g., COMM_model.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_parallel_pdaf(dim_ens, screen)

  USE mpi                         ! MPI
  USE mod_parallel_pdaf, &        ! PDAF parallelization variables
       ONLY: mype_world, npes_world, mype_model, npes_model, &
       COMM_model, mype_filter, npes_filter, COMM_filter, filterpe, &
       n_modeltasks, local_npes_model, task_id, COMM_couple, MPIerr
  USE parser, &                   ! Command line parser
       ONLY: parse

  IMPLICIT NONE    
  
! *** Arguments ***
  INTEGER, INTENT(inout) :: dim_ens !< Ensemble size or number of EOFs (only SEEK)
  !< Often dim_ens=0 when calling this routine, because the real ensemble size
  !< is initialized later in the program. For dim_ens=0 no consistency check
  !< for the ensemble size with the number of model tasks is performed.
  INTEGER, INTENT(in)    :: screen !< Whether screen information is shown

! *** local variables ***
  INTEGER :: i, j               ! Counters
  INTEGER :: COMM_ensemble      ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens ! rank and size in COMM_ensemble
  INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
  INTEGER :: pe_index           ! Index of PE
  INTEGER :: my_color, color_couple ! Variables for communicator-splitting 
  LOGICAL :: iniflag            ! Flag whether MPI is initialized
  CHARACTER(len=32) :: handle   ! handle for command line parser


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
  consist1: IF (n_modeltasks /= npes_world) THEN
     ! *** # parallel tasks is not identical to available PEs     ***
     ! *** This needs to hold for a model without parallelization ***
     IF (mype_world == 0) WRITE (*, '(3x, a)') &
          '!!! ERROR: dim_ens must equal number of processes!'
     CALL MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)
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
  ALLOCATE(local_npes_model(n_modeltasks))

  local_npes_model = FLOOR(REAL(npes_world) / REAL(n_modeltasks))
  DO i = 1, (npes_world - n_modeltasks * local_npes_model(1))
     local_npes_model(i) = local_npes_model(i) + 1
  END DO
  

  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
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


  CALL MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
       COMM_model, MPIerr)
  
  ! *** Re-initialize PE informations   ***
  ! *** according to model communicator ***
  CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
  CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

  if (screen > 1) then
    write (*,*) 'MODEL: mype(w)= ', mype_world, '; model task: ', task_id, &
         '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
  end if


  ! Init flag FILTERPE (all PEs of model task 1)
  IF (task_id == 1) THEN
     filterpe = .TRUE.
  ELSE
     filterpe = .FALSE.
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
  ENDIF


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter PEs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = mype_model + 1

  CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple, mype_world, &
       COMM_couple, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

  IF (screen > 0) THEN
     IF (mype_world == 0) THEN
        WRITE (*, '(/18x, a)') 'PE configuration:'
        WRITE (*, '(2x, a6, a9, a10, a14, a13, /2x, a5, a9, a7, a7, a7, a7, a7, /2x, a)') &
             'world', 'filter', 'model', 'couple', 'filterPE', &
             'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             '----------------------------------------------------------'
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
     IF (task_id == 1) THEN
        WRITE (*, '(2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
             mype_world, mype_filter, task_id, mype_model, color_couple, &
             mype_couple, filterpe
     ENDIF
     IF (task_id > 1) THEN
        WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
         mype_world, task_id, mype_model, color_couple, mype_couple, filterpe
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
