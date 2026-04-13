! Copyright (c) 2004-2026 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!> Initialize PDAF
!!
!! Initialization of PDAF. Performed are:\\
!!   * Initialization of filter independent parameters\\
!!   * Call to filter-specific routine for parameter initialization\\
!!   * Initialization of PDAF-internal parallelization\\
!!   * Call to filter-specific routine for allocation of arrays\\
!!   * Call to user-routine for ensemble/mode initialization.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-08 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
MODULE PDAF3init

CONTAINS

  SUBROUTINE PDAF3_init(filtertype, subtype, stepnull, param_int, dim_pint, &
       param_real, dim_preal, U_init_ens, in_screen, outflag)

    USE mpi
    USE PDAF_cb_procedures
    USE PDAF_timer, &
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount_ini
    USE PDAF_mod_core, &
         ONLY: dim_ens, dim_eof, dim_p, flag, &
         screen, step, step_obs, type_filter, filterstr, &
         subtype_filter, ensemblefilter, state, Ainv, ens, &
         debug
    USE PDAF_mod_parallel, &
         ONLY: mype, filterpe, PDAF_MPI_init, PDAF_init_parallel, COMM_pdaf, &
         isset_comm_pdaf, COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks
    USE PDAF_info, &
         ONLY: PDAF_print_version
    USE PDAF_DA, ONLY: &
         PDAF_print_filter_types
    USE PDAF_utils_filters, &
         ONLY: PDAF_init_filters, PDAF_alloc_filters, PDAF_options_filters

    IMPLICIT NONE

! *** Arguments ***
    ! For valid and default values see PDAF_mod_core.F90
    INTEGER, INTENT(in) :: filtertype     !< Type of filter
    INTEGER, INTENT(in) :: subtype        !< Sub-type of filter
    INTEGER, INTENT(in) :: stepnull       !< Initial time step of assimilation
    INTEGER, INTENT(in) :: dim_pint       !< Number of integer parameters
    INTEGER, INTENT(inout) :: param_int(dim_pint) !< Integer parameter array
    INTEGER, INTENT(in) :: dim_preal      !< Number of real parameter 
    REAL, INTENT(inout) :: param_real(dim_preal) !< Real parameter array
    INTEGER, INTENT(in) :: in_screen      !< Control screen output:
                                          !< (0) none, (1) some, default, (2) extensive
    INTEGER, INTENT(out):: outflag        !< Status flag, 0: no error, error codes:
                                   !< -1: Call with subtype=-1 for info display
                                   !<  1: No valid filter type
                                   !<  2: No valid sub type
                                   !<  3: Invalid dim_pint
                                   !<  4: Invalid dim_preal
                                   !<  5: Invalid state dimension
                                   !<  6: Invalid ensemble size
                                   !<  7: Invalid value for forgetting factor
                                   !<  8: Invalid other integer parameter value
                                   !<  9: Invalid other real parameter value
                                   !< 10: MPI information not initialized
                                   !< 20: error in allocation of array at PDAF init

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
    PROCEDURE(init_ens_cb) :: U_init_ens  !< User-supplied routine for ensemble initialization

! *** local variables ***
    INTEGER :: i                     ! Counter
    LOGICAL :: fixedbasis            ! Does the filter run with fixed error-space basis (EnOI mode)?


! ********************************************
! *** INITIALIZE VARIABLES FOR ALL FILTERS ***
! ********************************************

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_init -- START'

    ! Initialize MPI if not done before calling PDAF3_init
    CALL PDAF_MPI_init()

    ! set number of timers
    CALL PDAF_timeit(66, 'ini')

    ! Initialize memory counters
    CALL PDAF_memcount_ini(6)

    ! Call timer
    CALL PDAF_timeit(1, 'new')

    ! Set PDAF communicator if not set externally
    IF (.NOT. isset_comm_pdaf) THEN
       COMM_pdaf = MPI_COMM_WORLD

       IF (debug>0) &
            WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'Use MPI_COMM_WORLD for COMM_PDAF'
    ELSE
       IF (debug>0) &
            WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'Use user-defined communicator for COMM_PDAF'
    END IF

    ! Print version information
    CALL PDAF_print_version()

    info_assim: IF (subtype < 0) THEN

       ! ***********************************************************************
       ! *** For negative subtype only display information on filter options ***
       ! ***********************************************************************

       CALL PDAF_print_filter_types(1)
       CALL PDAF_options_filters(filtertype)

       subtype_filter = -1

       ! Set status flag
       flag = -1

    ELSE info_assim

       ! *** Check size of parameter arrays
       IF (dim_pint < 2) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(3): Invalid size of array of integer parameters!'
          flag = 3
       END IF
       IF (dim_preal < 1) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(4): Invalid size of array of real parameters!'
          flag = 4
       END IF

       ! *** Initialize variables for all PEs
       type_filter    = filtertype   ! Set filter type
       subtype_filter = subtype      ! Set sub-type of filter
       step_obs       = stepnull     ! Set initial time step
       step           = step_obs + 1 ! stepping index
       screen         = in_screen    ! Control verbosity

       dim_p          = param_int(1) ! PE-local state dimension
       IF (param_int(1) < 1) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(5): Invalid state dimension!'
          flag = 5
       END IF

       ! Ensemble size
       dim_ens = param_int(2)
       IF (param_int(2) < 1) THEN
          WRITE (*,'(/5x,a/)') 'PDAF-ERROR(6): Invalid ensemble size!'
          flag = 6
       END IF

       IF (debug>0 .AND. flag==0) THEN
          WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'param_int of size', dim_pint, &
               'values:', param_int(1:dim_pint)
          WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'param_real of size', dim_preal, &
               'values:', param_real(1:dim_preal)
          WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, &
               '  Note: If REAL values appear incorrect, please check if you provide them with the correct precision'
       END IF


! ********************************************
! *** Initialize filter-specific variables ***
! ********************************************

       IF (flag == 0) THEN
          CALL PDAF_init_filters(type_filter, subtype_filter, param_int, dim_pint, param_real, &
               dim_preal, filterstr, ensemblefilter, fixedbasis, screen, flag)
       END IF

  
! **********************************
! *** Initialize parallelization ***
! **********************************

       IF (flag == 0) THEN
          CALL PDAF_init_parallel(dim_ens, ensemblefilter, fixedbasis, &
               COMM_model, COMM_filter, COMM_couple, &
               n_modeltasks, task_id, screen, flag)
       END IF


! ********************************************
! *** Filter-specific allocation of arrays ***
! *** and screen output                    ***
! ********************************************

       IF (flag == 0) THEN
          CALL PDAF_alloc_filters(filterstr, subtype_filter, flag)
       END IF


! **********************************
! *** Initialize ensemble matrix ***
! **********************************

       filter_pe3: IF (filterpe .AND. flag == 0) THEN

          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(/a)') 'PDAF: Call ensemble initialization'

          CALL PDAF_timeit(39, 'new')

          typef: IF (ensemblefilter) THEN
! *** Initialize ensemble of ensemble-based filter      ***
! *** EnKF/SEIK/LSEIK/ETKF/LETKF                        ***
             CALL U_init_ens(type_filter, dim_p, dim_ens, state, Ainv, &
                  ens, flag)

             IF (debug>0) THEN
                DO i = 1, dim_ens
                   WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'ensemble member', i, &
                        ' values (1:min(dim_p,6)):', ens(1:MIN(dim_p,6),i)
                END DO
             END IF
          ELSE
! *** Mode-based filter (SEEK)                          ***
! *** Initialize rank reduced covariance matrix         ***
! *** factors Ainv and ens and estimated initial state ***
             CALL U_init_ens(type_filter, dim_p, dim_eof, state, Ainv, &
                  ens, flag)

             IF (debug>0) THEN
                DO i = 1, dim_ens
                   WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'covar mode', i, &
                        ' values (1:min(dim_p,6)):', ens(1:MIN(dim_p,6),i)
                END DO
                WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'mode weights (1:min(dim_eof,10)):', &
                     Ainv(1:MIN(dim_eof, 10),1:MIN(dim_eof, 10))
             END IF
          END IF typef

          CALL PDAF_timeit(39, 'old')

       END IF filter_pe3

    END IF info_assim


! ********************
! *** FINISHING UP ***
! ********************

    ! Store internal status flag
    outflag = flag

    CALL PDAF_timeit(1, 'old')

    IF (mype == 0 .AND. filterpe .AND. screen > 0) &
         WRITE (*, '(/a)') 'PDAF: Initialization completed'
    IF (mype == 0 .AND. filterpe .AND. screen > 1) &
         WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
         'PDAF', '--- duration of PDAF initialization:', PDAF_time_temp(1), 's'

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_init -- END'

  END SUBROUTINE PDAF3_init


!-------------------------------------------------------------------------------
!> Interface to control ensemble integration
!!
!! Interface routine called from the model before the 
!! forecast of each ensemble state to transfer data
!! from PDAF to the model.  For the parallelization 
!! this involves transfer from filter PEs to model PEs.
!!
!! At the beginning of a forecast phase sub-ensembles
!! are distributed to the model tasks. During the 
!! forecast phase each state vector of a sub-ensemble
!! is transferred to the model fields by U\_dist\_state.
!!
!! !  This is a core routine of PDAF and
!! should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-05 - Lars Nerger - Initial code as alias of PDAF_get_state
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_init_forecast(U_next_observation, U_distribute_state, &
       U_prepoststep, outflag)

    USE PDAF_cb_procedures
    USE PDAFget_state, ONLY: PDAF_get_state

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag !< Status flag

! Loca variables
    INTEGER :: steps   !< Flag and number of time steps
    REAL    :: time    !< current model time
    INTEGER :: doexit  !< Whether to exit from forecasts


! *** Argument procedures ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    PROCEDURE(distribute_cb) :: U_distribute_state   !< Routine to distribute a state vector
    PROCEDURE(prepost_cb) :: U_prepoststep           !< User supplied pre/poststep routine
    PROCEDURE(next_obs_cb) :: U_next_observation     !< Provide information on next forecast


    CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
         U_prepoststep, outflag)

  END SUBROUTINE PDAF3_init_forecast

!-------------------------------------------------------------------------------
!> Interface to set parallelization information
!!
!! Interface routine called from init_parallel_pdaf
!! providing MPI parallelization information to PDAF.
!!
!! !  This is a core routine of PDAF and
!! should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_set_parallel(in_COMM_pdaf, in_COMM_model, in_COMM_filter, in_COMM_couple, &
       in_task_id, in_n_modeltasks, in_filterpe, flag)

    USE mpi
    USE PDAF_mod_parallel, &
         ONLY: COMM_pdaf, COMM_model, COMM_filter, COMM_couple, &
         task_id, n_modeltasks, filterpe, isset_comm_pdaf, isset_parallel

    IMPLICIT NONE    

! *** Arguments ***
    INTEGER, INTENT(in) :: in_COMM_pdaf      !< MPI communicator for all PEs involved in PDAF
    INTEGER, INTENT(in) :: in_COMM_model     !< Model communicator
    INTEGER, INTENT(in) :: in_COMM_filter    !< Filter communicator
    INTEGER, INTENT(in) :: in_COMM_couple    !< Coupling communicator
    INTEGER, INTENT(in) :: in_task_id        !< Task ID of current PE
    INTEGER, INTENT(in) :: in_n_modeltasks   !< Number of model tasks
    LOGICAL, INTENT(in) :: in_filterpe       !< Is my PE a filter-PE?
    INTEGER, INTENT(inout):: flag            !< Status flag

! *** Initialize internal communicators ***
    COMM_pdaf = in_COMM_pdaf
    COMM_model = in_COMM_model
    COMM_filter = in_COMM_filter
    COMM_couple = in_COMM_couple

! *** Initialize internal parameters ***
    task_id = in_task_id
    n_modeltasks = in_n_modeltasks
    filterpe = in_filterpe
    isset_parallel = .TRUE.
    IF (COMM_PDAF /= MPI_COMM_WORLD) isset_comm_pdaf = .TRUE.

! *** Set status flag ***
    flag = 0

  END SUBROUTINE PDAF3_set_parallel



!-------------------------------------------------------------------------------
!>  Initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF. The subroutine is
!! called in the main program subsequently to the initialization of MPI. It
!! initializes MPI communicators for the model tasks, assimilation task and the
!! coupling between model and assimilation tasks. In addition some other variables 
!! for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF in the call to 
!! PDAF_set_parallel toward the end of this routine.
!!
!! 3 Communicators are generated:
!! * _COMM_assim_: Communicator in which the assimilation analysis is computed
!! * _COMM_model_: Communicators for parallel model forecasts
!! * _COMM_couple_: Communicator for coupling between model and assi. processes
!!
!! In addition there is the main communicator
!! * _COMM_ensemble_: The main communicator in which PDAF operates
!! COMM_ensemble is set to the communicator in which all model integration 
!! are computed. Typically, this is MPI_COMM_WORLD, but it can be defined
!! differently if the model only operators on a subset to MPI_COMM_WORLD.
!! This happens, e.g. if some processes are separated to operate an
!! I/O server or a model coupler for coupled model systems.
!!
!! Other variables that have to be initialized are:
!! * _assimpe_ - Logical: Does the Process execute the analysis step?
!! * _task_id_ - Integer: Index identifying the model task
!! * _my_ensemble_ - Integer: The index of the Process's model task
!! * _local_npes_model_ - Integer array holding numbers of Processs per model task
!!
!! For COMM_assim and COMM_model also the size of the communicators
!! (npes_assim and npes_model) and the rank of each process  (mype_assim,
!! mype_model) are initialized.
!!
!! __Revision history:__
!! * 2026-03 - Lars Nerger - Initial code moving functionality from user code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAF3_init_parallel(screen, type_parallel, online_coupling, dim_ens, n_modeltasks, &
       COMM_model, mype_model, npes_model, COMM_assim, mype_assim, npes_assim, &
       task_id)

    USE mpi
    USE PDAF_mod_parallel, &
         ONLY: COMM_pdaf, COMM_model_mod=>COMM_model, COMM_filter_mod=>COMM_filter, &
         COMM_couple, task_id_mod=>task_id, n_modeltasks_mod=>n_modeltasks, filterpe, &
         modelpe, isset_comm_pdaf, isset_parallel, mpi_init_by_pdaf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: screen            !< Whether screen information is shown

    ! Model variables for parallelization
    INTEGER, INTENT(in) :: type_parallel        !< Type of parallelization
                                                !< 0: common setup using task 1 for assimilation
                                                !< 1: setup using separate task 0 for assimilation
    INTEGER, INTENT(in) :: online_coupling      !< 1: online DA coupling, 0: offline DA coupling
    INTEGER, INTENT(in) :: dim_ens              !< Ensemble size / number of model tasks
    INTEGER, INTENT(inout) :: n_modeltasks      !< Number of model tasks
    INTEGER, INTENT(inout) :: COMM_model        !< Model MPI communicator for model tasks
    INTEGER, INTENT(out) :: npes_model          !< Number of Processs in COMM_model
    INTEGER, INTENT(out) :: mype_model          !< Process rank in COMM_model
    INTEGER, INTENT(out) :: COMM_assim          !< MPI communicator for assimilation processes 
    INTEGER, INTENT(out) :: npes_assim          !< Number of processes in COMM_assim
    INTEGER, INTENT(out) :: mype_assim          !< Process rank in COMM_assim
    INTEGER, INTENT(out) :: task_id             !< Index of my model task (1,...,n_modeltasks)

! *** Local variables ***
    INTEGER :: i, j                             ! Counters
    INTEGER :: pe_index                         ! Index of Process
    INTEGER :: MPIerr                           ! Error flag for MPI
    INTEGER :: my_color, color_couple           ! Variables for communicator-splitting 
    INTEGER, ALLOCATABLE :: local_npes_model(:) ! Number of processes per ensemble
    INTEGER :: mype_couple                      ! Rank in COMM_couple
    INTEGER :: npes_couple                      ! Size in COMM_couple
    INTEGER :: mype_pdaf                        ! Rank in COMM_pdaf
    INTEGER :: npes_pdaf                        ! Size of COMM_pdaf
    INTEGER :: dummy                            ! Dummy variable to avoid compiler warning
    LOGICAL :: iniflag                          ! Flag whether MPI is initialized
    INTEGER :: t_id                             ! Variable for storing task id for communicator splitting


! ************************************************
! *** Initialize communicators for ensemble DA ***
! ************************************************

    ! Dummy init to avoid compiler warning
    dummy = type_parallel

    ! *** Fix number or model tasks for offline DA ***
    IF (online_coupling==0) n_modeltasks = 1

    partype: IF (type_parallel /= 1) THEN
       
! **************************************************
! *** Common setup using task 1 for assimilation ***
! **************************************************

       ! All processes are model PEs
       modelpe = .TRUE.

! *** Initialize MPI if not yet initialized ***

       CALL MPI_Initialized(iniflag, MPIerr)
       IF (.NOT.iniflag) THEN
          CALL MPI_Init(MPIerr)

          mpi_init_by_pdaf = .TRUE.

          COMM_model = MPI_COMM_WORLD

          CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)
          IF (mype_model==0 .AND. screen>0) &
               WRITE (*, '(/a, 2x, a)') 'PDAF', 'MPI-initialization by PDAF'
       END IF

! ***                   COMM_PDAF                     ***
! *** This is the communicator in which PDAF operates ***

       COMM_pdaf = COMM_model

       ! *** Get rank and size of COMM_pdaf

       CALL MPI_Comm_Size(COMM_pdaf, npes_pdaf, MPIerr)
       CALL MPI_Comm_Rank(COMM_pdaf, mype_pdaf, MPIerr)


       ! Initial screen output
       IF (mype_pdaf == 0 .AND. screen>0) &
            WRITE (*, '(/a, 2x, a)') 'PDAF', '*** Initialize MPI communicators for assimilation with PDAF ***'

       ! *** Check consistency of number of parallel ensemble tasks ***
       IF (online_coupling==1) THEN
          consist1: IF (n_modeltasks > npes_pdaf) THEN
             ! *** # parallel tasks is set larger than available Processs ***
             n_modeltasks = npes_pdaf
             IF (mype_pdaf == 0) WRITE (*, '(a, 3x, a)') &
                  'PDAF', '!!! Resetting number of parallel ensemble tasks to total number of Processs!'
          END IF consist1
          IF (dim_ens > 0) THEN
             ! Check consistency with ensemble size
             consist2: IF (n_modeltasks > dim_ens) THEN
                ! # parallel ensemble tasks is set larger than ensemble size
                n_modeltasks = dim_ens
                IF (mype_pdaf == 0) WRITE (*, '(a, 5x, a)') &
                     'PDAF', '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
             END IF consist2
          END IF
       END IF


! *** Store # of processes per model task           ***
! *** used for info on Process 0 and for generation ***
! *** of model communicators on other Pes           ***

       ALLOCATE(local_npes_model(n_modeltasks))

       local_npes_model = FLOOR(REAL(npes_pdaf) / REAL(n_modeltasks))
       DO i = 1, (npes_pdaf - n_modeltasks * local_npes_model(1))
          local_npes_model(i) = local_npes_model(i) + 1
       END DO


! ***              COMM_MODEL               ***
! *** Generate communicators for model runs ***
! *** (Split COMM_PDAF)                     ***

       pe_index = 0
       doens1: DO i = 1, n_modeltasks
          DO j = 1, local_npes_model(i)
             IF (mype_pdaf == pe_index) THEN
                task_id = i
                EXIT doens1
             END IF
             pe_index = pe_index + 1
          END DO
       END DO doens1

       CALL MPI_Comm_split(COMM_pdaf, task_id, mype_pdaf, &
            COMM_model, MPIerr)
  
       ! *** Re-initialize Process information for COMM_model

       CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
       CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

       IF (screen > 1) THEN
          WRITE (*,*) 'PDAF: mype(w)= ', mype_pdaf, '; model task: ', task_id, &
               '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
       END IF


! *** Init flag for assim processes    ***
! *** (all processes of model task 1)  ***

       IF (task_id == 1) THEN
          filterpe = .TRUE.
       ELSE
          filterpe = .FALSE.
       END IF


! ***         COMM_ASSIM                  ***
! *** Generate communicator for analysis  ***

       IF (filterpe) THEN
          my_color = task_id
       ELSE
          my_color = MPI_UNDEFINED
       ENDIF

       CALL MPI_Comm_split(COMM_pdaf, my_color, mype_pdaf, &
            COMM_assim, MPIerr)

       ! *** Initialize Process information for COMM_assim

       IF (filterpe) THEN
          CALL MPI_Comm_Size(COMM_assim, npes_assim, MPIerr)
          CALL MPI_Comm_Rank(COMM_assim, mype_assim, MPIerr)
       ENDIF


! ***              COMM_COUPLE                 ***
! *** Generate communicators for communication ***
! *** between model and assim processes        ***
! *** (Split COMM_pdaf)                        ***

       color_couple = mype_model + 1

       CALL MPI_Comm_split(COMM_pdaf, color_couple, mype_pdaf, &
            COMM_couple, MPIerr)

       ! *** Initialize Process information for COMM_couple

       CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
       CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)


! *** Display process configuration ***

       IF (screen > 0) THEN
          IF (mype_pdaf == 0) THEN
             WRITE (*, '(a13, 2x, a)') 'PDAF    Pconf', 'Process configuration:'
             WRITE (*, '(a13, 2x, a6, a9, a12, a17, a15, /a13, 2x, a5, a9, a8, a9, a8, a9, a9, /a13, 2x, a)') &
                  'PDAF    Pconf', 'world', 'assim', 'model', 'couple', 'assimPE', &
                  'PDAF    Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
                  'PDAF    Pconf', '------------------------------------------------------------'
          END IF
          CALL MPI_Barrier(COMM_pdaf, MPIerr)
          IF (task_id == 1) THEN
             WRITE (*, '(a, 2x, i5, 4x, i5, 4x, i4, 4x, i5, 4x, i4, 4x, i5, 5x, l3)') &
                  'PDAF    Pconf', mype_pdaf, mype_assim, task_id, mype_model, color_couple, &
                  mype_couple, filterpe
          ENDIF
          IF (task_id > 1) THEN
             WRITE (*,'(a, 2x, i5, 13x, i4, 4x, i5, 4x, i4, 4x, i5, 5x, l3)') &
                  'PDAF    Pconf', mype_pdaf, task_id, mype_model, color_couple, mype_couple, filterpe
          END IF
          CALL MPI_Barrier(COMM_pdaf, MPIerr)

          IF (mype_pdaf == 0) WRITE (*, '(/a)') ''
       END IF


! *** Store parallelization information for internal use ***

       COMM_model_mod = COMM_model
       COMM_filter_mod = COMM_assim
       task_id_mod = task_id
       n_modeltasks_mod = n_modeltasks

       ! Set flags
       isset_parallel = .TRUE.
       IF (COMM_PDAF /= MPI_COMM_WORLD) isset_comm_pdaf = .TRUE.



    ELSE partype

! **********************************************************
! *** Setup using a separate task 0 for the assimilation ***
! **********************************************************

! *** Initialize MPI if not yet initialized ***

       CALL MPI_Initialized(iniflag, MPIerr)
       IF (.NOT.iniflag) THEN
          CALL MPI_Init(MPIerr)

          mpi_init_by_pdaf = .TRUE.

          COMM_model = MPI_COMM_WORLD

          CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)
          IF (mype_model==0 .AND. screen>0) &
               WRITE (*, '(/a, 2x, a)') 'PDAF', 'MPI-initialization by PDAF'
       END IF

! ***                   COMM_PDAF                     ***
! *** This is the communicator in which PDAF operates ***

       COMM_pdaf = COMM_model

       ! *** Get rank and size of COMM_pdaf

       CALL MPI_Comm_Size(COMM_pdaf, npes_pdaf, MPIerr)
       CALL MPI_Comm_Rank(COMM_pdaf, mype_pdaf, MPIerr)


       ! Initial screen output
       IF (mype_pdaf == 0 .AND. screen>0) THEN
          WRITE (*, '(/a, 2x, a)') 'PDAF', '*** Initialize MPI communicators for assimilation with PDAF ***'
          WRITE (*, '(a, 2x, a)') 'PDAF', '***   configure separate tasks for models and assimilation  ***'
       END IF

       ! *** Check consistency of number of parallel ensemble tasks ***
       IF (online_coupling==1) THEN
          consist1a: IF (n_modeltasks > npes_pdaf) THEN
             ! *** # parallel tasks is set larger than available Processs ***
             n_modeltasks = npes_pdaf
             IF (mype_pdaf == 0) WRITE (*, '(a, 3x, a)') &
                  'PDAF', '!!! Resetting number of parallel ensemble tasks to total number of Processs!'
          END IF consist1a
          IF (dim_ens > 0) THEN
             ! Check consistency with ensemble size
             consist2a: IF (n_modeltasks > dim_ens) THEN
                ! # parallel ensemble tasks is set larger than ensemble size
                n_modeltasks = dim_ens
                IF (mype_pdaf == 0) WRITE (*, '(a, 5x, a)') &
                     'PDAF', '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
             END IF consist2a
          END IF
       END IF


       ! *** Store # PEs per ensemble                 ***
       ! *** used for info on PE 0 and for generation ***
       ! *** of model communicators on other Pes      ***
       ALLOCATE(local_npes_model(n_modeltasks+1))

       ! This step takes into account that we have a 'task 0' 
       ! that computes the analysis step in PDAF.
       local_npes_model = FLOOR(REAL(npes_pdaf) / REAL(n_modeltasks+1))
       DO i = 1, (npes_pdaf - (n_modeltasks+1) * local_npes_model(1))
          local_npes_model(i) = local_npes_model(i) + 1
       END DO



! ***              COMM_MODEL               ***
! *** Generate communicators for model runs ***
! *** (Split COMM_PDAF)                     ***

       pe_index = 0
       doens1a: DO i = 1, n_modeltasks+1
          DO j = 1, local_npes_model(i)
             IF (mype_pdaf == pe_index) THEN
                task_id = i
                EXIT doens1a
             END IF
             pe_index = pe_index + 1
          END DO
       END DO doens1a
       IF (task_id==1) THEN
          t_id = MPI_UNDEFINED
       ELSE
          t_id = task_id
       END IF

       CALL MPI_Comm_split(COMM_pdaf, t_id, mype_pdaf, &
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
       ! *** Re-initialize PE information   ***
       ! *** according to model communicator ***
          CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
          CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

          IF (screen > 1) THEN
             WRITE (*,*) 'PDAF: MODEL - mype(w)= ', mype_pdaf, '; model task: ', task_id, &
                  '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
          END IF
       END IF

! ***            COMM_ASSIM               ***
! *** Generate communicator for analysis  ***

       IF (filterpe) THEN
          my_color = task_id
       ELSE
          my_color = MPI_UNDEFINED
       ENDIF

       CALL MPI_Comm_split(COMM_pdaf, my_color, mype_pdaf, &
            COMM_assim, MPIerr)

       ! *** Initialize PE information          ***
       ! *** according to coupling communicator ***
       IF (filterpe) THEN
          CALL MPI_Comm_Size(COMM_assim, npes_assim, MPIerr)
          CALL MPI_Comm_Rank(COMM_assim, mype_assim, MPIerr)

          IF (screen > 1) THEN
             WRITE (*,*) 'PDAF: ASSIM - mype(w)= ', mype_pdaf, '; mype(f)= ', mype_assim, '; npes(f)= ', npes_assim
          END IF
       ENDIF


       ! ***              COMM_COUPLE                 ***
       ! *** Generate communicators for communication ***
       ! *** between model and filter PEs             ***
       ! *** (Split COMM_PDAF)                    ***

       IF (modelpe) THEN
          color_couple = mype_model + 1
       ELSE
          color_couple = mype_assim + 1
       ENDIF

       CALL MPI_Comm_split(COMM_pdaf, color_couple, mype_pdaf, &
            COMM_couple, MPIerr)

       ! *** Initialize PE information          ***
       ! *** according to coupling communicator ***

       CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
       CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

! *** Display process configuration ***

       IF (screen > 0) THEN
          IF (mype_pdaf == 0) THEN
             WRITE (*, '(a13, 2x, a)') 'PDAF    Pconf', 'Process configuration:'
             WRITE (*, '(a13, 2x, a6, a9, a10, a14, 3x, a12, /a13, 2x, a5, a9, a7, a7, a7, a7, a8, /a13, 2x, a)') &
                  'PDAF    Pconf', 'world', 'assim', 'model', 'couple', 'model/assim', &
                  'PDAF    Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'M/A', &
                  'PDAF    Pconf', '----------------------------------------------------------'
          END IF
          CALL MPI_Barrier(COMM_pdaf, MPIerr)
          IF (filterpe) THEN
             WRITE (*, '(a, 2x, i4, 4x, i4, 4x, 3x, 4x, 3x, 4x, i3, 4x, i3, 5x, 3x, a1)') &
                  'PDAF    Pconf', mype_pdaf, mype_assim, color_couple, mype_couple, 'A'
          ENDIF
          IF (modelpe) THEN
             WRITE (*,'(a, 2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, 3x, a1)') &
                  'PDAF    Pconf', mype_pdaf, task_id, mype_model, color_couple, mype_couple, 'M'
          END IF
          CALL MPI_Barrier(COMM_pdaf, MPIerr)

          IF (mype_pdaf == 0) WRITE (*, '(/a)') ''

       END IF

! *** Store parallelization information for internal use ***

       COMM_model_mod = COMM_model
       COMM_filter_mod = COMM_assim
       task_id_mod = task_id
       n_modeltasks_mod = n_modeltasks

       ! Set flags
       isset_parallel = .TRUE.
       IF (COMM_PDAF /= MPI_COMM_WORLD) isset_comm_pdaf = .TRUE.

    END IF partype

  END SUBROUTINE PDAF3_init_parallel

END MODULE PDAF3init
