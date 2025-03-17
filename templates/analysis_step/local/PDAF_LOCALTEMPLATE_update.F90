!>Control analysis update of LOCALTEMPLATE
!!
!! This routine prepares the actual analysis update which
!! is computed then by calling PDAF_analysis_LOCALTEMPLATE.
!!
!! The analysis is performed by first preparing several
!! global quantities on the PE-local domain, like the
!! observed part of the state ensemble for all local
!! analysis domains on the PE-local state domain.
!! Then, the analysis (PDAF\_LOCALTEMPLATE\_analysis) is performed
!! within a loop over all local analysis domains in the PE-local 
!! state domain. In this loop, the local state and 
!! observation dimensions are initialized and the global 
!! state ensemble is restricted to the local analysis domain.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling outside of
!! the loop over the local domains to allow the user
!! to access the ensemble information.
!!
!! ADAPTING THE TEMPLATE
!! The structure of the operations included in this template is
!! typical for a (domain)-local ensemble DA method. One should
!! check if other operations are required in the global preparations
!! and the local steps which are included in this template. Further,
!! the call to PDAF_LOCALTEMPLATE_analysis and the particular
!! arguments of this subroutine need to be adapted. For the
!! different operations see the particular comments in the code.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on LETKF
!! * Later revisions - see repository log
!!
MODULE PDAF_LOCALTEMPLATE_update

CONTAINS
  SUBROUTINE PDAFLOCALTEMPLATE_update(step, dim_p, dim_obs_f, dim_ens, &
       state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       U_init_n_domains_p, U_init_dim_l, U_g2l_state, U_l2g_state, &
       U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, U_prodRinvA_l, &
       U_init_obsvar, U_init_obsvar_l, U_prepoststep, &
       screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

    USE PDAF_timer, &                 ! Routines for timings
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_memcounting, &           ! Routine for memory counting
         ONLY: PDAF_memcount
    USE PDAF_mod_parallel, &          ! Variables for parallelization
         ONLY: mype, dim_ens_l
    USE PDAF_analysis_utils, &        ! Routine for adaptive forgetting factor
         ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, &
         PDAF_incr_local_obsstats, PDAF_print_local_obsstats, &
         PDAF_set_forget, PDAF_set_forget_local
    USE PDAFobs, &                    ! Routines and variables for observations
         ONLY: PDAFobs_init, PDAFobs_init_local, PDAFobs_dealloc, PDAFobs_dealloc_local, &
         type_obs_init, observe_ens, HX_f => HX_p, HXbar_f => HXbar_p, obs_f => obs_p, &
         HX_l, HXbar_l, obs_l
    USE PDAFomi_obs_f, &              ! PDAF-OMI variables
         ONLY: omi_n_obstypes => n_obstypes, omi_obs_diag => obs_diag
! TEMPLATE: Include here variables from module of the DA method
    USE PDAF_LOCALTEMPLATE, &
         ONLY: localfilter, debug, forget, type_forget, &
         type_trans, inloop, forget_l, member_save
! TEMPLATE: Include here the name of the analysis routine
    USE PDAF_LOCALTEMPLATE_analysis, &
         ONLY: PDAF_LOCALTEMPLATE_ana
! TEMPLATE: If smoothing is use we  include the routine from the smoother module
    USE PDAF_smoother, &             ! Name of generic smoothing routine
         ONLY: PDAF_smoothing_local

    IMPLICIT NONE

! *** Arguments ***
! Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
!    suffix _f: Denotes a full variable of all observations required for the
!               analysis loop on the PE-local domain
    INTEGER, INTENT(in) :: step          ! Current time step
    INTEGER, INTENT(in) :: dim_p         ! PE-local dimension of model state
    INTEGER, INTENT(out) :: dim_obs_f    ! PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p)         ! PE-local model state
    REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens) ! Transform matrix
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)  ! PE-local ensemble matrix
    INTEGER, INTENT(in) :: screen        ! Verbosity flag
    INTEGER, INTENT(in) :: subtype       ! Filter subtype
    INTEGER, INTENT(in) :: dim_lag       ! Number of past time instances for smoother
    REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
    INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
    INTEGER, INTENT(inout) :: flag       ! Status flag

! ** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
    ! Routine for ensemble framework
    EXTERNAL :: U_prepoststep         !< User supplied pre/poststep routine
    ! Observation-related routines for analysis step
    EXTERNAL :: U_init_dim_obs, &     !< Initialize dimension of observation vector
         U_obs_op, &                  !< Observation operator
         U_init_dim_obs_l, &          !< Initialize dim. of obs. vector for local ana. domain
         U_init_obs, &                !< Initialize observation vector
         U_init_obs_l, &              !< Init. observation vector on local analysis domain
         U_g2l_obs, &                 !< Restrict full obs. vector to local analysis domain
         U_prodRinvA_l                !< Provide product R^-1 A on local analysis domain
    ! Routines for state localization
    EXTERNAL :: U_init_n_domains_p, & !< Provide number of local analysis domains
         U_init_dim_l, &              !< Init state dimension for local ana. domain
         U_init_obsvar, &             !< Initialize mean observation error variance
         U_init_obsvar_l, &           !< Initialize local mean observation error variance
         U_g2l_state, &               !< Get state on local ana. domain from full state
         U_l2g_state                  !< Init full state from state on local analysis domain

! *** Local variables ***
    INTEGER :: i, j, member           ! Counters
    INTEGER :: minusStep              ! Time step counter
    INTEGER :: domain_p               ! Counter for local analysis domain
    INTEGER :: n_domains_p            ! number of PE-local analysis domains
    REAL    :: forget_ana_l           ! forgetting factor supplied to analysis routine
    REAL    :: forget_ana             ! Possibly globally adaptive forgetting factor
    LOGICAL :: do_init_dim_obs        ! Flag for initializing dim_obs_p in PDAFobs_init
    INTEGER, SAVE :: allocflag = 0    ! Flag whether first time allocation is done
    REAL, ALLOCATABLE :: rndmat(:,:)  ! random rotation matrix for ensemble trans.
    ! Variables on local analysis domain
    INTEGER :: dim_l                  ! State dimension on local analysis domain
    INTEGER :: dim_obs_l              ! Observation dimension on local analysis domain
    REAL, ALLOCATABLE :: ens_l(:,:)   ! State ensemble on local analysis domain
    REAL, ALLOCATABLE :: state_l(:)   ! Mean state on local analysis domain
    REAL, ALLOCATABLE :: Ainv_l(:,:)  ! thread-local matrix Ainv


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

! +++ TEMPLATE:
! +++ For fixed-ensemble cases (like Ensemble OI) only the
! +++ central state is integrated by the model. Here, we
! +++ then need to add the ensemble perturbations

    CALL PDAF_timeit(3, 'new')
    CALL PDAF_timeit(51, 'new')

    fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
       ! Add mean/central state to ensemble members
       DO j = 1, dim_ens
          DO i = 1, dim_p
             ens_p(i, j) = ens_p(i, j) + state_p(i)
          END DO
       END DO
    END IF fixed_basis

    CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

! +++ TEMPLATE: The observation can be initialized here before
! +++ the call to U_prepoststep, or afterwards (see below)
! +++ The observation arrays (obs_p, HX_p, HXbar_p) are declared
! +++ in the module PDAFobs. The routine PDAFobs_initialize
! +++ allocates and fills theses arrays. In addition it can 
! +++ compute the ensemble mean state in state_p. The five logical 
! +++ options at the end of he argument list define which steps
! +++ are executed in the subroutine (see subroutine for its
! +++ documentation).

    IF (type_obs_init==0 .OR. type_obs_init==2) THEN
       ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
       ! It also compute the ensemble mean and stores it in state_p
       CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
            state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
            screen, debug, .true., .true., .true., .true., .true.)
    END IF
    CALL PDAF_timeit(3, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the forecast
! +++ ensemble is standard and should be kept

    CALL PDAF_timeit(5, 'new')
    minusStep = -step  ! Indicate forecast by negative time step number

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
       WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
    ENDIF
    CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
         state_p, Ainv, ens_p, flag)
    CALL PDAF_timeit(5, 'old')

    IF (mype == 0 .AND. screen > 0) THEN
       IF (screen > 1) &
            WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
            'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
    END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

! +++ The observation can be initialized here after the
! +++ call to U_prepoststep

    IF (type_obs_init>0) THEN
       CALL PDAF_timeit(3, 'new')

       IF (type_obs_init==1) THEN
          do_init_dim_obs=.true.
       ELSE
          ! Skip call to U_init_dim_obs when also called before prepoststep
          do_init_dim_obs=.false.   
       END IF

       ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
       ! It also compute the ensemble mean and stores it in state_p
       CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
            state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
            screen, debug, .true., do_init_dim_obs, .true., .true., .true.)

       CALL PDAF_timeit(3, 'old')
    END IF


! **************************************
! *** Preparation for local analysis ***
! **************************************

    IF (mype == 0 .AND. screen > 0) &
         WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

! +++ TEMPLATE: 
! +++ Before entering the local analysis loop we here initialize
! +++ - the number of local analysis domains
! +++ - potentially an adaptive value of the forgetting factor
! +++ - if requested a random rotation matrix

#ifndef PDAF_NO_UPDATE
    CALL PDAF_timeit(3, 'new')
    CALL PDAF_timeit(7, 'new')

! +++ TEMPLATE:
! +++ We recommend to include a list of the parameter values to support debugging
! +++ This needs to be adapted to the method-specific options

    IF (debug>0) THEN
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(3) dim_lag     ', dim_lag
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(4) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(5) type_forget ', type_forget
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(6) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(7) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_int(8) observe_ens           ', observe_ens

       WRITE (*,*) '++ PDAF-debug PDAF_LOCALTEMPLATE_update', debug, &
            'Configuration: param_real(1) forget     ', forget
    END IF

    ! *** Query number of analysis domains for the local analysis
    ! *** in the process-local domain
    CALL PDAF_timeit(42, 'new')
    CALL U_init_n_domains_p(step, n_domains_p)
    CALL PDAF_timeit(42, 'old')
  
    IF (screen > 0) THEN
       IF (mype == 0) THEN

! +++ TEMPLATE
! +++ If other subtypes exist, add them here
          IF (subtype == 0 .OR. subtype == 2) THEN
             WRITE (*, '(a, i7, 3x, a)') &
                  'PDAF ', step, 'Assimilating observations - LOCALTEMPLATE default analysis'
          END IF
       END IF
       IF (screen<3) THEN
          CALL PDAF_print_domain_stats(n_domains_p)
       ELSE
          WRITE (*, '(a, 5x, a, i6, a, i10)') &
               'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
       END IF
    END IF


! *** Local analysis: initialize global quantities ***

    CALL PDAF_timeit(51, 'new')

    ! *** Set forgetting factor globally 

! +++ TEMPLATE:
! +++ This is generic and could be kept if a globally-adaptive forgetting 
! +++ factor inflation is used in the DA method

    forget_ana = forget
    IF (type_forget == 1) THEN
       CALL PDAF_set_forget(step, localfilter, dim_obs_f, dim_ens, HX_f, &
            HXbar_f, obs_f, U_init_obsvar, forget, forget_ana, &
            screen)
    ELSE IF (type_forget == 0) THEN
       IF (mype == 0 .AND. screen > 0) THEN
          WRITE (*, '(a, 5x, a, F7.2)') &
               'PDAF', '--- apply multiplicative inflation with fixed forget', forget
       END IF
    ENDIF

    ! *** Initialize random transformation matrix

    CALL PDAF_timeit(33, 'new')
    ALLOCATE(rndmat(dim_ens, dim_ens))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

! +++ TEMPLATE: 
! +++ One might want to apply a random rotation to Ainv.
! +++ PDAF_general_rndmat provides a matrix with such random rotation
! +++ which is used in this example. For the local filter, the random
! +++ matrix is initialized here and then consistently used for all
! +++ local analysis domains

    IF (type_trans == 2) THEN
       ! Initialize random matrix
       CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)
    ELSE
       ! Initialize identity matrix
       rndmat = 0.0
       DO i = 1, dim_ens
          rndmat(1,1) = 1.0
       END DO
    END IF

    CALL PDAF_timeit(33, 'old')
    CALL PDAF_timeit(51, 'old')

    CALL PDAF_timeit(7, 'old')


! ******************************
! *** Perform local analysis ***
! ******************************

! +++ TEMPLATE: 
! +++ The lines before the start of the loop are usually kept unchanged
! +++ Do not change the OpenMP (!$OMP) lines unless you add shared variables 
! +++ that need to be private to a thread. Other changes might break
! +++ the OpenMP-parallelization.

    CALL PDAF_timeit(8, 'new')

    ! Initialize counters for statistics on local observations
    CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, Ainv_l, flag, forget_ana_l)

    ! Allocate ensemble transform matrix
    ALLOCATE(Ainv_l(dim_ens, dim_ens))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
    Ainv_l = 0.0

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
    localanalysis: DO domain_p = 1, n_domains_p

! +++ TEMPLATE: 
! +++ Inside the local analysis domain the typical operations are
! +++ - Determine the local state dimension (using U_init_dim_l)
! +++ - Initialize the local state ensemble and mean state (using U_g2l_state)
! +++ - Initialize local observations, observed ensemble and mean state (using PDAFobs_init_local)
! +++ - optionally set adaptive forgetting factor (using PDAFomi_init_local)
! +++ - Call the analysis routine of the local DA method (obviously method-specific)
! +++ - Initialize the global state ensemlbe from the local analysis states (using U_l2g_state)
! +++ - Potentially apply smoothing (e.g. using PAF_smoothing local)
! +++ - In addition, the code collect statistics about local observations (PDAF_*_local_obsstats)

! +++ TEMPLATE: The first part before the call to the analysis are usually generic and not changed

       ! Set flag that we are in the local analysis loop
       inloop = .true.

       ! Set forgetting factor to global standard value
       forget_l = forget_ana


       ! *************************************
       ! *** Initialize local state vector ***
       ! *************************************

       ! local state dimension
       CALL PDAF_timeit(45, 'new')
       CALL U_init_dim_l(step, domain_p, dim_l)
       CALL PDAF_timeit(45, 'old')

       ! Allocate arrays for local analysis domain
       ALLOCATE(ens_l(dim_l, dim_ens))
       ALLOCATE(state_l(dim_l))

       CALL PDAF_timeit(10, 'new')

       ! Get local state ensemble and mean state
       DO member = 1, dim_ens
          ! Store member index to make it accessible with PDAF_get_memberid
          member_save = member

          CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
               ens_l(:, member))
       END DO

       ! Store member index to make it accessible with PDAF_get_memberid
       member_save = 0

       CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
            state_l)

       CALL PDAF_timeit(10, 'old')


       ! *******************************************
       ! *** Initialize local observation arrays ***
       ! *******************************************

       CALL PDAF_timeit(11, 'new')
       CALL PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f, dim_ens, &
            U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, debug)
       CALL PDAF_timeit(11, 'old')


       ! ************************************************
       ! *** Compute local adaptive forgetting factor ***
       ! ************************************************

       ! Reset forget (can be reset with PDAF_reset_forget)
       forget_ana_l = forget_l

       IF (type_forget == 2 .AND. dim_obs_l > 0) THEN
          CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, &
               HX_l, HXbar_l, obs_l, U_init_obsvar_l, forget, forget_ana_l)
       ENDIF

! +++ TEMPLATE:
! +++ At this point we initialize the local ensemble, local ensemble mean 
! +++ as well as the local vector of observations, and local observed ensemble
! +++ and ensemble mean. These are then inputs to the local analysis step

       ! *********************
       ! *** Analysis step ***
       ! *********************

       ! Gather statistical information on local observations
       CALL PDAF_incr_local_obsstats(dim_obs_l)
     
       CALL PDAF_timeit(12, 'new')

       CALL PDAF_LOCALTEMPLATE_ana(domain_p, step, dim_l, dim_ens, &
            state_l, Ainv_l, ens_l, &
            dim_obs_l, HX_l, HXbar_l, obs_l, &
            rndmat, forget_ana_l, U_prodRinvA_l, &
            type_trans, screen, debug, flag)

       CALL PDAF_timeit(12, 'old')
       CALL PDAF_timeit(14, 'new')

       ! re-initialize full state ensemble on PE and mean state from local domain
       DO member = 1, dim_ens
          member_save = member

          CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
       END DO

       ! Initialize global mean state
       member_save = 0

       CALL U_l2g_state(step, domain_p, dim_l, state_l, dim_p, state_p)
    
       CALL PDAF_timeit(14, 'old')
       CALL PDAF_timeit(51, 'new')
       CALL PDAF_timeit(15, 'new')

! +++ TEMPLATE 
! +++ The smoother routine is generic. It can be used
! +++ as long as sens_l(:,lag) * (forget*Ainv + diag(invdimens))
! +++ yields the smoother update for the ensemble at lag 'lag'. 
! +++ The routine includes the projection onto the local
! +++ state vector and back.

       ! *** Perform smoothing of past ensembles ***
       CALL PDAF_smoothing_local(domain_p, step, dim_p, dim_l, dim_ens, &
            dim_lag, Ainv_l, ens_l, sens_p, cnt_maxlag, &
            U_g2l_state, U_l2g_state, forget_ana, screen)

       CALL PDAF_timeit(15, 'old')

       ! clean up
       DEALLOCATE(ens_l, state_l)
       CALL PDAFobs_dealloc_local()

       CALL PDAF_timeit(51, 'old')

    END DO localanalysis

! +++ TEMPLATE:
! +++ The code below is generic. Usually no changes are required.

    ! Set flag that we are not in the local analysis loop
    inloop = .false.

    CALL PDAF_timeit(51, 'new')

!$OMP CRITICAL
    ! Set Ainv - required for subtype=3
    Ainv = Ainv_l
!$OMP END CRITICAL

    DEALLOCATE(Ainv_l)
!$OMP END PARALLEL


    ! *** Print statistics for local analysis to the screen ***
    CALL PDAF_print_local_obsstats(screen)

    CALL PDAF_timeit(51, 'old')
    CALL PDAF_timeit(6, 'old')
    CALL PDAF_timeit(3, 'old')

    IF (mype == 0 .AND. screen > 0) THEN
       IF (screen > 1) THEN
          WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
               'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
       END IF
    END IF

! *** Clean up from local analysis update ***
    DEALLOCATE(rndmat)


! ******************************************************
! *** Initialize analysis observed ensemble and mean ***
! ******************************************************
! +++ TEMPLATE:
! +++ This additional call to PDAFomi_init initializes
! +++ the observed analysis ensemble and its mean. This
! +++ is used in the observation diagnostics of PDAF-OMI

  IF (omi_n_obstypes>0 .AND. omi_obs_diag>0) THEN
     ! This call initializes HX_p, HXbar_p in the module PDAFobs
     ! for the analysis ensemble
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

#else
    WRITE (*,'(/5x,a/)') &
         '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif


! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the analysis
! +++ ensemble is standard and should be kept

! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

    CALL PDAF_timeit(5, 'new')
    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
       WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
    ENDIF
    CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
         state_p, Ainv, ens_p, flag)
    CALL PDAF_timeit(5, 'old')
  
    IF (mype == 0 .AND. screen > 0) THEN
       IF (screen > 1) THEN
          WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
               'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
       END IF
       WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
    END IF


! ********************
! *** Finishing up ***
! ********************

    ! Set flag that allocation was already done once (used for memory counting)
    IF (allocflag == 0) allocflag = 1

! +++ TEMPLATE: The call to PDAFobs_dealloc is required to
! +++ deallocate the observation errors in the module PDAFobs

    ! Deallocate observation arrays
    CALL PDAFobs_dealloc()

  END SUBROUTINE PDAFLOCALTEMPLATE_update

END MODULE PDAF_LOCALTEMPLATE_update
