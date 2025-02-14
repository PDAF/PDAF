!> Control analysis update of the GLOBALTEMPLATE
!!
!! Routine to control the analysis update of the GLOBALTEMPLATE.
!! 
!! The analysis and ensemble tranformation are performed by
!! calling PDAF\_GLOBALTEMPLATE\_analysis. In addition, the routine
!! U\_prepoststep is called prior to the analysis and after
!! the resampling to allow the user to access the ensemble
!! information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAF_GLOBALTEMPLATE_update

CONTAINS
  SUBROUTINE PDAFGLOBALTEMPLATE_update(step, dim_p, dim_obs_p, dim_ens, &
       state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, &
       U_init_obs, U_prodRinvA, U_init_obsvar, U_prepoststep, &
       screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

! TEMPLATE: The first use-includes are generic
    USE PDAF_timer, &                ! Routines for timings
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_memcounting, &          ! Routine for memory counting
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &        ! Variables for parallelization
         ONLY: mype, dim_ens_l
    USE PDAF_analysis_utils, &       ! Routine for adaptive forgetting factor
         ONLY: PDAF_set_forget
    USE PDAFobs, &                   ! Routines and variables for observations
         ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
         observe_ens, HX_p, HXbar_p, obs_p
! TEMPLATE: Include here variables from module of the DA method
    USE PDAF_GLOBALTEMPLATE, &
         ONLY: localfilter, debug, forget, type_forget, type_trans
! TEMPLATE: Include here the name of the analysis routine
    USE PDAF_GLOBALTEMPLATE_analysis, &
         ONLY: PDAF_GLOBALTEMPLATE_ana
! TEMPLATE: If smoothing is use we include the routine from the smoother module
    USE PDAF_smoother, &             ! Name of generic smoothing routine
         ONLY: PDAF_smoothing

    IMPLICIT NONE

! +++ TEMPLATE: The argument variables are usually the same for each method

! *** Arguments ***
    INTEGER, INTENT(in) :: step          !< Current time step
    INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
    INTEGER, INTENT(out) :: dim_obs_p    !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens)!< Inverse of matrix U
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
    INTEGER, INTENT(in) :: screen        !< Verbosity flag
    INTEGER, INTENT(in) :: subtype       !< Filter subtype
    INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
    REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
    INTEGER, INTENT(inout) :: cnt_maxlag !< Count number of past time steps for smoothing
    INTEGER, INTENT(inout) :: flag       !< Status flag

! +++ TEMPLATE: The external subroutines used in a DA method are specific for the method

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_dim_obs, &      !< Initialize dimension of observation vector
         U_obs_op, &                   !< Observation operator
         U_init_obs, &                 !< Initialize observation vector
         U_init_obsvar, &              !< Initialize mean observation error variance
         U_prepoststep, &              !< User supplied pre/poststep routine
         U_prodRinvA                   !< Provide product R^-1 A for GLOBALTEMPLATE analysis

! *** local variables ***
    INTEGER :: i, j                    ! Counters
    INTEGER :: minusStep               ! Time step counter
    REAL :: forget_ana                 ! Forgetting factor actually used in analysis
    LOGICAL :: do_init_dim_obs         ! Flag for initializing dim_obs_p in PDAFobs_init


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
       CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
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
    CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
         state_p, Ainv, ens_p, flag)
    CALL PDAF_timeit(5, 'old')

    IF (mype == 0 .AND. screen > 1) THEN
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
       CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
            state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
            screen, debug, .true., do_init_dim_obs, .true., .true., .true.)

       CALL PDAF_timeit(3, 'old')
    END IF


! ***********************
! ***  Analysis step  ***
! ***********************

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
       WRITE (*, '(a, 1x, i7, 3x, a)') &
            'PDAF', step, 'Assimilating observations - GLOBALTEMPLATE'
    END IF

#ifndef PDAF_NO_UPDATE
    CALL PDAF_timeit(3, 'new')

! +++ TEMPLATE:
! +++ We recommend to include a list of the parameter values to support debugging
! +++ This needs to be adapted to the method-specific options

    IF (debug>0) THEN
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(3) dim_lag     ', dim_lag
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(4) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(5) type_forget ', type_forget
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(6) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(7) -not used-  '
       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_int(8) observe_ens ', observe_ens

       WRITE (*,*) '++ PDAF-debug PDAF_GLOBALTEMPLATE_update', debug, &
            'Configuration: param_real(1) forget     ', forget
    END IF


! *** Compute adaptive forgetting factor ***

! +++ TEMPLATE:
! +++ This is generic and could be kept if a globally-adaptive forgetting 
! +++ factor inflation is used in the DA method

    forget_ana = forget
    IF (type_forget == 1) THEN
       CALL PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, HX_p, &
            HXbar_p, obs_p, U_init_obsvar, forget, forget_ana, &
            screen)
    ELSE
       IF (mype == 0 .AND. screen > 0) THEN
          WRITE (*, '(a, 5x, a, F7.2)') &
               'PDAF', '--- apply multiplicative inflation with fixed forget', forget
       END IF
    END IF


! ***  Execute Analysis step  ***

! +++ TEMPLATE:
! +++ The call to PDAF_GLOBALTEMPLATE_analysis should be adapted for
! +++ the call-back routines and variables used by the method

    CALL PDAF_GLOBALTEMPLATE_ana(step, dim_p, dim_obs_p, dim_ens, &
         state_p, Ainv, ens_p, &
         HX_p, HXbar_p, obs_p, &
         forget_ana, U_prodRinvA, type_trans, screen, debug, flag)

! +++ TEMPLATE:
! +++ The smoother routine is generic. It can be used
! +++ as long as sens_l(:,lag) * (forget*Ainv + diag(invdimens))
! +++ yields the smoother update for the ensemble at lag 'lag'. 

    ! *** Perform smoothing of past ensembles ***
    IF (dim_lag>0) THEN
       CALL PDAF_timeit(15, 'new')
       CALL PDAF_timeit(51, 'new')
       CALL PDAF_smoothing(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
            cnt_maxlag, forget_ana, screen)
       CALL PDAF_timeit(51, 'old')
       CALL PDAF_timeit(15, 'old')
    END IF

    CALL PDAF_timeit(3, 'old')

    IF (mype == 0 .AND. screen > 1) THEN
       WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
            'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
    END IF

#else
    WRITE (*,'(/5x,a/)') &
         '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif
    
! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the analysis
! +++ ensemble is standard and should be kept

! *** Poststep for analysis ensemble ***
    CALL PDAF_timeit(5, 'new')
    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
       WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
    ENDIF
    CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
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

! +++ TEMPLATE:
! +++ This call to PDAFobs_dealloc is mandatory; do not change it

    ! Deallocate observation arrays
    CALL PDAFobs_dealloc()

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_GLOBALTEMPLATE_update -- END'

  END SUBROUTINE PDAFGLOBALTEMPLATE_update

END MODULE PDAF_GLOBALTEMPLATE_update
