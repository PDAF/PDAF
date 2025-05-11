!> Interface to PDAF for local filter for fully-parallel implementation
!!
!! Interface routines for oneline and offline coupling.
!!
!! The code is generic. The only part specific for a DA method
!! are the calls to the routines PDAF\_put\_state\_X and 
!! PDAF\_assim_offline_X where the analysis is computed, and
!! the names of the user-supplied call-back routines. These
!! specific call-back subroutines are specified in the calls
!! to the two core routines.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a filter, the only required changes to this routine
!! should be
!! - replace 'LOCALTEMPLATE' by the name of the new method
!! - potentially adapt the argument lists in PDAF\_assimilate\_LOCALTEMPLATE
!!   PDAF\_put\_state\_LOCALTEMPLATE and  PDAF\_assim_offline\_LOCALTEMPLATE
!!   PDAF\_LOCALTEMPLATE_update
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on LETKF
!! * Later revisions - see repository log
!!
MODULE PDAFassimilate_LOCALTEMPLATE

CONTAINS

!> Interface to PDAF analysis step for online coupling
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! forecast state is transferred to PDAF and the analysis
!! is computed by calling PDAF_put_state_LOCALTEMPLATE.
!! Subsequently, PDAF_get_state is called to initialize
!! the next forecast phase. 
!!
  SUBROUTINE PDAF_assimilate_LOCALTEMPLATE(U_collect_state, U_distribute_state, &
       U_init_dim_obs, U_obs_op, U_init_obs, &
       U_init_n_domains_p, U_init_dim_l, U_g2l_state, U_l2g_state, &
       U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, U_prodRinvA_l, &
       U_init_obsvar, U_init_obsvar_l, U_next_observation, U_prepoststep, outflag)

    USE PDAF_mod_core, &              ! Variables for framework functionality
         ONLY: cnt_steps, nsteps, assim_flag, reset_fcst_flag, use_PDAF_assim
    USE PDAF_mod_parallel, &          ! Variables for parallelization
         ONLY: mype_world
    USE PDAF_forecast, &              ! Routine for operations during forecast phase
         ONLY: PDAF_fcst_operations
    USE PDAFput_state_LOCALTEMPLATE, &  ! Put_state routine
         ONLY: PDAF_put_state_LOCALTEMPLATE

    IMPLICIT NONE

! TEMPLATE: 'outflag' is standard and should be kept

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag    !< Status flag
  
! TEMPLATE: The external subroutines depends on the DA method and should be adapted

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
    ! Routines for ensemble framework
    EXTERNAL :: U_collect_state, &     !< Write model fields into state vector
         U_next_observation, &         !< Provide time step, time and dimension of next observation
         U_distribute_state, &         !< Write state vector into model fields
         U_prepoststep                 !< User supplied pre/poststep routine
    ! Observation-related routines for analysis step
    EXTERNAL :: U_init_dim_obs, &      !< Initialize dimension of observation vector
         U_obs_op, &                   !< Observation operator
         U_init_obs, &                 !< Initialize PE-local observation vector
         U_init_obsvar, &              !< Initialize mean observation error variance
         U_init_dim_obs_l, &           !< Initialize dim. of obs. vector for local ana. domain
         U_init_obs_l, &               !< Init. observation vector on local analysis domain
         U_init_obsvar_l, &            !< Initialize local mean observation error variance
         U_prodRinvA_l, &              !< Provide product R^-1 A on local analysis domain
         U_g2l_obs                     !< Restrict full obs. vector to local analysis domain
    ! Routines for state localization
    EXTERNAL :: U_init_n_domains_p, &  !< Provide number of local analysis domains
         U_init_dim_l, &               !< Init state dimension for local ana. domain
         U_g2l_state, &                !< Get state on local ana. domain from full state
         U_l2g_state                   !< Init full state from state on local analysis domain

! TEMPLATE: The local variables are usually generic and don't need changes

! *** Local variables ***
    INTEGER :: steps     ! Number of time steps in next forecast phase
    INTEGER :: doexit    ! Exit flag; not used in PDAF_assimilate
    REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

! TEMPLATE: Generic - do not change

    ! Set flag for using PDAF_assimilate
    use_PDAF_assim = .TRUE.

    ! Increment time step counter
    cnt_steps = cnt_steps + 1

    ! *** Call generic routine for operations during time stepping.          ***
    ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

    CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
         U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! *****************************************************************
! *** At end of forecast phase call put_state for analysis step ***
! *****************************************************************

! TEMPLATE: Below the only non-generic part is the call to
! PDAF_put_state_LOCALTEMPLATE. Other lines should not be changed.

    IF (cnt_steps == nsteps) THEN

       ! Set flags for assimilation and forecast
       assim_flag = 0
       reset_fcst_flag = 1

       ! *** Call analysis step ***

! TEMPLATE: Specific call for DA method
       CALL PDAF_put_state_LOCALTEMPLATE(U_collect_state, &
            U_init_dim_obs, U_obs_op, U_init_obs, &
            U_init_n_domains_p, U_init_dim_l, U_g2l_state, U_l2g_state, &
            U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, U_prodRinvA_l, &
            U_init_obsvar, U_init_obsvar_l, U_prepoststep, outflag)

       ! *** Prepare start of next ensemble forecast ***

       IF (outflag==0) THEN
          CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
               U_prepoststep, outflag)
       END IF

       nsteps = steps

    ELSE
       ! *** During forecast phase ***
! TEMPLATE: These settings should not be changed
       ! Set flags for forecast
       assim_flag = 0  
       reset_fcst_flag = 0
       outflag = 0
    END IF

  END SUBROUTINE PDAF_assimilate_LOCALTEMPLATE


!-------------------------------------------------------------------------------
!> Interface to PDAF analysis step in offline coupling
!!
!! Interface routine called from the main program
!! for the PDAF offline mode.
!!
!! The code is very generic. Basically the only
!! filter-specific part is the call to the
!! update-routine PDAF\_X\_update where the analysis
!! is computed.  The filter-specific subroutines that
!! are specified in the call to PDAF\_assim\_offline\_X
!! are passed through to the update routine
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on put_state routine
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_assim_offline_LOCALTEMPLATE(U_init_dim_obs, U_obs_op, &
       U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
       U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
       U_init_obsvar, U_init_obsvar_l, outflag)

    USE PDAF_mod_core, &
         ONLY: dim_p, dim_ens, assim_flag, step_obs, &
         subtype_filter, screen, flag, offline_mode, &
         state, ens, Ainv, &
         sens, dim_lag, cnt_maxlag
    USE PDAF_mod_parallel, &
         ONLY: mype_world, filterpe
    USE PDAF_utils_filters, &
         ONLY: PDAF_configinfo_filters
    USE PDAFobs, &
         ONLY: dim_obs
    USE PDAF_LOCALTEMPLATE_update, &
         ONLY: PDAFLOCALTEMPLATE_update

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_dim_obs, &    !< Initialize dimension of observation vector
         U_obs_op, &                 !< Observation operator
         U_init_obs, &               !< Initialize observation vector
         U_init_obs_l, &             !< Init. observation vector on local analysis domain
         U_prodRinvA_l, &            !< Provide product R^-1 A on local analysis domain
         U_init_n_domains_p, &       !< Provide number of local analysis domains
         U_init_dim_l, &             !< Init state dimension for local ana. domain
         U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
         U_g2l_state, &              !< Get state on local ana. domain from full state
         U_l2g_state, &              !< Init full state from state on local analysis domain
         U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
         U_init_obsvar, &            !< Initialize mean observation error variance
         U_init_obsvar_l, &          !< Initialize local mean observation error variance
         U_prepoststep               !< User supplied pre/poststep routine

! *** local variables ***
    INTEGER :: i                     ! Counter


! *********************************************
! *** Perform analysis step in offline mode ***
! *********************************************

! TEMPLATE: Below the only non-generic part is the call to
! PDAFLOCALTEMPLATE_update. Other lines should not be changed.

    ! Set flag for assimilation
    assim_flag = 1

    ! Screen output
    IF (mype_world == 0 .AND. screen > 0) THEN
       ! Print configuration info (if not done before in PDAF_set_offline_mode)
       IF (.NOT.offline_mode) CALL PDAF_configinfo_filters(subtype_filter, 1)

       WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
       WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
       WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
    ENDIF

    ! Set flag for offline mode
    offline_mode = 1

    OnFilterPE: IF (filterpe) THEN
! TEMPLATE: Specific call for DA method
       CALL PDAFLOCALTEMPLATE_update(step_obs, dim_p, dim_obs, dim_ens, state, &
            Ainv, ens, U_init_dim_obs, U_obs_op, U_init_obs, &
            U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, &
            U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
            U_init_obsvar, U_init_obsvar_l, U_prepoststep, screen, &
            subtype_filter, dim_lag, sens, cnt_maxlag, flag)
    END IF OnFilterPE


! ********************
! *** finishing up ***
! ********************

    outflag = flag

  END SUBROUTINE PDAF_assim_offline_LOCALTEMPLATE

END MODULE PDAFassimilate_LOCALTEMPLATE
