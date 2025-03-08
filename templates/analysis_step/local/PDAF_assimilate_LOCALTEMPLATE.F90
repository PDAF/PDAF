!> Interface to PDAF for local filter for fully-parallel implementation
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! analysis is computed by calling PDAF_put_state_LOCALTEMPLATE.
!! Subsequently, PDAF_get_state is called to initialize
!! the next forecast phase. 
!!
!! The code is generic. The only filter-specific part is
!! the call to the routine PDAF\_put\_state\_X where the
!! analysis is computed. The filter-specific call-back subroutines 
!! are specified in the calls to the two core routines.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a filter, the only required changes to this routine
!! should be
!! - replace 'LOCALTEMPLATE' by the name of the new method
!! - potentially adapt the argument lists in PDAF\_assimilate\_LOCALTEMPLATE
!!   and PDAF\_put\_state\_LOCALTEMPLATE
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on LETKF
!! * Later revisions - see repository log
!!
MODULE PDAFassimilate_LOCALTEMPLATE

CONTAINS

  SUBROUTINE PDAF_assimilate_LOCALTEMPLATE(U_collect_state, U_distribute_state, &
       U_init_dim_obs, U_obs_op, U_init_obs, &
       U_init_n_domains_p, U_init_dim_l, U_g2l_state, U_l2g_state, &
       U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, U_prodRinvA_l, &
       U_init_obsvar, U_init_obsvar_l, U_next_observation, U_prepoststep, outflag)

    USE PDAF_mod_core, &              ! Variables for framework functionality
         ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
    USE PDAF_mod_parallel, &          ! Variables for parallelization
         ONLY: mype_world
    USE PDAF_forecast, &              ! Routine for operations during forecast phase
         ONLY: PDAF_fcst_operations

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
! PDAF_put_state_GLOBALTEMPLATE. Other lines should not be changed.

    IF (cnt_steps == nsteps) THEN

       IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

       ! Set flag for assimilation
       assim_flag = 1

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

       ! Set flag that no analysis was done
       assim_flag = 0  

       outflag = 0
    END IF

  END SUBROUTINE PDAF_assimilate_LOCALTEMPLATE

END MODULE PDAFassimilate_LOCALTEMPLATE
