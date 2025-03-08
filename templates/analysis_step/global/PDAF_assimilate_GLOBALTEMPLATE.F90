!> Interface to PDAF for global filter for fully-parallel implementation
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! analysis is computed by calling PDAF_put_state_GLOBALTEMPLATE.
!! Subsequently, PDAF_get_state is called to initialize
!! the next forecast phase. 
!!
!! The code is generic. The only part specific for a DA method are
!! the call to the routine PDAF\_put\_state\_X where the
!! analysis is computed and the names of the user-supplied 
!! call-back routines. These specific call-back subroutines 
!! are specified in the calls to the two core routines.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a filter, the only required changes to this routine
!! should be
!! - replace 'GLOBALTEMPLATE' by the name of the new method
!! - adapt the argument lists in PDAF\_assimilate\_GLOBALTEMPLATE
!!   and PDAF\_put\_state\_GLOBALTEMPLATE
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
MODULE PDAFassimilate_GLOBALTEMPLATE

CONTAINS

  SUBROUTINE PDAF_assimilate_GLOBALTEMPLATE(U_collect_state, U_distribute_state, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
       U_init_obsvar, U_next_observation, U_prepoststep, outflag)

    USE PDAF_mod_core,   &            ! Variables for framework functionality
         ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
    USE PDAF_mod_parallel, &          ! Variables for parallelization
         ONLY: mype_world
    USE PDAF_forecast, &              ! Routine for operations during forecast phase
         ONLY: PDAF_fcst_operations
    USE PDAFput_state_GLOBALTEMPLATE, &    ! Put_state routine for this DA method
         ONLY: PDAF_put_state_GLOBALTEMPLATE

    IMPLICIT NONE

! TEMPLATE: 'outflag' is standard and should be kept

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag
  
! TEMPLATE: The external subroutines depends on the DA method and should be adapted

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
    ! Routines for ensemble framework
    EXTERNAL :: U_collect_state, & !< Write model fields into state vector
         U_next_observation, &     !< Provide time step, time and dimension of next observation
         U_distribute_state, &     !< Write state vector into model fields
         U_prepoststep             !< User supplied pre/poststep routine
    ! Observation-related routines for analysis step
    EXTERNAL :: U_init_dim_obs, &  !< Initialize dimension of observation vector
         U_obs_op, &               !< Observation operator
         U_init_obs, &             !< Initialize observation vector
         U_init_obsvar, &          ! Initialize mean observation error variance
         U_prodRinvA               !< Provide product R^-1 A

! TEMPLATE: The local variables are usually generic and don't need changes

! *** Local variables ***
    INTEGER :: steps     ! Number of time steps in next forecast phase
    INTEGER :: doexit    ! Exit flag; not used in this variant
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


! ********************************
! *** At end of forecast phase ***
! ********************************

! TEMPLATE: Below the only non-generic part is the call to
! PDAF_put_state_GLOBALTEMPLATE. Other lines should not be changed.

    IF (cnt_steps == nsteps) THEN

       IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

       ! Set flag for assimilation
       assim_flag = 1

       ! *** Call analysis step ***

! TEMPLATE: Specific call for DA method
       CALL PDAF_put_state_GLOBALTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op, &
            U_init_obs, U_prodRinvA, U_init_obsvar, U_prepoststep, outflag)

       ! *** Prepare start of next ensemble forecast ***

       IF (outflag==0) THEN
          CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
               U_prepoststep, outflag)
       END IF

       nsteps = steps

    ELSE
       assim_flag = 0
       outflag = 0
    END IF

  END SUBROUTINE PDAF_assimilate_GLOBALTEMPLATE

END MODULE PDAFassimilate_GLOBALTEMPLATE
