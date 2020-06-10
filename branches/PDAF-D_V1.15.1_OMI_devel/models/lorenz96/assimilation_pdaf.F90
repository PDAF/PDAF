!$Id$
!BOP
!
! !ROUTINE: assimilation_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE assimilation_pdaf(time)

! !DESCRIPTION:
! This routine performs the ensemble forcasts.
! PDAF with domain-decomposition is used.
!
! The model gets state vectors to be evolved as well as
! the number of time steps and the current model time 
! from PDAF by calling PDAF\_get\_state.
! Each forecasted state is written back into the ensemble 
! matrix of PDAF by calling a filter-specific routine
! PDAF\_put\_state\_X. When all ensemble members are 
! evolved and hence the forecast phase is completed, 
! PDAF\_put\_state\_X executes the analysis step of the
! chosen filter algorithm.
!
! In this routine, the real names of most of the 
! user-supplied routines for PDAF are specified (see below)
!
! !REVISION HISTORY:
! 2004-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &     ! Parallelization
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(INOUT) :: time  ! Model time

! ! External subroutines 
! !  (subroutine names are passed over to PDAF in the calls to 
! !  PDAF_get_state and PDAF_put_state_X. This allows the user 
! !  to specify the actual name of a routine. However, the 
! !  PDAF-internal name of a subroutine might be different from
! !  the external name!)
!
! ! Subroutines used with all filters
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, &
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       collect_state_pdaf, &           ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &            ! Initialize dimension of observation vector
       obs_op_pdaf, &                  ! Implementation of the Observation operator
       init_obs_pdaf, &                ! Routine to provide vector of measurements
       distribute_stateinc_pdaf, &     ! Routine to add state increment for IA
       prepoststep_pdaf                ! User supplied pre/poststep routine for SEIK
! ! Subroutine used in ESTKF/SEIK/ETKF/LESTKF/LSEIK/LETKF
  EXTERNAL :: init_obsvar_pdaf         ! Initialize mean observation error variance
! ! Subroutine used in ESTKF/SEIK/ETKF/SEEK
  EXTERNAL :: prodRinvA_pdaf           ! Provide product R^-1 A for some matrix A
! ! Subroutines used in EnKF
  EXTERNAL :: add_obs_error_pdaf, &    ! Add obs. error covariance R to HPH in EnKF
       init_obscovar_pdaf              ! Initialize obs error covar R in EnKF
! ! Subroutines used in LSEIK
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some matrix A (for LSEIK)
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &         ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf      ! Get dimension of full obs. vector for PE-local domain
! ! Subroutines used for localization in LEnKF
  EXTERNAL :: localize_covar_pdaf       ! Apply localization to HP and HPH^T
! ! Subroutines used in NETF
  EXTERNAL :: likelihood_pdaf      ! Compute observation likelihood for an ensemble member
! ! Subroutines used in LNETF
  EXTERNAL :: likelihood_l_pdaf  ! Compute local observation likelihood for an ensemble member
! ! Subroutine used for generating observations
  EXTERNAL :: get_obs_f_pdaf, & ! Get vector of synthetic observations from PDAF
       init_obserr_f_pdaf       ! Initialize vector of observation errors (standard deviations)

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state
! Calls: integration
! Calls: PDAF_put_state_seek
! Calls: PDAF_put_state_seik
! Calls: PDAF_put_state_enkf
! Calls: PDAF_put_state_lseik
! Calls: PDAF_put_state_etkf
! Calls: PDAF_put_state_letkf
! Calls: PDAF_put_state_estkf
! Calls: PDAF_put_state_lestkf
! Calls: PDAF_put_state_netf
! Calls: PDAF_put_state_lnetf
! Calls: PDAF_put_state_generate_obs
! Calls: PDAF_put_state_pf
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: nsteps    ! Number of time steps to be performed in current forecast
  INTEGER :: doexit    ! Whether to exit forecasting (1=true)
  INTEGER :: status    ! Status flag for filter routines
  REAL :: timenow      ! Current model time


! *************************
! *** Perform forecasts ***
! *************************

  ! PDAF: External loop around model time stepper loop
  pdaf_modelloop: DO  

     ! *** PDAF: Get state and forecast information (nsteps,time)  ***
     CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
          distribute_state_pdaf, prepoststep_pdaf, status)

     ! Check whether forecast has to be performed
     checkforecast: IF (doexit /= 1 .AND. status == 0) THEN

        ! *** Forecast ensemble state ***
      
        IF (nsteps > 0) THEN

           ! Initialize current model time
           time = timenow

           ! *** call time stepper ***  
           CALL integration(time, nsteps)

        END IF

        ! *** PDAF: Send state forecast to filter;                           ***
        ! *** PDAF: Perform assimilation if ensemble forecast is completed   ***
        ! *** PDAF: Distinct calls due to different name of analysis routine ***
        IF (filtertype == 1) THEN
           CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 2) THEN
           CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, add_obs_error_pdaf, &
                init_obscovar_pdaf, status)
        ELSE IF (filtertype == 3) THEN
           CALL PDAF_put_state_lseik(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        ELSE IF (filtertype == 4) THEN
           CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 5) THEN
           CALL PDAF_put_state_letkf(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        ELSE IF (filtertype == 6) THEN
           CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 7) THEN
           CALL PDAF_put_state_lestkf(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        ELSE IF (filtertype == 8) THEN
           CALL PDAF_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, localize_covar_pdaf, add_obs_error_pdaf, &
                init_obscovar_pdaf, status)
        ELSE IF (filtertype == 9) THEN
           CALL PDAF_put_state_netf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, likelihood_pdaf, status)
        ELSE IF (filtertype == 10) THEN
           CALL PDAF_put_state_lnetf(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_l_pdaf, prepoststep_pdaf, &
                likelihood_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, status)
        ELSE IF (filtertype == 11) THEN
           CALL PDAF_put_state_generate_obs(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obserr_f_pdaf, get_obs_f_pdaf, &
                prepoststep_pdaf, status)
        ELSE IF (filtertype == 12) THEN
           CALL PDAF_put_state_pf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_pdaf, likelihood_pdaf, status)
        END IF

     ELSE checkforecast

        ! *** No more work, exit modeling loop
        EXIT pdaf_modelloop

     END IF checkforecast

  END DO pdaf_modelloop


! ************************
! *** Check error flag ***
! ************************

  IF (status /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilation_pdaf
