!$Id: assimilation_pdaf.F90 1673 2016-12-02 10:46:31Z lnerger $
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
  EXTERNAL :: next_observation, & ! Provide time step, model time, &
                                  ! and dimension of next observation
       distribute_state, &        ! Routine to distribute a state vector to model fields
       collect_state, &           ! Routine to collect a state vector from model fields
       init_dim_obs, &            ! Initialize dimension of observation vector
       obs_op, &                  ! Implementation of the Observation operator
       init_obs, &                ! Routine to provide vector of measurements
       distribute_stateinc        ! Routine to add state increment for IAU
! ! Subroutine used in SEIK
  EXTERNAL :: prepoststep_seik, & ! User supplied pre/poststep routine for SEIK
       init_obsvar                ! Initialize mean observation error variance
! ! Subroutine used in SEIK and SEEK
  EXTERNAL :: prodRinvA           ! Provide product R^-1 A for some matrix A
! ! Subroutines used in EnKF
  EXTERNAL :: add_obs_error, &    ! Add obs. error covariance R to HPH in EnKF
       init_obscovar              ! Initialize obs error covar R in EnKF
! ! Subroutines used in LSEIK
  EXTERNAL :: init_n_domains, &   ! Provide number of local analysis domains
       init_dim_local, &      ! Initialize state dimension for local ana. domain
       init_dim_obs_local,&   ! Initialize dim. of obs. vector for local ana. domain
       global2local_state, &  ! Get state on local ana. domain from global state
       local2global_state, &  ! Init global state from state on local analysis domain
       global2local_obs, &    ! Restrict a global obs. vector to local analysis domain
       init_obs_local, &      ! Provide vector of measurements for local ana. domain
       prodRinvA_local, &     ! Provide product R^-1 A for some matrix A (for LSEIK)
       init_obsvar_local, &   ! Initialize local mean observation error variance
       init_obs_full, &       ! Provide full vector of measurements for PE-local domain
       obs_op_full, &         ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_full      ! Get dimension of full obs. vector for PE-local domain
! ! Subroutines used in NETF
  EXTERNAL :: likelihood      ! Compute observation likelihood for an ensemble member
! ! Subroutines used in LNETF
  EXTERNAL :: likelihood_local  ! Compute local observation likelihood for an ensemble member

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
     CALL PDAF_get_state(nsteps, timenow, doexit, next_observation, &
          distribute_state, prepoststep_seik, status)

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
           CALL PDAF_put_state_seik(collect_state, init_dim_obs, obs_op, &
                init_obs, prepoststep_seik, prodRinvA, init_obsvar, status)
        ELSE IF (filtertype == 2) THEN
           CALL PDAF_put_state_enkf(collect_state, init_dim_obs, obs_op, &
                init_obs, prepoststep_seik, add_obs_error, &
                init_obscovar, status)
        ELSE IF (filtertype == 3) THEN
           CALL PDAF_put_state_lseik(collect_state, init_dim_obs_full, &
                obs_op_full, init_obs_full, init_obs_local, prepoststep_seik, &
                prodRinvA_local, init_n_domains, init_dim_local, &
                init_dim_obs_local, global2local_state, local2global_state, &
                global2local_obs, init_obsvar, init_obsvar_local, status)
        ELSE IF (filtertype == 4) THEN
           CALL PDAF_put_state_etkf(collect_state, init_dim_obs, obs_op, &
                init_obs, prepoststep_seik, prodRinvA, init_obsvar, status)
        ELSE IF (filtertype == 5) THEN
           CALL PDAF_put_state_letkf(collect_state, init_dim_obs_full, &
                obs_op_full, init_obs_full, init_obs_local, prepoststep_seik, &
                prodRinvA_local, init_n_domains, init_dim_local, &
                init_dim_obs_local, global2local_state, local2global_state, &
                global2local_obs, init_obsvar, init_obsvar_local, status)
        ELSE IF (filtertype == 6) THEN
           CALL PDAF_put_state_estkf(collect_state, init_dim_obs, obs_op, &
                init_obs, prepoststep_seik, prodRinvA, init_obsvar, status)
        ELSE IF (filtertype == 7) THEN
           CALL PDAF_put_state_lestkf(collect_state, init_dim_obs_full, &
                obs_op_full, init_obs_full, init_obs_local, prepoststep_seik, &
                prodRinvA_local, init_n_domains, init_dim_local, &
                init_dim_obs_local, global2local_state, local2global_state, &
                global2local_obs, init_obsvar, init_obsvar_local, status)
        ELSE IF (filtertype == 9) THEN
           CALL PDAF_put_state_netf(collect_state, init_dim_obs, obs_op, &
                init_obs, prepoststep_seik, likelihood, status)
        ELSE IF (filtertype == 10) THEN
           CALL PDAF_put_state_lnetf(collect_state, init_dim_obs_full, &
                obs_op_full, init_obs_local, prepoststep_seik, &
                likelihood_local, init_n_domains, init_dim_local, &
                init_dim_obs_local, global2local_state, local2global_state, &
                global2local_obs, status)
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
