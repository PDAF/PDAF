!$Id: assimilation_pdaf.F90 1248 2012-01-27 08:51:01Z lnerger $
!BOP
!
! !ROUTINE: assimilation_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE assimilation_pdaf(time)

! !DESCRIPTION:
! This routine performs the ensemble forecasts and includes
! the calls to PDAF to perform the analysis steps of the 
! data assimilation process.
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
! This implementation uses the full flexible interface of
! PDAF. Here, the real names of user-supplied routines are 
! specified (see code). This implementation allows the user
! to select the names of the user-supplied routines.
!
! This example is for the 1D dummy model. All 6 filters
! are implemented. If one only wants a single filter, one
! can delete all calls to PDAF_put_state_X except for that 
! for the chosen filter 'X'. 
!
! In this implementation we consider the case that the time
! stepper of the model is implemented as a subroutine. If this
! is not the case, the operations of this routine have to be
! added to the model code around the time stepping loop of
! the model.
!
! !REVISION HISTORY:
! 2004-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &     ! Parallelization
       ONLY: Comm_model, MPIerr, mype_world, abort_parallel
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
       distribute_stateinc_pdaf        ! Routine to add state increment for IAU
! ! Subroutines used in SEIK and ETKF
  EXTERNAL :: prepoststep_ens_pdaf, &  ! User supplied pre/poststep routine
       init_obsvar_pdaf                ! Initialize mean observation error variance
! ! Subroutine used in SEEK
  EXTERNAL :: prepoststep_seek_pdaf    ! User supplied pre/poststep routine
! ! Subroutines used in EnKF
  EXTERNAL :: add_obs_error_pdaf, &    ! Add obs. error covariance R to HPH in EnKF
       init_obscovar_pdaf              ! Initialize obs error covar R in EnKF
! ! Subroutine used in SEIK, ETKF, and SEEK
  EXTERNAL :: prodRinvA_pdaf           ! Provide product R^-1 A for some matrix A
! ! Subroutines used in LSEIK and LETKF
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state
! Calls: integration
! Calls: PDAF_put_state_seek
! Calls: PDAF_put_state_seik
! Calls: PDAF_put_state_enkf
! Calls: PDAF_put_state_lseik
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: nsteps    ! Number of time steps to be performed in current forecast
  INTEGER :: doexit    ! Whether to exit forecasting (1=true)
  INTEGER :: status    ! Status flag for filter routines
  REAL    :: timenow   ! Current model time


! *************************
! *** Perform forecasts ***
! *************************

  ! PDAF: External loop around model time stepper loop
  pdaf_ensembleloop: DO  

     ! *** PDAF: Get state and forecast information (nsteps,time)         ***
     ! *** PDAF: Distinct calls for ensemble-based filters and SEEK       ***
     IF (filtertype /= 0) THEN
        ! Ensemble-based filters (SEIK, EnKF, ETKF, LSEIK, LETKF)
        CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
             distribute_state_pdaf, prepoststep_ens_pdaf, status)
     ELSE 
        ! SEEK - distinct is only the routine prepoststep_seek
        CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
             distribute_state_pdaf, prepoststep_seek_pdaf, status)
     END IF

     ! Check whether forecast has to be performed
     checkforecast: IF (doexit /= 1 .AND. status == 0) THEN

        ! *** Forecast ensemble state ***
      
        IF (nsteps > 0) THEN

           ! Initialize current model time
           time = timenow

           ! *** call time stepper ***  
           CALL integration(time, nsteps)

        END IF

        ! *** PDAF: Send state forecast to filter;                                  ***
        ! *** PDAF: Perform assimilation if ensemble forecast is completed          ***
        ! *** PDAF: Distinct calls due to different names of user-supplied routines ***
        IF (filtertype == 0) THEN
           CALL PDAF_put_state_seek(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_seek_pdaf, prodRinvA_pdaf, status)
        ELSE IF (filtertype == 1) THEN
           CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 2) THEN
           CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_ens_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
                status)
        ELSE IF (filtertype == 3) THEN
           CALL PDAF_put_state_lseik(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        ELSE IF (filtertype == 4) THEN
           CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 5) THEN
           CALL PDAF_put_state_letkf(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        ELSE IF (filtertype == 6) THEN
           CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
                init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
        ELSE IF (filtertype == 7) THEN
           CALL PDAF_put_state_lestkf(collect_state_pdaf, init_dim_obs_f_pdaf, &
                obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
                prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
                init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
                g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
        END IF

        CALL MPI_barrier(COMM_model, MPIERR)

     ELSE checkforecast

        ! *** No more work, exit ensemble-forecast loop
        EXIT pdaf_ensembleloop

     END IF checkforecast

  END DO pdaf_ensembleloop


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
