!$Id$
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

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
  USE pdaf !, &   ! Interface definitions to PDAF core routines
!        ONLY: PDAFomi_put_state_local, PDAFomi_put_state_global, &
!        PDAFomi_put_state_lenkf, PDAFomi_put_state_generate_obs, &
!        PDAF_localfilter, &
!        PDAF_get_state, PDAF_get_assim_flag
  USE mod_parallel, &             ! Parallelization variables
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
!  REAL, INTENT(INOUT) :: time  ! Model time

! ! External subroutines 
! !  (subroutine names are passed over to PDAF in the calls to 
! !  PDAF_get_state and PDAF_put_state_X. This allows the user 
! !  to specify the actual name of a routine. However, the 
! !  PDAF-internal name of a subroutine might be different from
! !  the external name!)
!
  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &  ! Collect a state vector from model fields
       distribute_state_pdaf, &      ! Distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_pdaf              ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &     ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi        ! Apply localization to covariance matrix in LEnKF
! ! Subroutine used for generating observations
  EXTERNAL :: get_obs_f_pdaf         ! Get vector of synthetic observations from PDAF

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state
! Calls: integration
! Calls: PDAF_put_state_seek
! Calls: PDAF_put_state_seik
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: i           ! Counter
  INTEGER :: nsteps      ! Number of time steps to be performed in current forecast
  INTEGER :: doexit      ! Whether to exit forecasting (1=true)
  INTEGER :: status      ! Status flag for filter routines
  REAL :: timenow        ! Current model time
  INTEGER :: localfilter ! Flag for domain-localized filter (1=true)
  INTEGER :: assim

! *************************
! *** Perform forecasts ***
! *************************

  ! *** PDAF: Get state and forecast information (nsteps,time)  ***
!   CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
!        distribute_state_pdaf, prepoststep_pdaf, status)

  ! PDAF: External loop around model time stepper loop
!  pdaf_modelloop: DO  


     ! Check whether forecast has to be performed
!     checkforecast: IF (doexit /= 1 .AND. status == 0) THEN

        ! *** Forecast ensemble state ***

!        DO i = 1, nsteps

!            ! Initialize current model time
!            time = timenow
! 
!            ! *** call time stepper ***  
!            CALL integration(time, 1)


           ! *** PDAF: Perform assimilation if ensemble forecast is completed   ***

  IF (PDAF_localfilter() == 1) THEN
     CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, next_observation_pdaf, status)
  ELSE
     CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
          next_observation_pdaf, status)
  END IF

!END DO

!         ! *** PDAF: Get state and forecast information (nsteps,time)  ***
!         CALL PDAF_get_fcst_info(nsteps, timenow, doexit)
! 
!      ELSE checkforecast
! 
!         ! *** No more work, exit modeling loop
!         EXIT pdaf_modelloop
! 
!      END IF checkforecast
! 
!   END DO pdaf_modelloop


! ************************
! *** Check error flag ***
! ************************

  IF (status /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
