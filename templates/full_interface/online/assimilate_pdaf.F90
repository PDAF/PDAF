!$Id$
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

! !DESCRIPTION:
! This routine is called during the model integrations at each time 
! step. It check whether the forecast phase is completed. If so, 
! PDAF_put_state_X is called to perform the analysis step.
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF                     ! PDAF 
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype
  USE PDAF                     ! Include PDAF calls

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! Calls: PDAF_assimilate_X
!EOP

! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag


! ! External subroutines
! !  (subroutine names are passed over to PDAF in the calls to 
! !  PDAF_get_state and PDAF_assimilate_X. This allows the user 
! !  to specify the actual name of a routine. However, the 
! !  PDAF-internal name of a subroutine might be different from
! !  the external name!)
!
! ! Subroutines used with all filters
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &            ! Initialize dimension of observation vector
       obs_op_pdaf, &                  ! Implementation of the Observation operator
       init_obs_pdaf, &                ! Routine to provide vector of measurements
       prepoststep_ens_pdaf, &         ! User supplied pre/poststep routine
       prodRinvA_pdaf, &               ! Provide product R^-1 A for some matrix A
       init_obsvar_pdaf, &             ! Initialize mean observation error variance
       next_observation_pdaf, &        ! Provide time step, model time, and dimension of next observation
       distribute_state_pdaf           ! Routine to distribute a state vector to model fields
! ! Subroutines for local filters
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
! ! Subroutines used in EnKF
  EXTERNAL :: add_obs_error_pdaf, &    ! Add obs. error covariance R to HPH in EnKF
       init_obscovar_pdaf              ! Initialize obs error covar R in EnKF
! ! Subroutines used for localization in LEnKF
  EXTERNAL :: localize_covar_pdaf       ! Apply localization to HP and HPH^T
! ! Subroutines used in NETF
  EXTERNAL :: likelihood_pdaf          ! Compute observation likelihood for an ensemble member
! ! Subroutines used in LNETF
  EXTERNAL :: likelihood_l_pdaf        ! Compute local observation likelihood for an ensemble member
! ! Subroutines used in LKNETF
  EXTERNAL :: likelihood_hyb_l_pdaf, & ! Compute local likelihood awith hybrid weight for an ensemble member
       prodRinvA_hyb_l_pdaf            ! Provide product R^-1 A for some matrix A including hybrid weight


! *********************************
! *** Call assimilation routine ***
! *********************************

  IF (filtertype == 1) THEN
     CALL PDAF_assimilate_seik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
          prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 2) THEN
     CALL PDAF_assimilate_enkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
          next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 3) THEN
     CALL PDAF_assimilate_lseik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 4) THEN
     CALL PDAF_assimilate_etkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 5) THEN
     CALL PDAF_assimilate_letkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 6) THEN
     CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 7) THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 8) THEN
     CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, localize_covar_pdaf, add_obs_error_pdaf, &
          init_obscovar_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 9) THEN
     CALL PDAF_assimilate_netf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, &
          obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
          likelihood_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 10) THEN
     CALL PDAF_assimilate_lnetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
          likelihood_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 11) THEN
     CALL PDAF_assimilate_lknetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf, &
          prodRinvA_l_pdaf, prodRinvA_hyb_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, &
          likelihood_l_pdaf, likelihood_hyb_l_pdaf, next_observation_pdaf, status_pdaf)
  ELSE IF (filtertype == 12) THEN
     CALL PDAF_assimilate_pf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, &
          obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
          likelihood_pdaf, next_observation_pdaf, status_pdaf)
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_assimilate - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
