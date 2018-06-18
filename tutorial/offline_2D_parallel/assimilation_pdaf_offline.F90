!$Id: assimilation_pdaf_offline.F90 1366 2013-04-24 16:25:05Z lnerger $
!BOP
!
! !ROUTINE: assimilation_pdaf_offline - Control PDAF offline analysis
!
! !INTERFACE:
SUBROUTINE assimilation_pdaf_offline()

! !DESCRIPTION:
! This routine performs a single analysis step for
! PDAF in offline mode using PDAF with domain-decomposition.
!
! The analysis is performed by calling a filter-specific 
! routine PDAF\_put\_state\_X.
!
! In this routine, the real names of most of the 
! user-supplied routines for PDAF are specified (see below).
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code by restructuring
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &    ! Parallelization
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, & ! airables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
! ! External subroutines 
! !  (subroutine names are passed over to PDAF in the calls to 
! !  PDAF_get_state and PDAF_put_state_X. This allows the user 
! !  to specify the actual name of a routine. However, the 
! !  PDAF-internal name of a subroutine might be different from
! !  the external name!)
!
! ! Subroutines used with all filters
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &            ! Initialize dimension of observation vector
       obs_op_pdaf, &                  ! Implementation of the Observation operator
       init_obs_pdaf                   ! Routine to provide vector of measurements
! ! Subroutine used in SEIK and ETKF
  EXTERNAL :: prepoststep_ens_offline, & ! User supplied pre/poststep routine
       init_obsvar_pdaf                ! Initialize mean observation error variance
! ! Subroutines used in EnKF
  EXTERNAL :: add_obs_error_pdaf, &    ! Add obs. error covariance R to HPH in EnKF
       init_obscovar_pdaf              ! Initialize obs error covar R in EnKF
! ! Subroutine used in SEIK and ETKF
  EXTERNAL :: prodRinvA_pdaf           ! Provide product R^-1 A for some matrix A
! ! Subroutines used in LSEIK and LETKF
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
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state (possible, but not required!)
! Calls: PDAF_put_state_seik
! Calls: PDAF_put_state_enkf
! Calls: PDAF_put_state_lseik
! Calls: PDAF_put_state_etkf
! Calls: PDAF_put_state_letkf
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: status    ! Status flag for filter routines


! ************************
! *** Perform analysis ***
! ************************

! *** Note on PDAF_get_state for offline implementation: ***
! *** For the offline mode of PDAF the call to           ***
! *** PDAF_get_state is not required as no forecasting   ***
! *** is performed in this mode. However, it is save     ***
! *** to call PDAF_get_state, even it is not necessary.  ***
! *** The functionality of PDAF_get_state is deactived   ***
! *** for the offline mode.                              ***

  IF (filtertype == 1) THEN
     CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 2) THEN
     CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, add_obs_error_pdaf, init_obscovar_pdaf, &
          status)
  ELSE IF (filtertype == 3) THEN
     CALL PDAF_put_state_lseik( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  ELSE IF (filtertype == 4) THEN
     CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 5) THEN
     CALL PDAF_put_state_letkf( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  ELSE IF (filtertype == 6) THEN
     CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_offline, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 7) THEN
     CALL PDAF_put_state_lestkf( &
          collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_offline, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, status)
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilation_pdaf_offline
