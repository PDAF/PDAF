!>  Routine to call PDAF for analysis step
!!
!! This routine performs a single analysis step in the
!! offline implementation of PDAF. For this, it calls the
!! filter-specific assimilation routine of PDAF. For the
!! offline implementation this is PDAF_put_state_X.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! __Revision history:__
!! * 2009-11 - Lars Nerger - Initial code by restructuring
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf_offline()

  USE PDAF_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_put_state_global, &
       PDAFomi_put_state_lenkf, PDAF_get_localfilter
  USE PDAFlocal, &                ! Interface definitions for PDAFlocal
       ONLY: PDAFlocalomi_put_state
  USE mod_parallel_pdaf, &        ! Parallelization
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag
  INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to 
! PDAF_get_state and PDAF_put_state_X. This allows the user 
! to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different
! from the external name!

  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &   ! Collect a state vector from model fields
       distribute_state_pdaf, &       ! Distribute a state vector to model fields
       next_observation_pdaf, &       ! Provide time step of next observation
       prepoststep_ens_offline        ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &      ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi         ! Apply localization to covariance matrix in LEnKF


! *****************************
! *** Perform analysis step ***
! *****************************

! *** Note on PDAF_get_state for offline implementation: ***
! *** For the offline mode of PDAF the call to           ***
! *** PDAF_get_state is not required as no forecasting   ***
! *** is performed in this mode. However, it is save     ***
! *** to call PDAF_get_state, even it is not necessary.  ***
! *** The functionality of PDAF_get_state is deactivated ***
! *** for the offline mode.                              ***

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter == 1) THEN
     ! Call generic OMI interface routine for domain-localized filters
     CALL PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_pdafomi, &
          obs_op_pdafomi, prepoststep_ens_offline, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdafomi, status_pdaf)
  ELSE
     IF (filtertype /= 8) THEN
        ! Call generic OMI interface routine for global filters
        CALL PDAFomi_put_state_global(collect_state_pdaf, init_dim_obs_pdafomi, &
             obs_op_pdafomi, prepoststep_ens_offline, status_pdaf)
     ELSE
        ! LEnKF has its own OMI interface routine
        CALL PDAFomi_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdafomi, &
             obs_op_pdafomi, prepoststep_ens_offline, localize_covar_pdafomi, status_pdaf)
     END IF
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf_offline
