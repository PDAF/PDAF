!>  Routine to call PDAF for analysis step in flexible parallelization
!!
!! This routine is used in the case of the flexible ensemble parallelization
!! variant. It is called during the model integrations at the time
!! when an analysis step should be computed. It calls the filter-specific
!! assimilation routine of PDAF (PDAFomi_put_state_X) which computes the
!! analysis step inside PDAF.
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code for OMI
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_put_state_local, PDAFomi_put_state_global, &
       PDAFomi_put_state_lenkf, PDAF_get_localfilter, PDAFomi_generate_obs
  USE PDAFlocal, &                ! Interface definitions for PDAFlocal
       ONLY: PDAFlocalomi_put_state
  USE mod_parallel_pdaf, &        ! Parallelization variables
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
       prepoststep_pdaf               ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &      ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi         ! Apply localization to covariance matrix in LEnKF


! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter == 1) THEN
     ! Call generic OMI interface routine for domain-localized filters
     CALL PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_pdafomi, &
          obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdafomi, status_pdaf)
  ELSE
     IF (filtertype == 8) THEN
        ! LEnKF has its own OMI interface routine
        CALL PDAFomi_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdafomi, &
             obs_op_pdafomi, prepoststep_pdaf, localize_covar_pdafomi, status_pdaf)
     ELSE IF (filtertype == 100) THEN
        ! Observation generation can only be implemented using the fully-parallel
        ! implementation variant - but it is run with a single ensemble member
     ELSE
        ! Call generic OMI interface routine for global filters
        CALL PDAFomi_put_state_global(collect_state_pdaf, init_dim_obs_pdafomi, &
             obs_op_pdafomi, prepoststep_pdaf, status_pdaf)
     END IF
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAFomi_put_state - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
