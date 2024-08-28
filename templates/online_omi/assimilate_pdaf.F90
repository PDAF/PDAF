!>  Routine to call PDAF for analysis step in fully-parallel mode
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-specific assimilation routine of PDAF 
!! (PDAFomi_assimilate_X), which checks whether the forecast phase
!! is completed. If so, the analysis step is computed inside PDAF.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code for OMI
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE PDAF_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_global, &
       PDAFomi_assimilate_lenkf, PDAF_get_localfilter, PDAFomi_generate_obs
  USE PDAFlocal, &                ! Interface definitions for PDAFlocal
       ONLY: PDAFlocalomi_assimilate
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
       distribute_state_pdaf, &       ! Distribute a state vector to model fields
       next_observation_pdaf, &       ! Provide time step of next observation
       prepoststep_pdaf               ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &      ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi         ! Apply localization to covariance matrix in LEnKF
  ! Subroutine used for generating observations
  EXTERNAL :: get_obs_f_pdaf          ! Get vector of synthetic observations from PDAF


! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter == 1) THEN
     ! Call generic OMI interface routine for domain-localized filters
     CALL PDAFlocalomi_assimilate(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, next_observation_pdaf, status_pdaf)
  ELSE
     IF (filtertype == 8) THEN
        ! LEnKF has its own OMI interface routine
        CALL PDAFomi_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
             localize_covar_pdafomi, next_observation_pdaf, status_pdaf)
     ELSE IF (filtertype == 100) THEN
        ! Observation generation has its own OMI interface routine
        CALL PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_f_pdaf, &
             prepoststep_pdaf, next_observation_pdaf, status_pdaf)
     ELSE
        ! Call generic OMI interface routine for global filters
        CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
             next_observation_pdaf, status_pdaf)
     END IF
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAFomi_assimilate - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
