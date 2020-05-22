!$Id$
!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-speific assimilation routine of PDAF 
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_assimilate_estkf, PDAF_assimilate_lestkf
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf       ! PDAF status flag


  ! External subroutines
  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &  ! Collect a state vector from model fields
       distribute_state_pdaf, &      ! Distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_ens_pdaf          ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_f_pdafomi, &     ! Get dimension of full obs. vector for PE-local domain
       obs_op_f_pdafomi, &           ! Obs. operator for full obs. vector for PE-local domain
       init_obs_f_pdafomi, &         ! Provide full vector of measurements for PE-local domain
       init_dim_obs_l_pdafomi, &     ! Get dimension of obs. vector for local analysis domain
       g2l_obs_pdafomi, &            ! Get local observation vector from global observation vector
       init_obs_l_pdafomi, &         ! Provide vector of observations for local analysis domain
       prodRinvA_l_pdafomi, &        ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdafomi, &      ! Initialize local mean observation error variance
       init_obsvar_pdafomi, &        ! Initialize mean observation error variance
       prodRinvA_pdafomi             ! Provide product R^-1 A for some matrix A for global filter


! *********************************
! *** Call assimilation routine ***
! *********************************

  IF (filtertype == 6) THEN
     CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdafomi, obs_op_f_pdafomi, init_obs_f_pdafomi, prepoststep_ens_pdaf, &
          prodRinvA_pdafomi, init_obsvar_pdafomi, next_observation_pdaf, status_pdaf)
  ELSEIF (filtertype == 7) THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdafomi, obs_op_f_pdafomi, init_obs_f_pdafomi, init_obs_l_pdafomi, &
          prepoststep_ens_pdaf, prodRinvA_l_pdafomi, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdafomi, init_obsvar_pdafomi, init_obsvar_l_pdafomi, &
          next_observation_pdaf, status_pdaf)
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
