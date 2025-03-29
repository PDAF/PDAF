!>  Routine to call PDAF for analysis step in flexible parallelization
!!
!! This routine is used in the case of the flexible ensemble
!! parallelization variant. It is called at each time step during
!! the model integrations. It calls the filter-specific assimilation
!! routine of PDAF (PDAF3_assimilate_X) which check whether the 
!! forecast needs to be continued (more steps or another ensemble
!! state). When the forecast phase is complete, the analysis step
!! is computed inside PDAF.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code for PDAF3 using PDAF3_assimilate
!! * Other revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE PDAF, &                     ! PDAF interface definitions
       ONLY: PDAF3_assimilate
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, abort_parallel

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to 
! PDAF_get_state and PDAF_assimilate_X. This allows the user 
! to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different
! from the external name!

  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: distribute_state_pdaf, &  ! Distribute a state vector to model fields
       collect_state_pdaf, &            ! Collect a state vector from model fields
       prepoststep_pdaf, &              ! User supplied pre/poststep routine
       next_observation_pdaf            ! Provide time step of next observation
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &    ! Provide number of local analysis domains
       init_dim_l_pdaf                  ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, &   ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi           ! Get dimension of obs. vector for local analysis domain


! *********************************
! *** Call assimilation routine ***
! *********************************

! +++ Note: The universal routine PDAF3_assimilate can be used to
! +++ execute all filter methods. The specified routines for localization
! +++ are only executed if a local filter is used. If one uses
! +++ exclusively global filters or the LEnKF, one can use the specific
! +++ routine PDAF3_assimilate_global which does not include the
! +++ arguments for localization. This would avoid to include routines
! +++ that are never called for global filters. 

  ! Call universal PDAF3 interface routine
  CALL PDAF3_assimilate(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          prepoststep_pdaf, next_observation_pdaf, status_pdaf)


! *************************
! *** Check status flag ***
! *************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF3_assimilate - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
