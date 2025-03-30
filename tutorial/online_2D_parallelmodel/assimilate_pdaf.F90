!>  Routine to call PDAF for analysis step in fully-parallel mode
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-specific assimilation routine of PDAF 
!! (PDAFomi_assimilate_X), which checks whether the forecast phase
!! is completed. If so, the analysis step is computed inside PDAF.
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE PDAF, &                     ! PDAF interface definitions
       ONLY: PDAF3_assimilate
  USE mod_parallel_model, &       ! Parallelization variables
       ONLY: mype_world, abort_parallel

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag


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
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain


! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Call universal PDAF3 assimilation routine
  CALL PDAF3_assimilate(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       prepoststep_pdaf, next_observation_pdaf, status_pdaf)


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
