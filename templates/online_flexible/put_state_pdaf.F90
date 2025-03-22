!>  Routine to call PDAF for analysis step in flexible parallelization
!!
!! This routine is used in the case of the flexible ensemble
!! parallelization variant when PDAF*_put_state routines are
!! used. It is called during the model integrations at the time
!! when an analysis step should be computed. It calls the
!! assimilation routine for the DA method (PDAF3_put_state_X)
!! which computes the analysis step inside PDAF.
!!
!! In this routine, the real names of most of the 
!! user-supplied routines for PDAF are specified (see below).
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code for OMI
!! * 2025-03 - Lars Nerger - Adaption for PDAF3
!! * Other revisions - see repository log
!!
SUBROUTINE put_state_pdaf()

  USE PDAF, &                     ! PDAF interface definitions
       ONLY: PDAF3_put_state_local, PDAF3_put_state_global, &
       PDAF_localfilter
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, abort_parallel

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to 
! PDAF_get_state and PDAF3_put_state_X. This allows the user 
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
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain


! *********************************
! *** Call put_state routine    ***
! *********************************

  ! Call put_state routine for global or local filter
  IF (PDAF_localfilter() == 1) THEN
     ! Call generic PDAF3 interface routine for domain-localized filters
     CALL PDAF3_put_state_local(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          prepoststep_pdaf, status_pdaf)
  ELSE
     ! Call generic PDAF3 interface routine for global filters
     CALL PDAF3_put_state_global(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          prepoststep_pdaf, status_pdaf)
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF3_put_state - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE put_state_pdaf
