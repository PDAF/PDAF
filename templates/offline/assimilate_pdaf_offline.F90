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

  USE PDAF, &                     ! PDAF interface definitions
       ONLY: PDAF3_put_state
  USE mod_parallel_pdaf, &        ! Parallelization
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
       prepoststep_ens_offline        ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain


! *****************************
! *** Perform analysis step ***
! *****************************

! +++ Note on PDAF_get_state for offline implementation:
! +++ For the offline mode of PDAF the call to
! +++ PDAF_get_state is not required as no forecasting
! +++ is performed in this mode. However, it is save
! +++ to call PDAF_get_state, even it is not necessary.
! +++ The functionality of PDAF_get_state is deactivated
! +++ for the offline mode.

  ! Call universal PDAF3 put_state routine
  CALL PDAF3_put_state(collect_state_pdaf, &
       init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       prepoststep_ens_offline, status_pdaf)

! +++ Note: The universal routine PDAF3_put_state can be used to
! +++ execute all filter methods. The specified routines for localization
! +++ are only executed if a local filter is used. If one uses
! +++ exclusively global filters or the LEnKF, one can use the specific
! +++ routine PDAF3_put_state_global which does not include the
! +++ arguments for localization. This would avoid to include routines
! +++ that are never called for global filters. 


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
