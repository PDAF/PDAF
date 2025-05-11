!>  Routine to call PDAF for analysis step
!!
!! This routine performs a single analysis step in the
!! offline implementation of PDAF. For this, it calls the
!! filter-specific assimilation routine of PDAF. For the
!! offline implementation this is PDAF3_assim_offline.
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
       ONLY: PDAF3_assim_offline
  USE mod_parallel_pdaf, &        ! Parallelization
       ONLY: mype_world, abort_parallel

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the call to 
! PDAF3_assim_offline. This allows the user to specify
! the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different
! from the external name!

  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: prepoststep_ens_offline ! User supplied pre/poststep routine
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

  ! Call universal PDAF3 assim_offline routine
  CALL PDAF3_assim_offline(init_dim_obs_pdafomi, obs_op_pdafomi, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
       prepoststep_ens_offline, status_pdaf)

! +++ Note: The universal routine PDAF3_assim_offline can be used to
! +++ execute all filter methods. The specified routines for localization
! +++ are only executed if a local filter is used. If one uses
! +++ exclusively global filters or the LEnKF, one can use the specific
! +++ routine PDAF3_assim_offline_global which does not include the
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
