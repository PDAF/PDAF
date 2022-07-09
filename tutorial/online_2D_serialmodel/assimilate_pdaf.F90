!$Id$
!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-specific assimilation routine of PDAF 
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf(step)

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
       PDAFomi_assimilate_lenkf, PDAF_get_localfilter
  USE mod_parallel_pdaf, &        ! Parallelization
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: filtertype, async

  IMPLICIT NONE

! *** Arguments
  INTEGER, INTENT(in) :: step

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag
  INTEGER :: localfilter          ! Flag for domain-localized filter (1=true)
  INTEGER :: assim_stat           ! Flag whether analysis step was just computed
  INTEGER :: dim_obs_init         ! Dummy variable for number of observations

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
       prepoststep_ens_pdaf           ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf, &             ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &              ! Get state on local analysis domain from global state
       l2g_state_pdaf                 ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_async_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &      ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi         ! Apply localization to covariance matrix in LEnKF



  IF (async) CALL obs_op_async(step)

! *********************************
! *** Call assimilation routine ***
! *********************************

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter==1) THEN
     CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_async_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
          next_observation_pdaf, status_pdaf)
  ELSE
     IF (filtertype/=8) THEN
        ! All global filters, except LEnKF
        CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_async_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
             next_observation_pdaf, status_pdaf)
     ELSE
        ! localized EnKF has its own OMI interface routine
        CALL PDAFomi_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_async_pdafomi, obs_op_pdafomi, prepoststep_ens_pdaf, &
             localize_covar_pdafomi, next_observation_pdaf, status_pdaf)
     END IF
  END IF

  ! Asynchronous DA: After the analysis initialize next set of observations
  CALL PDAF_get_assim_flag(assim_stat)
  IF (async .AND. assim_stat==1) CALL init_dim_obs_pdafomi(step, dim_obs_init)


  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAFomi_assimilate - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf


SUBROUTINE obs_op_async(step)

  USE PDAF_interfaces_module, ONLY: PDAF_set_ens_pointer
  USE obs_A_pdafomi, ONLY: obs_times_A, thisobs, ostate_A
  USE mod_parallel_pdaf, ONLY: mype_world
  
  IMPLICIT NONE

! *** Arguments
  INTEGER, INTENT(in) :: step

! *** Local variables
  INTEGER :: cnt, i
  INTEGER :: status
  REAL, POINTER :: ens_pointer(:,:)

  cnt = 0
  DO i = 1, thisobs%dim_obs_p
     IF (step == obs_times_A(i)) THEN
        cnt = cnt+1
     END IF
  END DO
  if (cnt>0 .and. mype_world==0) write (*,*) 'Number of obs at step ', step, ' =', cnt

  IF (cnt>0) THEN
     ! In case of observations at this time

     CALL PDAF_set_ens_pointer(ens_pointer, status)
     CALL collect_state_pdaf(thisobs%dim_state_p, ens_pointer)
     CALL obs_op_async_pdafomi(step, dim_state_p, thisobs%dim_obs_p, ens_pointer, ostate_A)
  END IF


END SUBROUTINE obs_op_async
