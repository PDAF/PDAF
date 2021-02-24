!$Id$
!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-speific assimilation routine of PDAF 
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
       PDAF_get_localfilter
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, abort_parallel, task_id, mype_model, mype_submodel
  USE mod_assim_pdaf, &      ! Variables for assimilation
       ONLY: filtertype, step_null, restart, dim_state_p
  USE mod_assim_atm_pdaf, ONLY: dp
  USE mo_time_control,  ONLY: get_time_step

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf             ! PDAF status flag
  INTEGER :: localfilter             ! Flag for domain-localized filter (1=true)
  INTEGER :: seed_id
  REAL(dp), ALLOCATABLE :: sta_p(:)
  INTEGER :: istep


  ! External subroutines
  EXTERNAL :: collect_state_pdaf, &  ! Routine to collect a state vector from model fields
       distribute_state_pdaf, &      ! Routine to distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_pdaf              ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain


! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

  ! Store time step
  istep = get_time_step()

  IF (istep==step_null .AND. .NOT. restart) THEN

     ! Here we initialize the model state for ECHAM from the model fields of ECHAM
     ! we can add some randomness

     IF (mype_submodel==0 .and. task_id==1) WRITE(*,*) 'assmilate_pdaf: generate initial ensemble'

     ALLOCATE (sta_p(dim_state_p))

     CALL collect_state_pdaf(dim_state_p,sta_p)

!     seed_id=1000 * (task_id+1)+1
  
!     CALL distribute_state_ini_pdaf(dim_state_p, sta_p+ran(seed_id)) 

     CALL distribute_state_ini_pdaf(dim_state_p, sta_p) 
  
     DEALLOCATE(sta_p)
  END IF


! *********************************
! *** Call assimilation routine ***
! *********************************

!   if (mype_submodel==0 .and. task_id==1) write (*,'(2x,a,i2,a,i,a,i)') &
!        'ECHAM ',task_id,' ',mype_model,' assimilate_pdaf, step', istep-step_null

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter==1) THEN
     CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
          next_observation_pdaf, status_pdaf)
  ELSE
     ! All global filters except LEnKF
     CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
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
