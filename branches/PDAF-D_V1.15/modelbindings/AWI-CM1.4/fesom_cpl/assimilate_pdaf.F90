!$Id: assimilate_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf(istep)

! !DESCRIPTION:
! This routine is called during the model integrations at each time 
! step. It check whether the forecast phase is completed. If so, 
! PDAF_put_state_X is called to perform the analysis step.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, abort_parallel, task_id, mype_submodel
  USE mod_assim_pdaf, &        ! Variables for assimilation
       ONLY: filtertype, istep_asml

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: step
! CAlls: PDAF_assimilate_X
!EOP

  INTEGER, INTENT(in) :: istep

! Local variables
  INTEGER :: status_pdaf       ! PDAF status flag


  ! External subroutines
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector from model fields
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       next_observation_pdaf, &        ! Provide time step, model time, &
       init_obsvar_pdaf, &             ! Initialize mean observation error variance
       prepoststep_pdaf                ! User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain


! *********************************
! *** Call assimilation routine ***
! *********************************

  if(mype_submodel==0 .and. task_id==1) write (*,'(a,i1.1,a,i,a,i)') &
       'FESOM ',task_id,' ',mype_submodel,' assimilate_pdaf, step', istep

  istep_asml = istep

  IF (filtertype == 7) THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
          prepoststep_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
  ELSE
     ! No valid filter selected
     status_pdaf = 100
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
