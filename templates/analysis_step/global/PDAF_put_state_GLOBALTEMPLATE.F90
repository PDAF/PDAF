!> Interface to PDAF for GLOBALTEMPLATE DA method
!!
!! Interface routine called from the model or PDAF_assimilate after
!! the forecast of each ensemble state. The routine manages
!! the transfer of a model state from the model to PDAF.
!! For the parallelization this involves transfer from model
!! processes (processing elements, PEs) to filter processes.
!!
!! For the flexible parallelization variance the routine manages
!! the ensemble. Thus, during the forecast phase state vectors are 
!! re-initialized from the forecast model fields by U\_collect\_state. 
!! At the end of a forecast phase when all ensemble members have
!! been integrated by the model sub-ensembles are gathered from
!! the model tasks by the routine PDAF\_gather\_ens). Subsequently
!! the filter update is performed by PDAF\_GLOBALTEMPLATE\_update.
!!
!! The code is very generic. The filter-specific part is the call
!! to the update-routine PDAF\_GLOBALTEMPLATE\_update where the analysis
!! is computed. The filter-specific subroutines that are specified
!! as arguments in the call to PDAF\_put\_state\_GLOBALTEMPLATE are
!! passed though to the update routine PDAF\_GLOBALTEMPLATE\_update.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a new DA method, the only required changes to
!! this routine should be
!! - replace 'GLOBALTEMPLATE' by the name of the new method
!! - adapt the argument lists in PDAF\_put\_state\_GLOBALTEMPLATE
!!   and PDAF\_GLOBALTEMPLATE\_update
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_put_state_GLOBALTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, outflag)

  USE PDAF_timer, &                 ! Routines for timings
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &           ! Routine for memory counting
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &            ! Variables for framework functionality
       ONLY: screen, flag, initevol, offline_mode, &
       dim_p, dim_ens, local_dim_ens, nsteps, step_obs, step, &
       state, ens, Ainv, member, member_save, &
       subtype_filter, sens, dim_lag, cnt_maxlag
  USE PDAF_mod_filtermpi, &         ! Variables for parallelization
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model
  USE PDAF_communicate_ens, &       ! Ensemble gathering or scattering
       ONLY: PDAF_gather_ens
  USE PDAFobs, &                    ! Routines and variables for observations
       ONLY: dim_obs
  USE PDAF_GLOBALTEMPLATE_update, & ! Name of method-specific update routine
       ONLY: PDAFGLOBALTEMPLATE_update

  IMPLICIT NONE
  
! TEMPLATE: 'outflag' is standard and should be kept

! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! TEMPLATE: The external subroutines depends on the DA method and should be adapted

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  ! Routines for ensemble framework
  EXTERNAL :: U_collect_state, &   !< Write model fields into state vector
       U_prepoststep               !< User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &    !< Initialize dimension of observation vector
       U_obs_op, &                 !< Observation operator
       U_init_obs, &               !< Initialize observation vector
       U_init_obsvar, &            !< Initialize mean observation error variance
       U_prodRinvA                 !< Provide product R^-1 A

! TEMPLATE: The local variables are usually generic and don't need changes

! *** local variables ***
  INTEGER :: i                     ! Counter


! **************************************************
! *** Save forecasted state back to the ensemble ***
! *** Only done on the filter Pes                ***
! **************************************************

! TEMPLATE: This is generic as long as subtype_filter 2 and 3 are fixed ensemble cases
  doevol: IF (nsteps > 0 .OR. .NOT.offline_mode) THEN

     CALL PDAF_timeit(41, 'new')

     modelpes: IF (modelpe) THEN

        ! Store member index for PDAF_get_memberid
        member_save = member

! TEMPLATE: Only this IF-statement should be adapted if subtype_filter 2, 3 are used differently
        IF (subtype_filter /= 2 .AND. subtype_filter /= 3) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, ens(1 : dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1:dim_p))
        END IF
     END IF modelpes

     CALL PDAF_timeit(41, 'old')

! TEMPLATE: do NOT change the member counting as it will break PDAF's ensemble forecast handling

     member = member + 1
  ELSE
     member = local_dim_ens + 1
  END IF doevol

  IF (filter_no_model .AND. filterpe) THEN
     member = local_dim_ens + 1
  END IF


! ********************************************************
! *** When forecast phase is completed                 ***
! ***   - collect forecast sub_ensembles on filter PEs ***
! ***   - perform analysis step                        ***
! ***   - re-initialize forecast counters/flags        ***
! ********************************************************

! TEMPLATE: Everything below is generic except the call to PDAFGLOBALTEMPLATE_update

  completeforecast: IF (member == local_dim_ens + 1 &
       .OR. offline_mode) THEN

     ! ***********************************************
     ! *** Collect forecast ensemble on filter PEs ***
     ! ***********************************************

     doevolB: IF (nsteps > 0) THEN

        IF (.not.filterpe) THEN
           ! Non filter PEs only store a sub-ensemble
           CALL PDAF_gather_ens(dim_p, dim_ens_l, ens, screen)
        ELSE
           ! On filter PEs, the ensemble array has full size
           CALL PDAF_gather_ens(dim_p, dim_ens, ens, screen)
        END IF

     END IF doevolB

     ! *** call timer
     CALL PDAF_timeit(2, 'old')

     IF (.NOT.offline_mode .AND. mype_world == 0 .AND. screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of forecast phase:', PDAF_time_temp(2), 's'
     END IF


     ! **************************************
     ! *** Perform analysis on filter PEs ***
     ! **************************************

     ! Screen output
     IF (offline_mode .AND. mype_world == 0 .AND. screen > 0) THEN
        WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
        WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
        WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
     ENDIF

     OnFilterPE: IF (filterpe) THEN

! TEMPLATE: This needs to be adapted according use features of the DA method
!   Usually only the included call-back routines (U_*) are changed, but other
!   variables are kept unchanged
        CALL PDAFGLOBALTEMPLATE_update(step_obs, dim_p, dim_obs, dim_ens, &
             state, Ainv, ens, U_init_dim_obs, U_obs_op, &
             U_init_obs, U_prodRinvA, U_init_obsvar, U_prepoststep, &
             screen, subtype_filter, dim_lag, sens, cnt_maxlag, flag)

     END IF OnFilterPE


     ! ***********************************
     ! *** Set forecast counters/flags ***
     ! ***********************************
     initevol = 1
     member   = 1
     step     = step_obs + 1

  END IF completeforecast


! ********************
! *** finishing up ***
! ********************

  outflag = flag

END SUBROUTINE PDAF_put_state_GLOBALTEMPLATE
