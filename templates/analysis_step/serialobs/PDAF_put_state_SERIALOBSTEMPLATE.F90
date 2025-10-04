!> Interface to PDAF for SERIALOBSTEMPLATE DA method
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
!! the filter update is performed by PDAF\_SERIALOBSTEMPLATE\_update.
!!
!! The code is very generic. The filter-specific part is the call
!! to the update-routine PDAF\_SERIALOBSTEMPLATE\_update where the analysis
!! is computed. The filter-specific subroutines that are specified
!! as arguments in the call to PDAF\_put\_state\_SERIALOBSTEMPLATE are
!! passed though to the update routine PDAF\_SERIALOBSTEMPLATE\_update.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a new DA method, the only required changes to
!! this routine should be
!! - replace 'SERIALOBSTEMPLATE' by the name of the new method
!! - adapt the argument lists in PDAF\_put\_state\_SERIALOBSTEMPLATE
!!   and PDAF\_SERIALOBSTEMPLATE\_update
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
MODULE PDAFput_state_SERIALOBSTEMPLATE

CONTAINS

  SUBROUTINE PDAF_put_state_SERIALOBSTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op,  &
       U_init_obs, U_init_obsvars, U_localize_covar_serial, U_prepoststep, outflag)

    USE PDAF_timer, &                 ! Routines for timings
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_memcounting, &           ! Routine for memory counting
         ONLY: PDAF_memcount
    USE PDAF_mod_core, &              ! Variables for framework functionality
         ONLY: screen, flag, initevol, offline_mode, assim_flag, &
         dim_p, dim_ens, local_dim_ens, nsteps, step_obs, step, &
         state, ens, member, member_save, subtype_filter
    USE PDAF_mod_parallel, &         ! Variables for parallelization
         ONLY: mype_world, filterpe, &
         dim_ens_l, modelpe, filter_no_model
    USE PDAF_communicate_ens, &       ! Ensemble gathering or scattering
         ONLY: PDAF_gather_ens
    USE PDAFobs, &                    ! Routines and variables for observations
         ONLY: dim_obs
    USE PDAF_iau, &                   ! Step counter for incremental updating
         ONLY: step_cnt_iau
    USE PDAF_SERIALOBSTEMPLATE_update, &  ! Name of method-specific update routine
         ONLY: PDAFSERIALOBSTEMPLATE_update

    IMPLICIT NONE
  
! TEMPLATE: 'outflag' is standard and should be kept

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag   !< Status flag
  
! TEMPLATE: The external subroutines depends on the DA method and should be adapted

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    ! Routines for ensemble framework - generic and always needed
    EXTERNAL :: U_collect_state, &    !< Write model fields into state vector
         U_prepoststep                !< User supplied pre/poststep routine
    ! Observation-related routines for analysis step - generic and always needed
    EXTERNAL :: U_init_dim_obs, &     !< Initialize dimension of observation vector
         U_obs_op, &                  !< Observation operator
         U_init_obs                   !< Initialize observation vector
    ! Observation-related routines for analysis step - specific for the DA method
    EXTERNAL :: U_init_obsvars, &     !< Initialize vector of observation error variances
       U_localize_covar_serial        !< Apply localization for single-observation vectors

! TEMPLATE: The local variables are usually generic and don't need changes

! *** local variables ***
    INTEGER :: i, j                   ! Counters


! ***************************************************************
! *** Store forecasted state back to the ensemble array ens   ***
! *** and increment counter `member` for ensemble state index ***
! *** Only done on the filter processes                       ***
! ***************************************************************

! TEMPLATE: PDAF uses subtype_filter 10 and 11 for EnOI modes in which only the
!   state is integrated, but not the full ensemble. This is generic as long as
!   subtype_filter 10 and 11 are used in this way.
    doevol: IF (nsteps > 0 .OR. .NOT.offline_mode) THEN

       CALL PDAF_timeit(41, 'new')

       modelpes: IF (modelpe) THEN

          ! Store member index for PDAF_get_memberid
          member_save = member

          IF (subtype_filter /= 10 .AND. subtype_filter /= 11) THEN
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

       ! Reset step counter for IAU
       step_cnt_iau = 0
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

! TEMPLATE: Everything below is generic except the call to PDAFSERIALOBSTEMPLATE_update

    completeforecast: IF (member == local_dim_ens + 1 &
         .OR. offline_mode) THEN

       ! Set flag for assimilation
       assim_flag = 1


       ! ***********************************************
       ! *** Collect forecast ensemble on filter PEs ***
       ! ***********************************************

       doevolB: IF (nsteps > 0) THEN

          IF (.not.filterpe) THEN
             ! Non filter PEs only store a sub-ensemble
             CALL PDAF_gather_ens(dim_p, dim_ens_l, ens, state, screen)
          ELSE
             ! On filter PEs, the ensemble array has full size
             CALL PDAF_gather_ens(dim_p, dim_ens, ens, state, screen)
          END IF

       END IF doevolB


       ! **********************************************************
       ! *** For EnOI mode: add state to ensemble perturbations ***
       ! **********************************************************

       EnOI_mode: IF (subtype_filter==10 .OR. subtype_filter==11) THEN
          DO j = 1, dim_ens
             DO i = 1, dim_p
                ens(i, j) = ens(i, j) + state(i)
             END DO
          END DO
       END IF EnOI_mode

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
       ELSE IF (.NOT.offline_mode .AND. mype_world==0 .AND. screen > 0) THEN
          WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'
       ENDIF

       OnFilterPE: IF (filterpe) THEN

! TEMPLATE: This needs to be adapted according use features of the DA method
!   Usually only the included call-back routines (U_*) are changed, but other
!   variables are kept unchanged
          CALL  PDAFSERIALOBSTEMPLATE_update(step_obs, dim_p, dim_obs, dim_ens, state, &
               ens, U_init_dim_obs, U_obs_op, U_init_obs, &
               U_init_obsvars, U_localize_covar_serial, U_prepoststep, screen, &
               subtype_filter, flag)
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

  END SUBROUTINE PDAF_put_state_SERIALOBSTEMPLATE

END MODULE PDAFput_state_SERIALOBSTEMPLATE
