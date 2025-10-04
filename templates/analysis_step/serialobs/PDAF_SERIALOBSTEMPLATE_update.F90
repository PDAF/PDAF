!> Control analysis update of the GLOBALTEMPLATE
!!
!! Routine to control the analysis update of the GLOBALTEMPLATE.
!! 
!! The analysis and ensemble tranformation are performed by
!! calling PDAF\_GLOBALTEMPLATE\_analysis. In addition, the routine
!! U\_prepoststep is called prior to the analysis and after
!! the resampling to allow the user to access the ensemble
!! information.
!!
!! __Revision history:__
!! * 2025-10 - Lars Nerger - Initial template code based on ENSRF
!! * Later revisions - see repository log
!!
MODULE PDAF_SERIALOBSTEMPLATE_update

CONTAINS
SUBROUTINE PDAFSERIALOBSTEMPLATE_update(step, dim_p, dim_obs_p, dim_ens, state_p, &
     ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
     U_init_obsvars, U_localize_covar_serial, U_prepoststep, screen, &
     subtype, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAF_SERIALOBSTEMPLATE, &
       ONLY: debug, forget
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_init_obsvars, PDAFobs_dealloc, type_obs_init, &
       observe_ens, HX_p, HXbar_p, obs_p, var_obs_p
  USE PDAF_analysis_utils, &
       ONLY: PDAF_inflate_ens
  USE PDAF_SERIALOBSTEMPLATE_analysis, &
       ONLY: PDAF_SERIALOBSTEMPLATE_ana
  USE PDAFomi_obs_f, &
       ONLY: omi_n_obstypes => n_obstypes, omi_obs_diag => obs_diag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_p      !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens    !< Size of state ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
  INTEGER, INTENT(in) :: screen     !< Verbosity flag
  INTEGER, INTENT(in) :: subtype    !< Specification of filter subtype
  INTEGER, INTENT(inout) :: flag    !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &     !< Initialize dimension of observation vector
       U_obs_op, &                  !< Observation operator
       U_init_obs, &                !< Initialize observation vector
       U_init_obsvars, &            !< Initialize vector of observation error variances
       U_localize_covar_serial, &   !< Apply localization for single-observation vectors
       U_prepoststep                !< User supplied pre/poststep routine

! *** local variables ***
  INTEGER :: i                      ! Counter
  INTEGER :: minusStep              ! Time step counter
  INTEGER, SAVE :: allocflag = 0    ! Flag whether first time allocation is done
  LOGICAL :: do_init_dim_obs        ! Flag for initializing dim_obs_p in PDAFobs_init
  LOGICAL :: do_ensmean             ! Flag for computing ensemble mean state
  REAL :: Ainv(1, 1)                ! Unused array, but required in call to U_prepoststep


! **********************
! ***  Update phase  ***
! **********************


! ************************
! *** Inflate ensemble ***
! ************************

! +++ TEMPLATE: The observations can be initialized here before
! +++ the call to U_prepoststep, or afterwards (see below)
! +++ The observation arrays (obs_p, HX_p, HXbar_p) are declared
! +++ in the module PDAFobs. The routine PDAFobs_initialize
! +++ allocates and fills theses arrays. In addition it can 
! +++ compute the ensemble mean state in state_p. The five logical 
! +++ options at the end of he argument list define which steps
! +++ are executed in the subroutine (see subroutine for its
! +++ documentation).

  CALL PDAF_timeit(3, 'new')

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! We need to call the inflation of the forecast ensemble before
     ! the observed ensemble is initialized in PDAFobs_init

     CALL PDAF_timeit(51, 'new')

     do_ensmean = .true.
     IF (forget /= 1.0) THEN
        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(a, 5x, a, f10.3)') 'PDAF', 'Inflate forecast ensemble, forget=', forget
        
        ! Apply forgetting factor - this also computes the ensemble mean state_p
        CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

        ! PDAF_inflate_ens compute the ensmeble mean; thus don't do this in PDAFobs_init
        do_ensmean = .false.
     END IF

     CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! ***    optionally before call to U_prepoststep    ***
! *****************************************************

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, .true., .true., .true., .true.)

     ! Initialize vector var_obs_p in the module PDAFobs
     CALL PDAFobs_init_obsvars(step, dim_obs_p, U_init_obsvars)

  END IF

  CALL PDAF_timeit(3, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the forecast
! +++ ensemble is standard and should be kept

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
  END IF


! ************************
! *** Inflate ensemble ***
! ************************

! +++ The observation can be initialized here after the
! +++ call to U_prepoststep

  do_ensmean = .true.
  IF (type_obs_init==1) THEN
     ! If the observations are initnialize after prepoststep, we
     ! need to call the inflation of the forecast ensemble here.

     CALL PDAF_timeit(51, 'new')

     IF (forget /= 1.0) THEN
        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(a, 5x, a, f10.3)') 'PDAF', 'Inflate forecast ensemble, forget=', forget
        
        ! Apply forgetting factor - this also computes the ensemble mean state_p
        CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

        ! PDAF_inflate_ens compute the ensemble mean; thus don't do this in PDAFobs_init
        do_ensmean = .false.
     END IF

     CALL PDAF_timeit(51, 'old')
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! ***     optionally after call to U_prepoststep    ***
! *****************************************************

  IF (type_obs_init>0) THEN
     CALL PDAF_timeit(3, 'new')

     IF (type_obs_init==1) THEN
        do_init_dim_obs=.true.
     ELSE
        ! Skip call to U_init_dim_obs when also called before prepoststep
        do_init_dim_obs=.false.   
     END IF

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, do_init_dim_obs, .true., .true., .true.)

     ! Initialize vector var_obs_p in the module PDAFobs
     CALL PDAFobs_init_obsvars(step, dim_obs_p, U_init_obsvars)

     CALL PDAF_timeit(3, 'old')
  END IF


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE

  CALL PDAF_timeit(3, 'new')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_SERIALOBSTEMPLATE_update', debug, &
          'Configuration: param_int(3) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_SERIALOBSTEMPLATE_update', debug, &
          'Configuration: param_int(8) observe_ens ', observe_ens

     WRITE (*,*) '++ PDAF-debug PDAF_SERIALOBSTEMPLATE_update', debug, &
          'Configuration: param_real(1) forget     ', forget
  END IF

! ***  Execute Analysis step  ***

! +++ TEMPLATE:
! +++ The call to PDAF_GLOBALTEMPLATE_analysis should be adapted for
! +++ the call-back routines and variables used by the method

  CALL PDAF_SERIALOBSTEMPLATE_ana(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, HX_p, HXbar_p, obs_p, var_obs_p, &
       U_localize_covar_serial, screen, debug, flag)

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF


! ******************************************************
! *** Initialize analysis observed ensemble and mean ***
! ******************************************************

! +++ TEMPLATE:
! +++ This additional call to PDAFomi_init initializes
! +++ the observed analysis ensemble and its mean. This
! +++ is used in the observation diagnostics of PDAF-OMI

  IF (omi_n_obstypes>0 .AND. omi_obs_diag>0) THEN
     ! This call initializes HX_p, HXbar_p in the module PDAFobs
     ! for the analysis ensemble
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif


! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the analysis
! +++ ensemble is standard and should be kept

! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF


! ********************
! *** Finishing up ***
! ********************

! +++ TEMPLATE:
! +++ This call to PDAFobs_dealloc is mandatory; do not change it

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAFSERIALOBSTEMPLATE_update

END MODULE PDAF_SERIALOBSTEMPLATE_update
