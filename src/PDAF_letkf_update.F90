! Copyright (c) 2004-2025 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!> Control analysis update of the LETKF
!!
!! Routine to control the analysis update of the LETKF filter.
!!
!! The analysis is performed by first preparing several
!! global quantities on the PE-local domain, like the
!! observed part of the state ensemble for all local
!! analysis domains on the PE-local state domain.
!! Then the analysis and ensemble tranformations are 
!! performed within a loop over all local analysis domains
!! in the PE-local state domain in the subroutine
!! (PDAF\_letkf\_analysis). In this loop, the local state and 
!! observation dimensions are initialized and the global 
!! state ensemble is restricted to the local analysis domain.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling outside of
!! the loop over the local domains to allow the user
!! to access the ensemble information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_letkf_update

CONTAINS
SUBROUTINE  PDAFletkf_update(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Ainv, ens_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prodRinvA_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, U_prepoststep, screen, &
     subtype, dim_lag, sens_p, cnt_maxlag, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_letkf, &
       ONLY: localfilter, debug, forget, type_forget, &
       type_trans, inloop, forget_l, member_save
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAF_analysis_utils, &
       ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, &
       PDAF_incr_local_obsstats, PDAF_print_local_obsstats, PDAF_set_forget, &
       PDAF_set_forget_local
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_init_local, PDAFobs_dealloc, PDAFobs_dealloc_local, &
       type_obs_init, observe_ens, HX_f => HX_p, HXbar_f => HXbar_p, obs_f => obs_p, &
       HX_l, HXbar_l, obs_l
  USE PDAF_analysis_utils, &
       ONLY: PDAF_generate_rndmat
  USE PDAF_smoother, &
       ONLY: PDAF_smoothing_local
  USE PDAFomi_obs_f, &
       ONLY: omi_n_obstypes => n_obstypes, omi_obs_diag => obs_diag
  USE PDAF_letkf_analysis, &
       ONLY: PDAF_letkf_ana
  USE PDAF_letkf_analysis_T, &
       ONLY: PDAF_letkf_ana_T
  USE PDAF_letkf_analysis_fixed, &
       ONLY: PDAF_letkf_ana_fixed

  IMPLICIT NONE

! *** Arguments ***
! Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
!    suffix _f: Denotes a full variable of all observations required for the
!               analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f  !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens)  !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  INTEGER, INTENT(in) :: subtype     !< Filter subtype
  INTEGER, INTENT(in) :: dim_lag     !< Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag     !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prepoststep         !< User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &     !< Initialize dimension of observation vector
       U_obs_op, &                  !< Observation operator
       U_init_dim_obs_l, &          !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                !< Initialize observation vector
       U_init_obs_l, &              !< Init. observation vector on local analysis domain
       U_g2l_obs, &                 !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l                !< Provide product R^-1 A on local analysis domain
  ! Routines for state localization
  EXTERNAL :: U_init_n_domains_p, & !< Provide number of local analysis domains
       U_init_dim_l, &              !< Init state dimension for local ana. domain
       U_init_obsvar, &             !< Initialize mean observation error variance
       U_init_obsvar_l, &           !< Initialize local mean observation error variance
       U_g2l_state, &               !< Get state on local ana. domain from full state
       U_l2g_state                  !< Init full state from state on local analysis domain

! *** local variables ***
  INTEGER :: i, j, member            ! Counters
  INTEGER :: domain_p                ! Counter for local analysis domain
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER :: minusStep               ! Time step counter
  INTEGER :: n_domains_p             ! number of PE-local analysis domains
  REAL    :: forget_ana_l            ! forgetting factor supplied to analysis routine
  REAL    :: forget_ana              ! Possibly globally adaptive forgetting factor
  LOGICAL :: storerndmat = .FALSE.   ! Store and reuse random rotation matrix
  LOGICAL :: do_init_dim_obs         ! Flag for initializing dim_obs_p in PDAFobs_init
  REAL, ALLOCATABLE :: rndmat(:,:)   ! random rotation matrix for ensemble trans.
  REAL, SAVE, ALLOCATABLE :: rndmat_save(:,:) ! Stored rndmat
  ! Variables on local analysis domain
  INTEGER :: dim_l                   ! State dimension on local analysis domain
  INTEGER :: dim_obs_l               ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: ens_l(:,:)    ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:)    ! Mean state on local analysis domain
  REAL, ALLOCATABLE :: Ainv_l(:,:)   ! thread-local matrix Ainv
 

! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- START'

  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(51, 'new')

  fixed_basis: IF (subtype == 10 .OR. subtype == 11) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis

  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .true., .true., .true., .true.)
  END IF
  CALL PDAF_timeit(3, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  CALL PDAF_timeit(5, 'new')
  minusStep = - step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
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
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., do_init_dim_obs, .true., .true., .true.)

     CALL PDAF_timeit(3, 'old')
  END IF


! **************************************
! *** Preparation for local analysis ***
! **************************************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(7, 'new')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(3) dim_lag     ', dim_lag
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(4) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(5) type_forget ', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(6) type_trans  ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(7) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_int(8) observe_ens           ', observe_ens

     WRITE (*,*) '++ PDAF-debug PDAF_letkf_update', debug, &
          'Configuration: param_real(1) forget     ', forget
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_n_domains'

  ! Query number of analysis domains for the local analysis
  ! in the PE-local domain
  CALL PDAF_timeit(42, 'new')
  CALL U_init_n_domains_p(step, n_domains_p)
  CALL PDAF_timeit(42, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  n_domains_p', n_domains_p
  
  IF (screen > 0) THEN
     IF (mype == 0) THEN
        IF (subtype == 0 .OR. subtype == 10) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - LETKF analysis using T-matrix'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - LETKF following Hunt et al. (2007)'
        ELSE IF (subtype == 10 .OR. subtype == 11) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'LETKF analysis for fixed covariance matrix'
        END IF
     END IF
     IF (screen<3) THEN
        CALL PDAF_print_domain_stats(n_domains_p)
     ELSE
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
     END IF
  END IF


! *** Local analysis: initialize global quantities ***

  CALL PDAF_timeit(51, 'new')

  ! *** Set forgetting factor globally
  forget_ana = forget
  IF (type_forget == 1) THEN
     CALL PDAF_set_forget(step, localfilter, dim_obs_f, dim_ens, HX_f, &
          HXbar_f, obs_f, U_init_obsvar, forget, forget_ana, &
          screen)
  ELSE IF (type_forget == 0) THEN
     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a, F7.2)') &
             'PDAF', '--- apply multiplicative inflation with fixed forget', forget
     END IF
  ENDIF

  ! *** Initialize random transformation matrix
  CALL PDAF_timeit(33, 'new')
  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  rnd_store: IF (.NOT. storerndmat .OR. (storerndmat .AND. allocflag == 0)) THEN

     IF (type_trans == 2) THEN
        ! Initialize random matrix
        CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  rndmat', rndmat

     ELSE
        rndmat = 0.0
     END IF

     IF (storerndmat) THEN
        ALLOCATE(rndmat_save(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        rndmat_save = rndmat
     END IF

  ELSE rnd_store
     ! Re-use stored rndmat
     if (mype == 0 .AND. screen > 0) &
          write (*,'(a, 5x, a)') 'PDAF', '--- Use stored random rotation matrix'
     rndmat = rndmat_save

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  rndmat', rndmat

  END IF rnd_store

  CALL PDAF_timeit(33, 'old')
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(7, 'old')


! ************************
! *** Perform analysis ***
! ************************

  CALL PDAF_timeit(8, 'new')

  ! Initialize counters for statistics on local observations
  CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, Ainv_l, flag, forget_ana_l)

  forget_ana_l = forget_ana

  ! Allocate ensemble transform matrix
  ALLOCATE(Ainv_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  Ainv_l = 0.0

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- Enter local analysis loop'

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysis: DO domain_p = 1, n_domains_p

     ! Set flag that we are in the local analysis loop
     inloop = .true.

     ! Set forgetting factor to global standard value
     forget_l = forget_ana


     ! *************************************
     ! *** Initialize local state vector ***
     ! *************************************

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_letkf_update -- local analysis for domain_p', domain_p
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_dim_l'
     END IF

     ! local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  dim_l', dim_l
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_dim_obs_l'
     END IF

     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))

     CALL PDAF_timeit(10, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_memberid
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update -- call g2l_state for ensemble member', member
           if (member==1) &
                WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update --    Note: if ens_l is incorrect check user-defined indices in g2l_state!'
        END IF

        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  ens_l', ens_l(:,member)
     END DO

     ! Store member index to make it accessible with PDAF_get_memberid
     member_save = 0

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_letkf_update -- call g2l_state for ensemble mean'

     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  meanens_l', state_l

     CALL PDAF_timeit(10, 'old')


     ! *******************************************
     ! *** Initialize local observation arrays ***
     ! *******************************************

     CALL PDAF_timeit(11, 'new')
     CALL PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f, dim_ens, &
          U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, debug)
     CALL PDAF_timeit(11, 'old')


     ! ************************************************
     ! *** Compute local adaptive forgetting factor ***
     ! ************************************************

     ! Reset forget (can be reset with PDAF_reset_forget)
     forget_ana_l = forget_l

     IF (type_forget == 2 .AND. dim_obs_l > 0) THEN
        CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, &
             HX_l, HXbar_l, obs_l, U_init_obsvar_l, forget, forget_ana_l)
     ENDIF


     ! *********************
     ! *** Analysis step ***
     ! *********************

     ! Gather statistical information on local observations
     CALL PDAF_incr_local_obsstats(dim_obs_l)
     
     CALL PDAF_timeit(12, 'new')

     IF (subtype == 0 .OR. subtype == 10) THEN
        ! *** LETKF analysis using T-matrix ***
        CALL PDAF_letkf_ana_T(domain_p, step, dim_l, dim_obs_l, dim_ens, &
             state_l, Ainv_l, ens_l, HX_l, HXbar_l, obs_l, &
             rndmat, forget_ana_l, U_prodRinvA_l, type_trans, screen, debug, flag)
     ELSE IF (subtype == 1) THEN
        ! *** ETKF analysis following Hunt et al. (2007) ***
        CALL PDAF_letkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, &
             state_l, Ainv_l, ens_l, HX_l, HXbar_l, obs_l, &
             rndmat, forget_ana_l, U_prodRinvA_l, type_trans, screen, debug, flag)
     ELSE IF (subtype == 11) THEN
        ! Analysis with state update but no ensemble transformation
        CALL PDAF_letkf_ana_fixed(domain_p, step, dim_l, dim_obs_l, dim_ens, &
             state_l, Ainv_l, ens_l, HX_l, HXbar_l, obs_l, &
             forget_ana_l, U_prodRinvA_l, screen, debug, flag)
     END IF

     CALL PDAF_timeit(12, 'old')
     CALL PDAF_timeit(14, 'new')

     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update -- call l2g_state for ensemble member', member
           WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  ens_l', ens_l(:,member)
        END IF

        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
     IF (subtype == 11) THEN
        ! Initialize global state for ETKF with fixed covariance matrix
        member_save = 0

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update -- call l2g_state for ensemble mean'
           WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  meanens_l', state_l
        END IF

        CALL U_l2g_state(step, domain_p, dim_l, state_l, dim_p, state_p)
     END IF

     CALL PDAF_timeit(14, 'old')
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_timeit(15, 'new')

     ! *** Perform smoothing of past ensembles ***
     CALL PDAF_smoothing_local(domain_p, step, dim_p, dim_l, dim_ens, &
          dim_lag, Ainv_l, ens_l, sens_p, cnt_maxlag, &
          U_g2l_state, U_l2g_state, forget_ana, screen)

     CALL PDAF_timeit(15, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l)
     CALL PDAFobs_dealloc_local()

     CALL PDAF_timeit(51, 'old')

  END DO localanalysis

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- End of local analysis loop'

  ! Set flag that we are not in the local analysis loop
  inloop = .false.

!$OMP CRITICAL
  ! Set Ainv - required for subtype=11
  Ainv = Ainv_l
!$OMP END CRITICAL

  DEALLOCATE(Ainv_l)
!$OMP END PARALLEL

  CALL PDAF_timeit(51, 'new')

  ! *** Print statistics for local analysis to the screen ***
  CALL PDAF_print_local_obsstats(screen)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(8, 'old')
  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF

! *** Clean up from local analysis update ***
  DEALLOCATE(rndmat)


! ******************************************************
! *** Initialize analysis observed ensemble and mean ***
! ******************************************************

  IF (omi_n_obstypes>0 .AND. omi_obs_diag>0) THEN
     ! This call initializes HX_p, HXbar_p in the module PDAFobs
     ! for the analysis ensemble
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif


! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
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

  IF (allocflag == 0) allocflag = 1

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- END'

END SUBROUTINE PDAFletkf_update

END MODULE PDAF_letkf_update
