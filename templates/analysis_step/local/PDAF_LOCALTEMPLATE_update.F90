! Copyright (c) 2004-2024 Lars Nerger
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
!
!>Control analysis update of LOCALTEMPLATE
!!
!! This routine prepares the actual analysis update which
!! is computed the in PDAF_analysis_LOCALTEMPLATE.
!!
!! The analysis is performed by first preparing several
!! global quantities on the PE-local domain, like the
!! observed part of the state ensemble for all local
!! analysis domains on the PE-local state domain.
!! Then the analysis (PDAF\_LOCALTEMPLATE\_analysis) is performed
!! within a loop over all local analysis domains in the PE-local 
!! state domain. In this loop, the local state and 
!! observation dimensions are initialized and the global 
!! state ensemble is restricted to the local analysis domain.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling outside of
!! the loop over the local domains to allow the user
!! to access the ensemble information.
!!
!! ADAPTING THE TEMPLATE
!! The structure of the operations included in this template is
!! typical for a (domain)-local ensemble DA method. One should
!! check if other operations are required in the global preparations
!! and the local steps which are included in this template. Further,
!! the call to PDAF_LOCALTEMPLATE_analysis and the particular
!! arguments of this subroutine need to be adapted. For the
!! different operations see the particular comments in the code.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on LETKF
!! * Later revisions - see repository log
!!
SUBROUTINE  PDAF_LOCALTEMPLATE_update(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, &
     U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_prepoststep, screen, subtype, incremental, &
     dim_lag, sens_p, cnt_maxlag, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_filter, &
       ONLY: forget, type_trans, filterstr, obs_member, &
       inloop, member_save, debug
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_analysis_utils, &
       ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, PDAF_incr_local_obsstats, &
       PDAF_print_local_obsstats

  IMPLICIT NONE

! *** Arguments ***
! Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
!    suffix _f: Denotes a full variable of all observations required for the
!               analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_p         ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)         ! PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens) ! Transform matrix
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)  ! PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen        ! Verbosity flag
  INTEGER, INTENT(in) :: subtype       ! Filter subtype
  INTEGER, INTENT(in) :: incremental   ! Control incremental updating
  INTEGER, INTENT(in) :: dim_lag       ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag       ! Status flag

! ** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
  ! Routine for ensemble framework
  EXTERNAL :: U_prepoststep         ! User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &     ! Initialize dimension of observation vector
       U_obs_op, &                  ! Observation operator
       U_init_dim_obs_l, &          ! Initialize dim. of obs. vector for local ana. domain
       U_init_obs_l, &              ! Init. observation vector on local analysis domain
       U_g2l_obs, &                 ! Restrict full obs. vector to local analysis domain
       U_prodRinvA_l                ! Provide product R^-1 A on local analysis domain
  ! Routines for state localization
  EXTERNAL :: U_init_n_domains_p, & ! Provide number of local analysis domains
       U_init_dim_l, &              ! Init state dimension for local ana. domain
       U_g2l_state, &               ! Get state on local ana. domain from full state
       U_l2g_state                  ! Init full state from state on local analysis domain

! *** Local variables ***
  INTEGER :: i, j, member, row      ! Counters
  INTEGER :: minusStep              ! Time step counter
  INTEGER :: domain_p               ! Counter for local analysis domain
  INTEGER :: n_domains_p            ! number of PE-local analysis domains
  REAL    :: invdimens              ! Inverse global ensemble size
  INTEGER, SAVE :: allocflag = 0    ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HX_f(:,:)    ! HX for PE-local ensemble
  REAL, ALLOCATABLE :: HXmean_f(:)  ! PE-local observed mean state
  REAL, ALLOCATABLE :: obs_f(:)     ! PE-local observation vector
  REAL, ALLOCATABLE :: rndmat(:,:)  ! random rotation matrix for ensemble trans.
  ! Variables on local analysis domain
  INTEGER :: dim_l                  ! State dimension on local analysis domain
  INTEGER :: dim_obs_l              ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: ens_l(:,:)   ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:)   ! Mean state on local analysis domain
  REAL, ALLOCATABLE :: stateinc_l(:) ! State increment on local analysis domain
  REAL, ALLOCATABLE :: Ainv_l(:,:)  ! thread-local matrix Ainv


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

! +++ TEMPLATE:
! +++ For fixed-ensemble cases (like Ensemble OI) only the
! +++ central state is integrated by the model. Here, we
! +++ then need to add the ensemble perturbations

  CALL PDAF_timeit(51, 'new')

  fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis

  CALL PDAF_timeit(51, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the forecast
! +++ ensemble is standard and should be kept

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) &
          WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF


! *************************************************
! *** Non-local preparations for local analysis ***
! *************************************************

! +++ TEMPLATE: 
! +++ Before entering the local analysis loop we here initialize
! +++ - the number of local analysis domains
! +++ - the ensemble mean state
! +++ - the (parallel) full dimension of observations
! +++ - the observed ensemble and the mean of it
! +++ - if requested a random rotation matrix

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(4, 'new')

  ! *** Query number of analysis domains for the local analysis
  ! *** in the process-local domain
  CALL PDAF_timeit(42, 'new')
  CALL U_init_n_domains_p(step, n_domains_p)
  CALL PDAF_timeit(42, 'old')
  
  IF (screen > 0) THEN
     IF (mype == 0) THEN
! +++ TEMPLATE
! +++ If other subtypes exist, add them here
        IF (subtype == 0 .OR. subtype == 2) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - LOCALTEMPLATE default analysis'
        END IF
     END IF
     IF (screen<3) THEN
        CALL PDAF_print_domain_stats(n_domains_p)
     ELSE
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
     END IF
  END IF


! *** Local analysis: initialize non-local quantities ***

  ! *** Compute mean forecast state

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')
  CALL PDAF_timeit(51, 'old')

  ! *** Get observation dimension for all observations required 
  ! *** for the loop of local analyses on the PE-local domain.

  CALL PDAF_timeit(43, 'new')
  CALL U_init_dim_obs(step, dim_obs_f)
  CALL PDAF_timeit(43, 'old')

  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a, i6, a, i10)') &
          'PDAF', '--- PE-Domain:', mype, ' dimension of process-local full obs. vector', dim_obs_f
  END IF

  ! *** Apply observation operator to all ensemble states
  ! *** for full DIM_OBS_F region on PE-local domain

  ALLOCATE(HX_f(dim_obs_f, dim_ens))
  ALLOCATE(HXmean_f(dim_obs_f))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f * dim_ens + dim_obs_f)

  CALL PDAF_timeit(44, 'new')
  ENS: DO member = 1,dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! Call observation operator
     CALL U_obs_op(step, dim_p, dim_obs_f, ens_p(:, member), HX_f(:, member))
  END DO ENS
  CALL PDAF_timeit(44, 'old')


  ! *** Compute mean state of ensemble on PE-local observation space 
  CALL PDAF_timeit(51, 'new')

  HXmean_f = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_obs_f
        HXmean_f(row) = HXmean_f(row) + invdimens * HX_f(row, member)
     END DO
  END DO

  ! *** Initialize random transformation matrix

  CALL PDAF_timeit(33, 'new')
  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

! +++ TEMPLATE: 
! +++ One might want to apply a random rotation to Ainv.
! +++ PDAF_general_rndmat provides a matrix with such random rotation
! +++ which is used in this example. For the local filter, the random
! +++ matrix is initialized here and then consistently used for all
! +++ local analysis domains

  IF (type_trans == 2) THEN
     ! Initialize random matrix
     CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)
  ELSE
     ! Initialize identity matrix
     rndmat = 0.0
     DO i = 1, dim_ens
        rndmat(1,1) = 1.0
     END DO
  END IF

  CALL PDAF_timeit(33, 'old')
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(4, 'old')


! ******************************
! *** Perform local analysis ***
! ******************************

  CALL PDAF_timeit(6, 'new')

  ! Initialize counters for statistics on local observations
  CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, stateinc_l, Ainv_l, flag)

  ! Allocate ensemble transform matrix
  ALLOCATE(Ainv_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  Ainv_l = 0.0

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysis: DO domain_p = 1, n_domains_p

! +++ TEMPLATE: 
! +++ Inside the local analysis domain the typical operations are
! +++ - Determine the local state dimension
! +++ - Determine the local number of observations 
! +++ - Initialize the local state ensemble and mean state
! +++ - Call the analysis routine of the local DA method (obviously method-specific)
! +++ - Initialize the global state ensmelbe and mean state from the local analysis states
! +++ - Potentially apply smoothing

     ! Set flag that we are in the local analysis loop
     inloop = .true.

     ! *** Get local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     ! *** Get observation dimension for local domain
     CALL PDAF_timeit(9, 'new')
     dim_obs_l = 0
     CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
     CALL PDAF_timeit(9, 'old')

     CALL PDAF_timeit(51, 'new')

     ! Gather statistical information on local observations
     CALL PDAF_incr_local_obsstats(dim_obs_l)
     
     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     ALLOCATE(stateinc_l(dim_l))

     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(15, 'new')

     ! *** Intialize state ensemble and mean state on local analysis domain
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_memberid
        member_save = member

        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))
     END DO

     ! Store member index to make it accessible with PDAF_get_memberid
     member_save = 0

     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     CALL PDAF_timeit(15, 'old')

     CALL PDAF_timeit(7, 'new')

! +++ TEMPLATE:
! +++ The call to PDAF_LOCALTEMPLATE_analysis should be adapted for
! +++ the call-back routines and variables used by the method

     ! *** Analysis step
     CALL PDAF_LOCALTEMPLATE_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
          dim_ens, state_l, Ainv_l, ens_l, HX_f, &
          HXmean_f, rndmat, forget, &
          U_g2l_obs, U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, &
          screen, flag)

     CALL PDAF_timeit(7, 'old')

     CALL PDAF_timeit(16, 'new')
 
     ! *** re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        member_save = member

        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
     IF (subtype == 3) THEN
        ! Initialize global state for ETKF with fixed covariance matrix
        member_save = 0

        CALL U_l2g_state(step, domain_p, dim_l, state_l, dim_p, state_p)
     END IF

     CALL PDAF_timeit(16, 'old')

     CALL PDAF_timeit(51, 'new')
     CALL PDAF_timeit(17, 'new')

! +++ TEMPLATE 
! +++ The smoother routine is generic. It can be used
! +++ as long as sens_l(:,lag) * (forget*Ainv + diag(invdimens))
! +++ yields the smoother update for the ensemble at lag 'lag'. 
! +++ The routine includes the projection onto the local
! +++ state vector and back.

     ! *** Perform smoothing of past ensembles ***
     CALL PDAF_smoother_local(domain_p, step, dim_p, dim_l, dim_ens, &
          dim_lag, Ainv_l, ens_l, sens_p, cnt_maxlag, &
          U_g2l_state, U_l2g_state, forget, screen)

     CALL PDAF_timeit(17, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l, stateinc_l)
     CALL PDAF_timeit(51, 'old')

  END DO localanalysis

  ! Set flag that we are not in the local analysis loop
  inloop = .false.

  CALL PDAF_timeit(51, 'new')

!$OMP CRITICAL
  ! Set Ainv - required for subtype=3
  Ainv = Ainv_l
!$OMP END CRITICAL

  DEALLOCATE(Ainv_l)
!$OMP END PARALLEL


  ! *** Print statistics for local analysis to the screen ***
  CALL PDAF_print_local_obsstats(screen)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(6, 'old')
  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF

! *** Clean up from local analysis update ***
  DEALLOCATE(HX_f, HXmean_f, rndmat)
#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif


! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the analysis
! +++ ensemble is standard and should be kept

  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
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

END SUBROUTINE PDAF_LOCALTEMPLATE_update
