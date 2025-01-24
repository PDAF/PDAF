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
!$Id$
!BOP
!
! !ROUTINE: PDAF_lknetf_step_update --- Control analysis update of the 2-step LKNETF
!
! !INTERFACE:
SUBROUTINE  PDAF_lknetf_step_update(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Uinv, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prodRinvA_hyb_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, U_likelihood_l, U_likelihood_hyb_l, &
     U_prepoststep, screen, subtype, incremental, type_forget, &
     dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the 2-step LKNETF.
!
! The analysis is performed by first preparing several
! global quantities on the PE-local domain, like the
! observed part of the state ensemble for all local
! analysis domains on the PE-local state domain.
! Then the analysis is done in two steps local analysis 
! domains in the PE-local state domain.  Either the NETF is 
! computed before the LETKF (subtype=0) or LETKF is 
! computed before NETF (subtype=1). 
! In the local analysis loops, the local state and 
! observation dimensions are initialized and the global 
! state ensemble is restricted to the local analysis domain.
! In addition, the routine U\_prepoststep is called prior
! to the analysis and after the resampling outside of
! the loop over the local domains to allow the user
! to access the ensemble information.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2018-01 - Lars Nerger - Initial code based on PDAF_lknetf_update
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: obs_member, type_trans, filterstr, forget, inloop, &
       member_save, type_hyb, hyb_g, hyb_k, &
       skewness, kurtosis, store_rndmat, debug
  USE PDAFomi, &
       ONLY: omi_omit_obs => omit_obs
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l, npes_filter, COMM_filter, MPIerr
  USE PDAF_analysis_utils, &
       ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, &
       PDAF_incr_local_obsstats, PDAF_print_local_obsstats
  USE PDAFobs, &
       ONLY: PDAFobs_initialize, PDAFobs_dealloc, &
       HX_f => HX_p, HXbar_f => HXbar_p, obs_f => obs_p

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
! !   suffix _f: Denotes a full variable of all observations required for the
! !              analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_p         ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)         ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens, dim_ens) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)  ! PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)     ! PE-local state analysis increment
  INTEGER, INTENT(in) :: screen        ! Verbosity flag
  INTEGER, INTENT(in) :: subtype       ! Filter subtype
  INTEGER, INTENT(in) :: incremental   ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget   ! Type of forgetting factor
  INTEGER, INTENT(in) :: dim_lag       ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag       ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_obs_op, &    ! Observation operator
       U_init_n_domains_p, & ! Provide number of local analysis domains
       U_init_dim_l, &       ! Init state dimension for local ana. domain
       U_init_dim_obs, &     ! Initialize dimension of observation vector
       U_init_dim_obs_l, &   ! Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &         ! Initialize observation vector
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_obsvar, &      ! Initialize mean observation error variance
       U_init_obsvar_l, &    ! Initialize local mean observation error variance
       U_g2l_state, &        ! Get state on local ana. domain from global state
       U_l2g_state, &        ! Init full state from state on local analysis domain
       U_g2l_obs, &          ! Restrict full obs. vector to local analysis domain
       U_prodRinvA_hyb_l, &  ! Compute product of R^(-1) with HV with hybrid weight
       U_likelihood_l, &     ! Compute likelihood
       U_likelihood_hyb_l, & ! Compute likelihood with hybrid weight
       U_prepoststep         ! User supplied pre/poststep routine

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_lknetf
! Calls: U_prepoststep
! Calls: U_init_n_domains_p
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_init_dim_l
! Calls: U_init_dim_obs_l
! Calls: U_g2l_state
! Calls: U_l2g_state
! Calls: U_g2l_obs
! Calls: PDAF_set_forget
! Calls: PDAF_inflate_ens
! Calls: PDAF_generate_rndmat
! Calls: PDAF_lknetf_compute_gamma
! Calls: PDAF_lknetf_ana_lnetf
! Calls: PDAF_lknetf_ana_letkfT
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: MPI_reduce
!EOP

! *** local variables ***
  INTEGER :: i, member, row        ! Counters
  INTEGER :: domain_p              ! Counter for local analysis domain
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  REAL    :: invdimens             ! Inverse global ensemble size
  INTEGER :: minusStep             ! Time step counter
  INTEGER :: n_domains_p           ! number of PE-local analysis domains
  REAL    :: forget_ana_l          ! forgetting factor supplied to analysis routine
  REAL    :: forget_ana            ! Possibly globally adaptive forgetting factor
  LOGICAL :: storerndmat = .FALSE. ! Store and reuse random rotation matrix
  REAL, ALLOCATABLE :: rndmat(:,:) ! random rotation matrix for ensemble trans.
  REAL, SAVE, ALLOCATABLE :: rndmat_save(:,:) ! Stored rndmat
  REAL :: invforget                ! inverse forgetting factor
  ! Variables on local analysis domain
  INTEGER :: dim_l                 ! State dimension on local analysis domain
  INTEGER :: dim_obs_l             ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: HX_l(:,:)   ! HX for local analysis domain
  REAL, ALLOCATABLE :: HXbar_l(:)  ! local observed mean state
  REAL, ALLOCATABLE :: obs_l(:)    ! local observation vector
  REAL, ALLOCATABLE :: resid_l(:)  ! local residual
  REAL, ALLOCATABLE :: ens_l(:,:)  ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:)  ! Mean state on local analysis domain
  REAL, ALLOCATABLE :: stateinc_l(:)  ! State increment on local analysis domain
  REAL, ALLOCATABLE :: n_eff(:)    ! Effective sample size for each local domain (in alpha)
  REAL, ALLOCATABLE :: n_eff_b(:)  ! Effective sample size for each local domain (hybrid)
  LOGICAL, ALLOCATABLE :: MASK(:)  ! Mask for effective sample sizes > 0
  REAL :: max_n_eff_l, min_n_eff_l ! PE-local min/max. effective ensemble sizes
  REAL :: max_n_eff, min_n_eff     ! Global min/max. effective ensemble sizes
  REAL :: max_gamma_l, min_gamma_l ! PE-local min/max. hybrid weight
  REAL :: max_gamma, min_gamma     ! Global min/max. hybrid weight
  REAL :: sum_gamma_l, mean_gamma  ! Local alpha sum; global mean alpha
  REAL :: sum_n_eff_l, mean_n_eff  ! Local sum of N_eff; global mean N_eff
  REAL :: max_stats_l(2), min_stats_l(2)   ! PE-local min/max of skewness and kurtosis
  REAL :: max_stats(2), min_stats(2)  ! Global min/max of skewness and kurtosis
  REAL :: sum_stats_l(2)           ! PE-local sum of skewness and kurtosis for averaging
  REAL :: mean_stats(2)            ! Global average skewness and kurtosis
  REAL, ALLOCATABLE :: Uinv_l(:,:) ! thread-local matrix Uinv
  REAL :: state_inc_p_dummy        ! Dummy variable to avoid compiler warning
  REAL, ALLOCATABLE :: gamma(:)    ! Hybrid weight for state update
  INTEGER :: cnt_small_svals       ! Counter for small values
  INTEGER :: n_domains_with_obs_p  ! Domain-local number of local domains with observations
  INTEGER :: n_domains_with_obs    ! Global number of local domains with observations


! **********************
! *** Initialization ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- START'

  ! Initialize variable to prevent compiler warning
  state_inc_p_dummy = state_inc_p(1)


! ***********************************
! *** Compute mean forecast state ***
! ***********************************

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


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  CALL PDAF_timeit(5, 'new')
  minusStep = - step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Uinv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF ', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
  END IF


! *********************************************
! *** Apply covariance inflation to global  ***
! *** ensemble if prior inflation is chosen ***            
! *********************************************

  CALL PDAF_timeit(51, 'new')

  IF (type_forget==0 .OR. type_forget==1) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: PDAF_letkf_update', debug, &
          'Inflate forcast ensemble: type_forget, forget', type_forget, forget

     CALL PDAF_timeit(14, 'new')
     CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget)
     CALL PDAF_timeit(14, 'old')
  ENDIF

  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  CALL PDAFobs_initialize(step, dim_p, dim_ens, dim_obs_f, &
       state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       screen, debug, .true., .true., .true., .true.)


! **************************************
! *** Preparation for local analysis ***
! **************************************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(4, 'new')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_int(3) -not used-     '
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_int(4) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_int(5) type_forget', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_int(6) type_trans ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_int(7) type_hyb   ', type_hyb

     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_real(1) forget   ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_real(2) hyb_gamma', hyb_g
     WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update', debug, &
          'Configuration: param_real(3) hyb_kappa', hyb_k
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_n_domains'

  ! Query number of analysis domains for the local analysis
  ! in the PE-local domain
  CALL PDAF_timeit(42, 'new')
  CALL U_init_n_domains_p(step, n_domains_p)
  CALL PDAF_timeit(42, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  n_domains_p', n_domains_p

  ! Initialize effective sample sizes and MASK
  ALLOCATE(n_eff(n_domains_p))
  ALLOCATE(n_eff_b(n_domains_p))
  n_eff = 0.0
  n_eff_b = 0.0
  ALLOCATE(MASK(n_domains_p))

  ! Allocate arrays for skewness and kurtosis
  IF(ALLOCATED(skewness)) DEALLOCATE(skewness)
  IF(ALLOCATED(kurtosis)) DEALLOCATE(kurtosis)
  ALLOCATE(skewness(n_domains_p))
  ALLOCATE(kurtosis(n_domains_p))
  skewness = 0.0
  kurtosis = 0.0

  ! Allocate arrays for hybrid weights
  ALLOCATE(gamma(n_domains_p))
  gamma = 0.0
  
  IF (screen > 0) THEN
     IF (mype == 0) THEN
        IF (subtype == 0) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - 2-step LKNETF-HNK: NETF before LETKF'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - 2-step LKNETF-HKN: LETKF before NETF'
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
  IF (type_forget == 5) THEN
     CALL PDAF_set_forget(step, filterstr, dim_obs_f, dim_ens, HX_f, &
          HXbar_f, obs_f, U_init_obsvar, forget, forget_ana)
  ENDIF

  ! *** Initialize random transformation matrix
  CALL PDAF_timeit(33, 'new')
  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  rnd_store: IF (.NOT. storerndmat .OR. (storerndmat .AND. allocflag == 0)) THEN

     IF (type_trans == 0) THEN
        ! Initialize random matrix
        IF (screen > 0 .AND. mype == 0) &
             WRITE (*, '(a, 6x, a)') 'PDAF', '--- Initialize random transformation'
        CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  rndmat', rndmat
     ELSE
        IF (screen > 0 .AND. mype == 0) &
             WRITE (*, '(a, 6x, a)') 'PDAF', '--- Initialize deterministic transformation'
        rndmat = 0.0
        DO i = 1, dim_ens
           rndmat(i,i) = 1.0
        END DO
     END IF

     IF (store_rndmat) THEN
        IF (allocflag == 0) ALLOCATE(rndmat_save(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        rndmat_save = rndmat
     END IF

  ELSE rnd_store
     ! Re-use stored rndmat
     if (mype == 0 .AND. screen > 0) &
          write (*,'(a, 5x, a)') 'PDAF', '--- Use stored random rotation matrix'
     rndmat = rndmat_save

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  stored rndmat', rndmat
  END IF rnd_store

  CALL PDAF_timeit(33, 'old')
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(4, 'old')


! *************************************
! *** Perform analysis - first loop ***
! *************************************

  CALL PDAF_timeit(6, 'new')

  ! Initialize counters for statistics on local observations
  CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, obs_l, HX_l, HXbar_l, resid_l, ens_l) &
!$OMP private(state_l, stateinc_l, Uinv_l, forget_ana_l)

  forget_ana_l = forget_ana

  ! Allocate ensemble transform matrix
  ALLOCATE(Uinv_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  Uinv_l = 0.0

  ! initialize number of small singular values
  cnt_small_svals = 0

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- Enter first local analysis loop'

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysis: DO domain_p = 1, n_domains_p

     ! Set flag that we are in the local analysis loop
     inloop = .true.

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- First local analysis for domain_p', domain_p
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_dim_l'
     END IF

     ! local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  dim_l', dim_l
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_dim_obs_l'
     END IF

     ! Get observation dimension for local domain
     CALL PDAF_timeit(9, 'new')
     dim_obs_l = 0
     CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
     CALL PDAF_timeit(9, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  dim_obs_l', dim_obs_l

     CALL PDAF_timeit(51, 'new')

     ! Gather statistical information on local observations
     CALL PDAF_incr_local_obsstats(dim_obs_l)
     
     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     ALLOCATE(stateinc_l(dim_l))
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(15, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_memberid
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update -- call g2l_state for ensemble member', member
           if (member==1) &
                WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update --    Note: if ens_l is incorrect check user-defined indices in g2l_state!'
        END IF

        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  ens_l', ens_l(:,member)
     END DO

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_lknetf_update -- call g2l_state for ensemble mean'

     member_save = 0
     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  meanens_l', state_l

     CALL PDAF_timeit(15, 'old')

     IF (dim_obs_l > 0) THEN

        ALLOCATE(HX_l(dim_obs_l, dim_ens))
        ALLOCATE(HXbar_l(dim_obs_l))
        ALLOCATE(obs_l(dim_obs_l))

        ! Restrict full observation to local observation
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- call g2l_obs', dim_ens, 'times'

        CALL PDAF_timeit(46, 'new')
        DO member = 1, dim_ens
           obs_member = member
           CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:,member), HX_l(:,member))
        END DO

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call g2l_obs for mean'

        obs_member = 0
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)
        CALL PDAF_timeit(46, 'old')

        ! get local observation vector
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_obs_l'

        CALL PDAF_timeit(47, 'new')
        CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
        CALL PDAF_timeit(47, 'old')

        ! Omit observations with too high innovation
        IF (omi_omit_obs)  THEN
           CALL PDAF_timeit(51, 'new')

           ALLOCATE(resid_l(dim_obs_l))

           resid_l = obs_l - HXbar_l
           CALL PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, resid_l, obs_l)

           DEALLOCATE(resid_l)

           CALL PDAF_timeit(51, 'old')
        END IF

        CALL PDAF_timeit(7, 'new')


        ! *** 2-step LKNETF analysis - STEP 1 ***

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- Compute hybrid gamma'

        CALL PDAF_timeit(53,'new')
        CALL PDAF_lknetf_compute_gamma(domain_p, step, dim_obs_l, dim_ens, &
             HX_l, HXbar_l, obs_l, type_hyb, hyb_g, hyb_k, &
             gamma(domain_p), n_eff(domain_p), skewness(domain_p), kurtosis(domain_p), &
             U_likelihood_l, screen, flag)
        CALL PDAF_timeit(53,'old')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  gamma', gamma(domain_p)

        IF (subtype == 0) THEN

           ! 2-step LKNETF with NETF before LETKF 
           CALL PDAF_lknetf_ana_lnetf(domain_p, step, dim_l, dim_obs_l, &
                dim_ens, ens_l, HX_l, rndmat, obs_l, U_likelihood_hyb_l, &
                cnt_small_svals, n_eff_b(domain_p), gamma(domain_p), screen, flag)
        ELSE IF (subtype == 1) THEN

           ! 2-step LKNETF with LETKF before NETF
           CALL PDAF_lknetf_ana_letkfT(domain_p, step, dim_l, dim_obs_l, &
                dim_ens, state_l, Uinv_l, ens_l, HX_l, &
                HXbar_l, stateinc_l, rndmat, forget_ana_l, &
                obs_l, U_prodRinvA_hyb_l, U_init_obsvar_l, U_init_n_domains_p, &
                gamma(domain_p), screen, incremental, type_forget, flag)
        END IF

        DEALLOCATE(obs_l, HX_l, HXbar_l)

        CALL PDAF_timeit(7, 'old')

     ELSE

        ! UNOBSERVED DOMAIN

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- dim_obs_l = 0; omit call to local analysis function'

        IF (type_forget==1) THEN
           ! prior inflation NOT on unobserved domains - take it back!
           invforget = 1.0/forget
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_p, ens_l, invforget)
        ENDIF

     END IF

     CALL PDAF_timeit(16, 'new')
 
     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update -- call l2g_state for ensemble member', member
           WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  ens_l', ens_l(:,member)
        END IF

        member_save = member
        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
    
     ! Initialize global state increment
!      IF (incremental == 1) THEN
!         IF (debug>0) THEN
!            WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- init gobal state increment'
!            WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  stateinc_l', stateinc_l
!         END IF
!
!         member_save = member
!         CALL U_l2g_state(step, domain_p, dim_l, stateinc_l, dim_p, state_inc_p)
!      END IF

     CALL PDAF_timeit(16, 'old')

     CALL PDAF_timeit(51, 'new')
     CALL PDAF_timeit(17, 'new')

     ! *** Perform smoothing of past ensembles ***
!     CALL PDAF_smoother_local(domain_p, step, dim_p, dim_l, dim_ens, &
!          dim_lag, Uinv_l, ens_l, sens_p, cnt_maxlag, &
!          U_g2l_state, U_l2g_state, forget_ana, screen)

     CALL PDAF_timeit(17, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l, stateinc_l)
     CALL PDAF_timeit(51, 'old')

  END DO localanalysis

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- End of first local analysis loop'

  ! Set flag that we are not in the local analysis loop
  inloop = .false.


! ***************************************
! *** Prepare second step of analysis ***
! ***************************************

  CALL PDAF_timeit(4, 'new')

  ! HX = [Hx_1 ... Hx_(r+1)] for full DIM_OBS_F region on PE-local domain
  CALL PDAF_timeit(44, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call obs_op', dim_ens, 'times'

  ENSA: DO member = 1,dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! Call observation operator
     CALL U_obs_op(step, dim_p, dim_obs_f, ens_p(:, member), HX_f(:, member))
  END DO ENSA
  CALL PDAF_timeit(44, 'old')

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(11, 'new')

  ! *** Compute mean forecast state
  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')

  ! *** Compute mean state of ensemble on PE-local observation space 
  HXbar_f = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_obs_f
        HXbar_f(row) = HXbar_f(row) + invdimens * HX_f(row, member)
     END DO
  END DO
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(4, 'old')


! **************************************
! *** Perform analysis - second loop ***
! **************************************

  CALL PDAF_timeit(6, 'new')

  forget_ana_l = forget_ana

  ! Allocate ensemble transform matrix
  Uinv_l = 0.0

  ! initialize number of small singular values
  cnt_small_svals = 0

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- Enter second local analysis loop'

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysisA: DO domain_p = 1, n_domains_p

     ! Set flag that we are in the local analysis loop
     inloop = .true.

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- Second local analysis for domain_p', domain_p
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_dim_l'
     END IF

     ! local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  dim_l', dim_l
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_dim_obs_l'
     END IF

     ! Get observation dimension for local domain
     CALL PDAF_timeit(9, 'new')
     dim_obs_l = 0
     CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
     CALL PDAF_timeit(9, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  dim_obs_l', dim_obs_l

     CALL PDAF_timeit(51, 'new')
     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     ALLOCATE(stateinc_l(dim_l))
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(15, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update -- call g2l_state for ensemble member', member
           if (member==1) &
                WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update --    Note: if ens_l is incorrect check user-defined indices in g2l_state!'
        END IF

        member_save = member
        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  ens_l', ens_l(:,member)
     END DO

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_lknetf_update -- call g2l_state for ensemble mean'

     member_save = 0
     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  meanens_l', state_l

     CALL PDAF_timeit(15, 'old')

     IF (dim_obs_l > 0) THEN

        ! get local observation vector
        ALLOCATE(HX_l(dim_obs_l, dim_ens))
        ALLOCATE(HXbar_l(dim_obs_l))
        ALLOCATE(obs_l(dim_obs_l))

        ! Restrict full observation to local observation
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- call g2l_obs', dim_ens, 'times'

        CALL PDAF_timeit(46, 'new')
        DO member = 1, dim_ens
           obs_member = member
           CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:,member), HX_l(:,member))
        END DO

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call g2l_obs for mean'

        obs_member = 0
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)
        CALL PDAF_timeit(46, 'old')

        ! get local observation vector
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- call init_obs_l'

        CALL PDAF_timeit(47, 'new')
        CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
        CALL PDAF_timeit(47, 'new')

        CALL PDAF_timeit(7, 'new')


        ! *** 2-step LKNETF analysis - STEP 2 ***

        IF (subtype == 0) THEN

           ! 2-step LKNETF with NETF before LETKF 
           CALL PDAF_lknetf_ana_letkfT(domain_p, step, dim_l, dim_obs_l, &
                dim_ens, state_l, Uinv_l, ens_l, HX_l, &
                HXbar_l, stateinc_l, rndmat, forget_ana_l, &
                obs_l, U_prodRinvA_hyb_l, U_init_obsvar_l, U_init_n_domains_p, &
                gamma(domain_p), screen, incremental, type_forget, flag)
        ELSE IF (subtype == 1) THEN
           CALL PDAF_lknetf_ana_lnetf(domain_p, step, dim_l, dim_obs_l, &
                dim_ens, ens_l, HX_l, rndmat, obs_l, U_likelihood_hyb_l, &
                cnt_small_svals, n_eff_b(domain_p), gamma(domain_p), screen, flag)
        END IF

        IF (type_forget==3) THEN
           IF (debug>0) &
                WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- apply inflation to analysis ensemble'

           ! Apply forgetting factor to posterior ensemble - only observed domains
           CALL PDAF_timeit(14, 'new')
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_l, ens_l, forget)
           CALL PDAF_timeit(14, 'old')
        ENDIF

        CALL PDAF_timeit(7, 'old')

        DEALLOCATE(obs_l, HX_l, HXbar_l)

     ELSE
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_update -- dim_obs_l = 0; omit call to local analysis function'
     END IF

     CALL PDAF_timeit(16, 'new')
 
     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lknetf_update -- call l2g_state for ensemble member', member
           WRITE (*,*) '++ PDAF-debug PDAF_lknetf_update:', debug, '  ens_l', ens_l(:,member)
        END IF

        member_save = member
        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
    
     ! Initialize global state increment
!      IF (incremental == 1) THEN
!         member_save = member
!         CALL U_l2g_state(step, domain_p, dim_l, stateinc_l, dim_p, state_inc_p)
!      END IF

     CALL PDAF_timeit(16, 'old')

     CALL PDAF_timeit(51, 'new')
     CALL PDAF_timeit(17, 'new')

     ! *** Perform smoothing of past ensembles ***
!     CALL PDAF_smoother_local(domain_p, step, dim_p, dim_l, dim_ens, &
!          dim_lag, Uinv_l, ens_l, sens_p, cnt_maxlag, &
!          U_g2l_state, U_l2g_state, forget_ana, screen)

     CALL PDAF_timeit(17, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l, stateinc_l)
     CALL PDAF_timeit(51, 'old')

  END DO localanalysisA

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- End of second local analysis loop'

  ! Set flag that we are not in the local analysis loop
  inloop = .false.

  CALL PDAF_timeit(51, 'new')

  IF (type_forget==2) THEN
     ! Apply forgetting factor to posterior ensemble - all domains
     CALL PDAF_timeit(14, 'new')
     CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget)
     CALL PDAF_timeit(14, 'old')
  ENDIF

!$OMP CRITICAL
  ! Set Uinv - required for subtype=3
  Uinv = Uinv_l
!$OMP END CRITICAL

  DEALLOCATE(Uinv_l)
!$OMP END PARALLEL

  CALL PDAF_timeit(6, 'old')

  CALL PDAF_timeit(3, 'old')

  ! Initialize mask array for computing global effective ensemble size
  MASK = (n_eff_b > 0)

  ! *** Print statistics for local analysis to the screen ***
  CALL PDAF_print_local_obsstats(screen, n_domains_with_obs_p)

  IF (npes_filter>1) THEN
     ! Number globally observed domains
     CALL MPI_Reduce(n_domains_with_obs_p, n_domains_with_obs, 1, MPI_INTEGER, MPI_SUM, &
          0, COMM_filter, MPIerr)

     ! Min/max effective sample sizes
     max_n_eff_l = MAXVAL(n_eff_b)
     CALL MPI_Reduce(max_n_eff_l, max_n_eff, 1, MPI_REALTYPE, MPI_MAX, &
          0, COMM_filter, MPIerr)
     min_n_eff_l = MINVAL(n_eff_b, MASK)
     CALL MPI_Reduce(min_n_eff_l, min_n_eff, 1, MPI_REALTYPE, MPI_MIN, &
          0, COMM_filter, MPIerr)
     sum_n_eff_l = SUM(n_eff_b)
     CALL MPI_Reduce(sum_n_eff_l, mean_n_eff, 1, MPI_REALTYPE, MPI_SUM, &
          0, COMM_filter, MPIerr)
     mean_n_eff = mean_n_eff / REAL(n_domains_with_obs)

     ! Min/max hybrid weight
     max_gamma_l = MAXVAL(gamma)
     CALL MPI_Reduce(max_gamma_l, max_gamma, 1, MPI_REALTYPE, MPI_MAX, &
          0, COMM_filter, MPIerr)
     min_gamma_l = MINVAL(gamma, MASK)
     CALL MPI_Reduce(min_gamma_l, min_gamma, 1, MPI_REALTYPE, MPI_MIN, &
          0, COMM_filter, MPIerr)
     sum_gamma_l = SUM(gamma)
     CALL MPI_Reduce(sum_gamma_l, mean_gamma, 1, MPI_REALTYPE, MPI_SUM, &
          0, COMM_filter, MPIerr)
     mean_gamma = mean_gamma / REAL(n_domains_with_obs)

     ! Min/max skewness and kurtosis
     max_stats_l(1) = MAXVAL(skewness)/SQRT(REAL(dim_ens))
     max_stats_l(2) = MAXVAL(kurtosis)/REAL(dim_ens)
     CALL MPI_Reduce(max_stats_l, max_stats, 2, MPI_REALTYPE, MPI_MAX, &
          0, COMM_filter, MPIerr)

     min_stats_l(1) = MINVAL(skewness, MASK)/SQRT(REAL(dim_ens))
     min_stats_l(2) = MINVAL(kurtosis, MASK)/REAL(dim_ens)
     CALL MPI_Reduce(min_stats_l, min_stats, 2, MPI_REALTYPE, MPI_MIN, &
          0, COMM_filter, MPIerr)

     sum_stats_l(1) = SUM(skewness)/SQRT(REAL(dim_ens))
     sum_stats_l(2) = SUM(kurtosis)/REAL(dim_ens)
     CALL MPI_Reduce(sum_stats_l, mean_stats, 2, MPI_REALTYPE, MPI_SUM, &
          0, COMM_filter, MPIerr)
     mean_stats = mean_stats / REAL(n_domains_with_obs)

  ELSE
     ! Min/max effective ensemble sizes
     max_n_eff = MAXVAL(n_eff_b)
     min_n_eff = MINVAL(n_eff_b, MASK)
     mean_n_eff = SUM(n_eff_b) / n_domains_with_obs_p

     ! Min/max effective ensemble sizes
     max_gamma = MAXVAL(gamma)
     min_gamma = MINVAL(gamma, MASK)
     mean_gamma = SUM(gamma) / n_domains_with_obs_p

     ! Min/max skewness and kurtosis
     max_stats(1) = MAXVAL(skewness)/SQRT(REAL(dim_ens))
     max_stats(2) = MAXVAL(kurtosis)/REAL(dim_ens)
     min_stats(1) = MINVAL(skewness, MASK)/SQRT(REAL(dim_ens))
     min_stats(2) = MINVAL(kurtosis, MASK)/REAL(dim_ens)
     mean_stats(1) = SUM(skewness)/SQRT(REAL(dim_ens)) / REAL(n_domains_with_obs_p)
     mean_stats(2) = SUM(kurtosis)/REAL(dim_ens)/ REAL(n_domains_with_obs_p)
  END IF

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 8x, a, 23x, 3f10.3)') &
         'PDAF', 'Minimal/Maximal/Mean N_eff:', &
         min_n_eff, max_n_eff, mean_n_eff
     WRITE (*, '(a, 8x, a, 15x, 3f10.3)') &
         'PDAF', 'Minimal/Maximal/Mean hybrid weight:', &
         min_gamma, max_gamma, mean_gamma
     WRITE (*, '(a, 8x, a, 1x, 3f10.3)') &
         'PDAF', 'Minimal/Maximal/Mean abs. skewness/SQRT(dim_ens):', &
         min_stats(1), max_stats(1), mean_stats(1)
     WRITE (*, '(a, 8x, a, 7x, 3f10.3)') &
         'PDAF', 'Minimal/Maximal/Mean abs. kurtosis/dim_ens:', &
         min_stats(2), max_stats(2), mean_stats(2)

     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF
  CALL PDAF_timeit(51, 'old')

! *** Clean up from local analysis update ***
  DEALLOCATE(HX_f, HXbar_f, rndmat)
#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

! *** Poststep for analysis ensemble ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Uinv, ens_p, flag)
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

  DEALLOCATE(n_eff, n_eff_b, MASK, gamma)
  DEALLOCATE(skewness, kurtosis)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_update -- END'

END SUBROUTINE PDAF_lknetf_step_update
