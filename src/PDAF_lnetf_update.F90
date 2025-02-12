! Copyright (c) 2014-2025 Paul Kirchgessner
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
!> Control analysis update of the LNETF
!!
!! Routine to control the analysis update of the LNETF.
!!
!! The analysis is performed by first preparing several
!! global quantities on the PE-local domain, like the
!! observed part of the state ensemble for all local
!! analysis domains on the PE-local state domain.
!! Then the analysis (PDAF\_lnetf\_analysis) is performed within
!! a loop over all local analysis domains in the PE-local 
!! state domain. In this loop, the local state and 
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
!! * 2014-05 - Paul Kirchgessner - Initial code based on LETKF
!! * Later revisions - see repository log
!!
MODULE PDAF_lnetf_update

CONTAINS
SUBROUTINE  PDAFlnetf_update(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Ainv, ens_p, &
     U_obs_op, U_init_dim_obs, U_init_obs, U_init_obs_l, U_likelihood_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, &
     U_l2g_state, U_g2l_obs, U_prepoststep, screen, subtype, &
     dim_lag, sens_p, cnt_maxlag, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_lnetf, &
       ONLY: type_trans, type_noise, noise_amp, type_winf, limit_winf, &
       forget, type_forget, inloop, member_save, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l, npes_filter, COMM_filter, MPIerr
  USE PDAF_analysis_utils, &
       ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, &
       PDAF_incr_local_obsstats, PDAF_print_local_obsstats, &
       PDAF_add_particle_noise, PDAF_inflate_ens
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_init_local, PDAFobs_dealloc, &
       PDAFobs_dealloc_local, type_obs_init, HX_f => HX_p, &
       HX_l, obs_l
  USE PDAF_lnetf_analysis, &
       ONLY: PDAF_lnetf_ana, PDAF_lnetf_smootherT, PDAF_smoother_lnetf

  IMPLICIT NONE

! *** Arguments ***
! Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
!    suffix _f: Denotes a full variable of all observations required for the
!               analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f    !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)         !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)  !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: subtype       !< Filter subtype
  INTEGER, INTENT(inout) :: dim_lag    !< Status flag
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_obs_op, &          !< Observation operator
       U_init_n_domains_p, &       !< Provide number of local analysis domains
       U_init_dim_l, &             !< Init state dimension for local ana. domain
       U_init_dim_obs, &           !< Initialize dimension of observation vector
       U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &               !< Initialize PE-local observation vector
       U_init_obs_l, &             !< Init. observation vector on local analysis domain
       U_g2l_state, &              !< Get state on local ana. domain from global state
       U_l2g_state, &              !< Init full state from state on local analysis domain
       U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
       U_likelihood_l, &           !< Compute observation likelihood for an ensemble member
       U_prepoststep               !< User supplied pre/poststep routine

! *** local variables ***
  INTEGER :: i, j, member          ! Counters
  INTEGER :: domain_p              ! Counter for local analysis domain
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  INTEGER, SAVE :: allocflag_l = 0 ! Flag whether first time allocation is done
  INTEGER :: minusStep             ! Time step counter
  INTEGER :: n_domains_p           ! number of PE-local analysis domains
  LOGICAL :: do_init_dim_obs       ! Flag for initializing dim_obs_p in PDAFobs_init
  LOGICAL :: do_ensmean            ! Flag for computing ensemble mean state
  REAL, ALLOCATABLE :: TA_l(:,:)   ! Local ensemble transform matrix
  REAL, ALLOCATABLE :: HX_noinfl_f(:,:) ! HX for smoother (without inflation)
  REAL, ALLOCATABLE :: TA_noinfl_l(:,:) ! TA for smoother (without inflation)
  REAL, ALLOCATABLE :: rndmat(:,:) ! random rotation matrix for ensemble trans.
  ! Variables on local analysis domain
  INTEGER :: dim_l                 ! State dimension on local analysis domain
  INTEGER :: dim_obs_l             ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: ens_l(:,:)  ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:)  ! Mean state on local analysis domain
  REAL :: invforget                ! inverse forgetting factor
  REAL, ALLOCATABLE :: n_eff(:)    ! Effective sample size for each local domain
  LOGICAL, ALLOCATABLE :: MASK(:)  ! Mask for effective sample sizes > 0
  REAL :: max_n_eff_l, min_n_eff_l ! PE-local min/max. effective ensemble sizes
  REAL :: max_n_eff, min_n_eff     ! Global min/max. effective ensemble sizes
  INTEGER :: cnt_small_svals       ! Counter for small values
  INTEGER :: subtype_dummy         ! Dummy variable to avoid compiler warning
  REAL :: avg_n_eff_l, avg_n_eff   ! Average effective sample size


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- START'

  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(51, 'new')

  fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF
  CALL PDAF_timeit(51, 'old')


! ************************
! *** Inflate ensemble ***
! ************************

  CALL PDAF_timeit(3, 'new')

  do_ensmean = .true.
  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! We need to call the inflation of the forecast ensemble before
     ! the observed ensemble is initialized in PDAFobs_init

     IF (dim_lag==0) THEN

        ! We can apply the inflation here if no smoothing is done
        ! In case of smoothing the inflation is done later
        IF ((type_forget==0 .OR. type_forget==1) .AND. (forget /= 1.0)) THEN

           CALL PDAF_timeit(51, 'new')

           IF (mype == 0 .AND. screen > 0) WRITE (*, '(a, 5x, a, i2, a, f10.3)') &
                'PDAF', 'Inflate forecast ensemble, type_forget=',type_forget,', forget=', forget

           ! Apply forgetting factor
           CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

           ! PDAF_inflate_ens compute the ensmeble mean; thus don't do this in PDAFobs_init
           do_ensmean = .false.

           CALL PDAF_timeit(51, 'old')
        ENDIF
     END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! ***    optionally before call to U_prepoststep    ***
! *****************************************************

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, .true., .true., .true., .true.)
  END IF

  CALL PDAF_timeit(3, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  ! Initialize variable to prevent compiler warning
  subtype_dummy = subtype

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
             'PDAF ', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
  END IF


! ************************
! *** Inflate ensemble ***
! ************************

  CALL PDAF_timeit(3, 'new')

  do_ensmean = .true.
  IF (type_obs_init==1) THEN
     ! We need to call the inflation of the forecast ensemble before
     ! the observed ensemble is initialized in PDAFobs_init

     IF (dim_lag==0) THEN

        ! We can apply the inflation here if no smoothing is done
        ! In case of smoothing the inflation is done later
        IF ((type_forget==0 .OR. type_forget==1) .AND. (forget /= 1.0)) THEN

           CALL PDAF_timeit(51, 'new')

           IF (mype == 0 .AND. screen > 0) WRITE (*, '(a, 5x, a, i2, a, f10.3)') &
                'PDAF', 'Inflate forecast ensemble, type_forget=',type_forget,', forget=', forget

           ! Apply forgetting factor
           CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

           ! PDAF_inflate_ens compute the ensmeble mean; thus don't do this in PDAFobs_init
           do_ensmean = .false.

           CALL PDAF_timeit(51, 'old')
        ENDIF
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
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, do_init_dim_obs, .true., .true., .true.)

     CALL PDAF_timeit(3, 'old')
  END IF


! **************************************
! *** Preparation for local analysis ***
! **************************************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')  ! Time for assimilation 
  CALL PDAF_timeit(7, 'new')  ! global preparation

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(3) dim_lag     ', dim_lag
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(4) type_noise  ', type_noise
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(5) type_forget ', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(6) type_trans  ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(7) type_winf   ', type_winf

     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(1) forget     ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(2) limit_winf ', limit_winf
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(3) noise amp. ', noise_amp
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- call init_n_domains'

  ! Query number of analysis domains for the local analysis
  ! in the PE-local domain
  CALL PDAF_timeit(42, 'new')
  CALL U_init_n_domains_p(step, n_domains_p)
  CALL PDAF_timeit(42, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  n_domains_p', n_domains_p

  ! Initialize effective sample size
  ALLOCATE(n_eff(n_domains_p))
  n_eff = 0.0  !initialize 
  ALLOCATE(MASK(n_domains_p))

  IF (screen > 0) THEN
     IF (mype == 0) THEN
        WRITE (*, '(a, i7, 3x, a)') &
             'PDAF ', step, 'LNETF analysis using T-matrix'
     END IF
     IF (screen<3) THEN
        CALL PDAF_print_domain_stats(n_domains_p)
     ELSE
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
     END IF
  END IF


! *** Local analysis: initialize global quantities ***

  IF (dim_lag>0) THEN
     ! In case of smoothing we first initialize the ensemble without
     ! inflation. Then we apply inflation here and initialize
     ! the observed inflated ensemble

     CALL PDAF_timeit(51, 'new')

     ! For the smoother get observed uninflated ensemble
     ALLOCATE(HX_noinfl_f(dim_obs_f, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f * dim_ens)

     ! This is already initialized by PDAFobs_init
     HX_noinfl_f = HX_f

     CALL PDAF_timeit(51, 'old')

     ! Apply covariance inflation to global ensemlbe
     ! if prior inflation is chosen
     IF (type_forget==0 .OR. type_forget==1) THEN

        CALL PDAF_timeit(51, 'new')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: PDAF_lnetf_update', debug, &
             'Inflate ensemble: type_forget, forget', type_forget, forget

        CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, .false.)

        CALL PDAF_timeit(51, 'old')

        CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_f, &
             state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs_l, &
             screen, debug, .false., .false., .true., .false., .false.)
     END IF

  END IF

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(33, 'new')

  ! Generate orthogonal random matrix with eigenvector (1,...,1)^T

  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  IF (type_trans==0) THEN
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize random transformation'
     CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  rndmat', rndmat
  ELSE
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize deterministic transformation'
     rndmat = 0.0
     DO i = 1, dim_ens
        rndmat(i,i) = 1.0
     END DO
  END IF

  CALL PDAF_timeit(33, 'old')
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(7, 'old')


! ************************
! *** Perform analysis ***
! ************************

  CALL PDAF_timeit(8, 'new')

  ! Initialize counters for statistics on local observations
  CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, TA_l, TA_noinfl_l, flag)

  CALL PDAF_timeit(51, 'new')

  ! Allocate ensemble transform matrix
  ALLOCATE(TA_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)
  TA_l = 0.0

  ! For smoother: Allocate ensemble transform matrix without inflation
  IF (dim_lag>0) THEN
     ALLOCATE(TA_noinfl_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

     TA_noinfl_l = 0.0
  ELSE
     ALLOCATE(TA_noinfl_l(1, 1))
  END IF

  ! initialize number of small singular values
  cnt_small_svals = 0

  CALL PDAF_timeit(51, 'old')

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- Enter local analysis loop'

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysis: DO domain_p = 1, n_domains_p    

     ! Set flag that we are in the local analysis loop
     inloop = .true.


     ! *************************************
     ! *** Initialize local state vector ***
     ! *************************************

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lnetf_update -- local analysis for domain_p', domain_p
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- call init_dim_l'
     END IF

     ! local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  dim_l', dim_l
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- call init_dim_obs_l'
     END IF

     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     IF (allocflag_l == 0) CALL PDAF_memcount(3, 'r', dim_l*dim_ens + dim_l)

     CALL PDAF_timeit(10, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_memberid
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lnetf_update -- call g2l_state for ensemble member', member
           if (member==1) &
                WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lnetf_update --    Note: if ens_l is incorrect check user-defined indices in g2l_state!'
        END IF

        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  ens_l', ens_l(:,member)
     END DO

     ! Store member index to make it accessible with PDAF_get_memberid
     member_save = 0

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_lnetf_update -- call g2l_state for ensemble mean'

     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  meanens_l', state_l

     CALL PDAF_timeit(10, 'old')


     ! *******************************************
     ! *** Initialize local observation arrays ***
     ! *******************************************

     CALL PDAF_timeit(11, 'new')
     CALL PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f, dim_ens, &
          U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, debug)
     CALL PDAF_timeit(11, 'old')


     ! *********************
     ! *** Analysis step ***
     ! *********************

     ! Gather statistical information on local observations
     CALL PDAF_incr_local_obsstats(dim_obs_l)

     IF (dim_obs_l > 0) THEN
        ! OBSERVED DOMAIN

        CALL PDAF_timeit(12, 'new')

        ! Compute NETF analysis
        CALL PDAF_lnetf_ana(domain_p, step, dim_l, dim_obs_l, &
             dim_ens, ens_l, HX_l, obs_l, rndmat, &
             U_likelihood_l, type_forget, forget, &
             type_winf, limit_winf, cnt_small_svals, n_eff(domain_p), TA_l, &
             screen, debug, flag)

        CALL PDAF_timeit(12, 'old')

        ! Compute transform matrix for smoother
        IF (dim_lag>0) THEN
           CALL PDAF_timeit(15, 'new')

           CALL PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l, &
                dim_ens, HX_noinfl_f, rndmat, U_g2l_obs, U_init_obs_l, U_likelihood_l, &
                screen, TA_noinfl_l, flag)

           CALL PDAF_timeit(15, 'old')
        END IF

     ELSE
        ! UNOBSERVED DOMAIN

        CALL PDAF_timeit(51, 'new')

        ! Depending on type_forget, inflation on unobserved domain has to be inverted or applied here
        IF (type_forget==1) THEN
           ! prior inflation NOT on unobserved domains - take it back!
           invforget = 1.0/forget
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_l, ens_l, invforget, .true.)
        ELSEIF (type_forget==2) THEN 
           ! analysis inflation ALSO on unobserved domains - add it!
           invforget = forget
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_l, ens_l, invforget, .true.)
        ENDIF

        CALL PDAF_timeit(51, 'old')

     ENDIF 

     CALL PDAF_timeit(14, 'new')

     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_lnetf_update -- call l2g_state for ensemble member', member
           WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update:', debug, '  ens_l', ens_l(:,member)
        END IF

        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
    
     CALL PDAF_timeit(14, 'old')

     ! *** Perform smoothing of past ensembles ***
     IF (dim_lag>0) THEN
        CALL PDAF_timeit(15, 'new')

        CALL PDAF_smoother_lnetf(domain_p, step, dim_p, dim_l, dim_ens, &
             dim_lag, TA_noinfl_l, ens_l, sens_p, cnt_maxlag, &
             U_g2l_state, U_l2g_state, screen)

        CALL PDAF_timeit(15, 'old')
     END IF


     ! clean up
     DEALLOCATE(ens_l, state_l)
     CALL PDAFobs_dealloc_local()

     ! Set allocflag
     allocflag_l = 1

  END DO localanalysis

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- End of local analysis loop'

  ! Set flag that we are not in the local analysis loop
  inloop = .false.

  ! Clean up arrays allocated in parallel
  DEALLOCATE(TA_l)
  DEALLOCATE(TA_noinfl_l)

!$OMP END PARALLEL

  CALL PDAF_timeit(8, 'old')


  ! *****************************************
  ! *** Perturb particles by adding noise ***
  ! *****************************************

  CALL PDAF_timeit(19, 'new')

  IF (type_noise>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_lnetf_update -- add noise to particles'

     CALL PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p, type_noise, noise_amp, screen)
  END IF

  CALL PDAF_timeit(19, 'old')

  ! Initialize mask array for effective ensemble size
  MASK = (n_eff > 0.0)

  ! *** Print statistics for local analysis to the screen ***
  CALL PDAF_print_local_obsstats(screen)

  IF (npes_filter>1) THEN
     ! Min/max effective sample sizes
     max_n_eff_l = MAXVAL(n_eff)
     CALL MPI_Reduce(max_n_eff_l, max_n_eff, 1, MPI_REALTYPE, MPI_MAX, &
          0, COMM_filter, MPIerr)
     min_n_eff_l = MINVAL(n_eff, MASK)
     CALL MPI_Reduce(min_n_eff_l, min_n_eff, 1, MPI_REALTYPE, MPI_MIN, &
          0, COMM_filter, MPIerr)

     ! Average effective sample size
     avg_n_eff_l = SUM(n_eff)/n_domains_p
     CALL MPI_Reduce(avg_n_eff_l, avg_n_eff, 1, MPI_REALTYPE, MPI_SUM, &
          0, COMM_filter, MPIerr)
     avg_n_eff = avg_n_eff / REAL(npes_filter)
  ELSE
     ! Min/max effective ensemble sizes
     max_n_eff = MAXVAL(n_eff)
     min_n_eff = MINVAL(n_eff, MASK)

     ! Average effective sample sizes
     avg_n_eff = SUM(n_eff)/n_domains_p
  END IF
 
  CALL PDAF_timeit(3, 'old')
 
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
         'PDAF', 'Minimal  effective ensemble size:', min_n_eff
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
          'PDAF', 'Maximal  effective ensemble size:', max_n_eff
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
          'PDAF', 'Average  effective ensemble size:', avg_n_eff
     WRITE (*, '(a, 8x, a, 11x, i6)') &
          'PDAF', 'Number of small singular values:', &
          cnt_small_svals 

     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF
 
! *** Clean up from local analysis update ***
  DEALLOCATE(rndmat)
  IF (dim_lag>0) DEALLOCATE(HX_noinfl_f)
#else
  WRITE (*,'(/5x,a/)') &
       '!!! WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!' 
#endif
 

! *** Poststep for analysis ensemble ***
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

#ifndef PDAF_NO_UPDATE
  DEALLOCATE(n_eff)
#endif

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_update -- END'

END SUBROUTINE PDAFlnetf_update

END MODULE PDAF_lnetf_update
