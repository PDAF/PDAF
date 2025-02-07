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
!> Module holding observation-related arrays and methods
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
!!
!! __Revision history:__
!! * 2025-01 - Lars Nerger - Initial code from revising observation handling
!! * Later revisions - see repository log
!!
MODULE PDAFobs

  IMPLICIT NONE
  SAVE

  REAL, ALLOCATABLE :: HX_p(:,:)      !< PE-local or full observed ensemble
  REAL, ALLOCATABLE :: HXbar_p(:)     !< PE-local or full observed state
  REAL, ALLOCATABLE :: obs_p(:)       !< PE-local or_full observation vector
  INTEGER :: type_obs_init=1          !< Set at which time the observations are initialized
                                      !< (0) before, (1) after, (2) before and after call to U_prepoststep
  LOGICAL :: observe_ens=.false.      !< (F) to apply H to ensemble mean to compute innovation
                                      !< or (T) apply H to X, compute mean of HX and then residual

  ! Variables for domain-local filters
  REAL, ALLOCATABLE :: HX_l(:,:)      !< Local observed ensemble 
  REAL, ALLOCATABLE :: HXbar_l(:)     !< Local observed ensemble mean 
  REAL, ALLOCATABLE :: obs_l(:)       !< Local observation vector 

!$OMP THREADPRIVATE(HX_l, HXbar_l, obs_l)

!-------------------------------------------------------------------------------
  
CONTAINS
!> Initialize full observations
!!
!! This routine collects the operations to initialize the observations,
!! observed ensemble and observed ensemble mean
!!
  SUBROUTINE PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
       state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       screen, debug, &
       do_ens_mean, do_init_dim, do_HX, do_HXbar, do_init_obs)

    USE PDAF_timer, &
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_mod_filtermpi, &
         ONLY: mype
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filter, &
         ONLY: obs_member, localfilter
    USE PDAFomi, &
         ONLY: omi_n_obstypes => n_obstypes, omi_omit_obs => omit_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step        !< Current time step
    INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble
    INTEGER, INTENT(inout) :: dim_obs_p  !< PE-local dimension of observation vector
    REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
    INTEGER, INTENT(in) :: screen      !< Verbosity flag
    INTEGER, INTENT(in) :: debug       !< Flag for writing debug output
    LOGICAL, INTENT(in) :: do_ens_mean !< Whether to compute ensemble mean
    LOGICAL, INTENT(in) :: do_init_dim !< Whether to call U_init_dim_obs
    LOGICAL, INTENT(in) :: do_HX       !< Whether to initialize HX_p
    LOGICAL, INTENT(in) :: do_HXbar    !< Whether to initialize HXbar
    LOGICAL, INTENT(in) :: do_init_obs !< Whether to initialize obs_p

! *** External subroutines 
! ***  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_dim_obs, &      !< Initialize dimension of observation vector
         U_obs_op, &                   !< Observation operator
         U_init_obs                    !< Initialize observation vector

! *** local variables ***
    INTEGER :: i, member, row          ! Counters
    INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
    REAL :: invdimens                  ! Inverse global ensemble size
    REAL, ALLOCATABLE :: resid_p(:)    ! PE-local observation residual


! ***********************************
! *** Compute ensemble mean state ***
! ***********************************

    IF (do_ens_mean) THEN
       CALL PDAF_timeit(51, 'old')
       CALL PDAF_timeit(9, 'new')

       state_p = 0.0
       invdimens = 1.0 / REAL(dim_ens)
       DO member = 1, dim_ens
          DO row = 1, dim_p
             state_p(row) = state_p(row) + invdimens * ens_p(row, member)
          END DO
       END DO
  
       CALL PDAF_timeit(9, 'old')
       CALL PDAF_timeit(51, 'old')
    END IF


! *********************************
! *** Get observation dimension ***
! *********************************

    CALL PDAF_timeit(6, 'new')

    IF (mype == 0 .AND. screen > 0 .AND. do_init_dim) &
         WRITE (*, '(a, 55a)') 'PDAF Prepare observations ', ('-', i = 1, 43)

    IF (debug>0) THEN
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- START'
       WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, '  dim_p', dim_p
       IF (do_init_dim) WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call init_dim_obs'
    END IF

    IF (do_init_dim) THEN
       CALL PDAF_timeit(43, 'new')
       CALL U_init_dim_obs(step, dim_obs_p)
       CALL PDAF_timeit(43, 'old')
    END IF

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, '  dim_obs_p', dim_obs_p
  
    IF (screen > 2) THEN
       WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
            'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
    END IF


! **********************************************************************
! *** Initialize                                                     ***
! *** - observed ensemble mean (HXbar_p)                             ***
! *** - observed ensemble (HX_p)                                     ***
! *** - observation vector (obs_p)                                   ***
! *** also initialize omission of observations with large innovation ***
! **********************************************************************

    haveobs: IF (dim_obs_p > 0) THEN
       ! Full initialization if observations exist on this process domain

       IF (.NOT.ALLOCATED(HX_p)) ALLOCATE(HX_p(dim_obs_p, dim_ens))
       IF (.NOT.ALLOCATED(HXbar_p)) ALLOCATE(HXbar_p(dim_obs_p))
       IF (.NOT.ALLOCATED(obs_p)) ALLOCATE(obs_p(dim_obs_p))
       IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_obs_p + dim_obs_p * dim_ens)

       ! *** Get observed ensemble ***

       IF (debug>0) THEN
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- observe_ens', observe_ens
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call obs_op', dim_ens, 'times'
       END IF

       IF (do_HX) THEN
          CALL PDAF_timeit(44, 'new')
          ENS1: DO member = 1, dim_ens
             ! Store member index to make it accessible with PDAF_get_obsmemberid
             obs_member = member

             ! [Hx_1 ... Hx_N]
             CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
          END DO ENS1
          CALL PDAF_timeit(44, 'old')
       END IF

       IF (do_HXbar) THEN
          IF (.NOT.observe_ens) THEN

             ! Directly obtain observed ensemble mean
             IF (debug>0) &
                  WRITE (*,*) '++ PDAF-debug: ', debug, &
                  'PDAFobs_initialize -- call obs_op for ensemble mean'

             obs_member = 0 ! Store member index (0 for central state)
             CALL PDAF_timeit(44, 'new')
             CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)
             CALL PDAF_timeit(44, 'old')
          ELSE
             ! Compute observed mean from observed ensemble
             ! This is more accurate if H is nonlinear

             CALL PDAF_timeit(51, 'new')
             HXbar_p = 0.0
             invdimens = 1.0 / REAL(dim_ens)
             DO member = 1, dim_ens
                DO row = 1, dim_obs_p
                   HXbar_p(row) = HXbar_p(row) + invdimens * HX_p(row, member)
                END DO
             END DO
             CALL PDAF_timeit(51, 'old')
          END IF
       END IF

       ! *** get vector of observations ***

       IF (do_init_obs) THEN
          IF (debug>0) &
               WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call init_obs'

          CALL PDAF_timeit(50, 'new')
          CALL U_init_obs(step, dim_obs_p, obs_p)
          CALL PDAF_timeit(50, 'old')
       END IF

       ! *** Omit observations with too high innovation ***

       IF (omi_omit_obs .AND. do_HXbar .AND. do_init_obs &
            .AND. localfilter/=1)  THEN

          ALLOCATE(resid_p(dim_obs_p))
          IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

          CALL PDAF_timeit(51, 'new')
          resid_p = obs_p - HXbar_p

          IF (debug>0) THEN
             WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, &
                  'innovation d(1:min(dim_obs_p,10))', resid_p(1:min(dim_obs_p,10))
             WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, &
                  'MIN/MAX of innovation', MINVAL(resid_p), MAXVAL(resid_p)
          END IF

          CALL PDAFomi_omit_by_inno_cb(dim_obs_p, resid_p, obs_p)
          CALL PDAF_timeit(51, 'old')

          DEALLOCATE(resid_p)
       END IF

    ELSE IF (dim_obs_p == 0) THEN

       ! For OMI we need to call observation operator also for dim_obs_p=0
       ! in order to initialize the pointer to the observation types
       ! Further the observation operator has to be executed in cases
       ! in which the operation include a global communication
       IF (omi_n_obstypes>0) THEN
          IF (.NOT.ALLOCATED(HX_p)) ALLOCATE(HX_p(1,1))
          IF (.NOT.ALLOCATED(HXbar_p)) ALLOCATE(HXbar_p(1))
          IF (.NOT.ALLOCATED(obs_p)) ALLOCATE(obs_p(1))

          ! [Hx_1 ... Hx_N]
          IF (do_HXbar) THEN
             obs_member = 0
             CALL PDAF_timeit(44, 'new')
             CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)
             CALL PDAF_timeit(44, 'old')
          END IF

          IF (do_HX) THEN
             CALL PDAF_timeit(44, 'new')
             DO member = 1, dim_ens
                ! Store member index to make it accessible with PDAF_get_obsmemberid
                obs_member = member

                ! [Hx_1 ... Hx_N]
                CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
             END DO
             CALL PDAF_timeit(44, 'old')
          END IF

       END IF
    END IF haveobs
    CALL PDAF_timeit(6, 'old')

    IF (mype == 0 .AND. screen > 1) THEN
       WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
            'PDAF', '--- duration of observation preparation:', PDAF_time_temp(6), 's'
    END IF
    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- END'

  END SUBROUTINE PDAFobs_init



!-------------------------------------------------------------------------------
!> Initialize local observation arrays
!!
!! This routine collects the operations to initialize the 
!! local observations, observed ensemble and observed ensemble mean
!! for domain-local filters
!!
  SUBROUTINE PDAFobs_init_local(domain_p, step, dim_obs_l, dim_obs_f, dim_ens, &
       U_init_dim_obs_l, U_g2l_obs, U_init_obs_l, debug)

    USE PDAF_timer, &
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filter, &
         ONLY: obs_member
    USE PDAFomi, &
         ONLY: omi_omit_obs => omit_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: domain_p    !< Current local analysis domain
    INTEGER, INTENT(in) :: step        !< Current time step
    INTEGER, INTENT(out) :: dim_obs_l  !< Size of local observation vector
    INTEGER, INTENT(out) :: dim_obs_f  !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble 
    INTEGER, INTENT(in) :: debug       !< Flag for writing debug output

! *** External subroutines 
! ***  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_obs_l, &        !< Init. observation vector on local analysis domain
         U_g2l_obs, &                  !< Restrict full obs. vector to local analysis domain
         U_init_n_domains_p            !< Provide number of local analysis domains

! *** Local variables ***
    INTEGER :: member                  ! Counter
    INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
    REAL, ALLOCATABLE :: innov_l(:)    ! Local innovation


! *******************************************
! *** Initialize local observation arrays ***
! *******************************************

    ! *** Get observation dimension for local domain ***
    CALL PDAF_timeit(38, 'new')
    dim_obs_l = 0
    CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
    CALL PDAF_timeit(38, 'old')

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  dim_obs_l', dim_obs_l

     haveobs: IF (dim_obs_l > 0) THEN

        ALLOCATE(obs_l(dim_obs_l))
        ALLOCATE(HXbar_l(dim_obs_l))
        ALLOCATE(HX_l(dim_obs_l, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_obs_l + dim_obs_l*dim_ens)

        ! *** Get local observed ensemble mean state ***

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call g2l_obs for mean'

        CALL PDAF_timeit(46, 'new')
        obs_member = 0
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_p, HXbar_l)
        CALL PDAF_timeit(46, 'old')


        ! *** Get local observed state ensemble ***

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_analysis -- call g2l_obs', dim_ens, 'times'

        CALL PDAF_timeit(46, 'new')
        ENS: DO member = 1, dim_ens
           ! Store member index
           obs_member = member

           ! [Hx_1 ... Hx_N] for local analysis domain
           CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_p(:, member), &
                HX_l(:, member))
        END DO ENS
        CALL PDAF_timeit(46, 'old')


        ! *** Get local observation vector ***

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_obs_l'

        CALL PDAF_timeit(47, 'new')
        CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
        CALL PDAF_timeit(47, 'old')
        
        ! *** Optional case: omit obserations with too high innovation

        IF (omi_omit_obs) THEN

           CALL PDAF_timeit(51, 'new')

           ALLOCATE(innov_l(dim_obs_l))

           ! Compute local innovation
           innov_l = obs_l - HXbar_l

           ! Omit observations with too high innovation
           CALL PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, innov_l, obs_l)

           DEALLOCATE(innov_l)

           CALL PDAF_timeit(51, 'old')
        END IF

     ELSE

        ALLOCATE(obs_l(1))
        ALLOCATE(HX_l(1,1))
        ALLOCATE(HXbar_l(1))

     END IF haveobs


! ********************
! *** Finishing up ***
! ********************

     IF (allocflag == 0) allocflag = 1

   END SUBROUTINE PDAFobs_init_local

!-------------------------------------------------------------------------------
!> Deallocate observation arrays
!!
!! This routine deallocates the observation-related arrays
!! that were allocated in PDAFomi_initialize
!!
  SUBROUTINE PDAFobs_dealloc()

    IMPLICIT NONE

    IF (ALLOCATED(HX_p)) DEALLOCATE(HX_p)
    IF (ALLOCATED(HXbar_p)) DEALLOCATE(HXbar_p)
    IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)

  END SUBROUTINE PDAFobs_dealloc

!-------------------------------------------------------------------------------
!> Deallocate local observation arrays
!!
!! This routine deallocates the observation-related arrays
!! that were allocated in PDAFomi_initialize
!!
  SUBROUTINE PDAFobs_dealloc_local()

    IMPLICIT NONE

    IF (ALLOCATED(HX_l)) DEALLOCATE(HX_l)
    IF (ALLOCATED(HXbar_l)) DEALLOCATE(HXbar_l)
    IF (ALLOCATED(obs_l)) DEALLOCATE(obs_l)

  END SUBROUTINE PDAFobs_dealloc_local

END MODULE PDAFobs
