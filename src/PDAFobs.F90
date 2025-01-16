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
!> Module holding observation-related arrays and methods
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
!!
!! __Revision history:__
!! * 2025-01 - Lars Nerger - Initial code from revising observation handling
!! * Later revisions - see repository log
!!
MODULE PDAFobs

  IMPLICIT NONE
  SAVE

  REAL, ALLOCATABLE :: HX_p(:,:)      ! PE-local observed ensemble
  REAL, ALLOCATABLE :: HXbar_p(:)     ! PE-local observed state
  REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation vector


!-------------------------------------------------------------------------------
  
CONTAINS
!> Initialize full observations
!!
!! This routine collects the operations to initialize the observations,
!! observed ensemble and observed ensemble mean
!!
  SUBROUTINE PDAFobs_initialize(step, dim_p, dim_ens, dim_obs_p, &
       state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       screen, debug)

    USE PDAF_timer, &
         ONLY: PDAF_timeit, PDAF_time_temp
    USE PDAF_mod_filtermpi, &
         ONLY: mype
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filter, &
         ONLY: obs_member, observe_ens
    USE PDAFomi, &
         ONLY: omi_n_obstypes => n_obstypes, omi_omit_obs => omit_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step        ! Current time step
    INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
    INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
    REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
    INTEGER, INTENT(in) :: screen      ! Verbosity flag
    INTEGER, INTENT(in) :: debug       ! Flag for writing debug output

! *** External subroutines 
! ***  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_dim_obs, &      ! Initialize dimension of observation vector
         U_obs_op, &                   ! Observation operator
         U_init_obs                    ! Initialize observation vector

! *** local variables ***
    INTEGER :: i, member, row          ! Counters
    INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
    REAL :: invdimens                  ! Inverse global ensemble size
    REAL, ALLOCATABLE :: resid_p(:)    ! PE-local observation residual


! *********************************
! *** Get observation dimension ***
! *********************************

    IF (mype == 0 .AND. screen > 0) &
         WRITE (*, '(a, 55a)') 'PDAF Prepare observations ', ('-', i = 1, 43)

    IF (debug>0) THEN
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- START'
       WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, '  dim_p', dim_p
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call init_dim_obs'
    END IF

    CALL PDAF_timeit(15, 'new')
    CALL U_init_dim_obs(step, dim_obs_p)
    CALL PDAF_timeit(15, 'old')

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, '  dim_obs_p', dim_obs_p
  
    IF (screen > 2) THEN
       WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
            'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
    END IF


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

    CALL PDAF_timeit(12, 'new')
    haveobs: IF (dim_obs_p > 0) THEN
       ! The residual only exists for domains with observations

       CALL PDAF_timeit(30, 'new')

       ALLOCATE(HX_p(dim_obs_p, dim_ens))
       ALLOCATE(HXbar_p(dim_obs_p))
       ALLOCATE(obs_p(dim_obs_p))
       IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_obs_p + dim_obs_p * dim_ens)

       ! *** Get observed ensemble ***

       IF (debug>0) THEN
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- observe_ens', observe_ens
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call obs_op', dim_ens, 'times'
       END IF

       CALL PDAF_timeit(44, 'new')
       ENS1: DO member = 1, dim_ens
          ! Store member index to make it accessible with PDAF_get_obsmemberid
          obs_member = member

          ! [Hx_1 ... Hx_N]
          CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
       END DO ENS1
       CALL PDAF_timeit(44, 'old')
       CALL PDAF_timeit(30, 'old')

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

       ! *** get vector of observations ***

       IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- call init_obs'

       CALL PDAF_timeit(50, 'new')
       CALL U_init_obs(step, dim_obs_p, obs_p)
       CALL PDAF_timeit(50, 'old')

       ! *** Omit observations with too high innovation ***
       IF (omi_omit_obs)  THEN

          ALLOCATE(resid_p(dim_obs_p))
          IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

          CALL PDAF_timeit(51, 'new')
          resid_p = obs_p - HXbar_p
          CALL PDAF_timeit(51, 'old')

          IF (debug>0) THEN
             WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, &
                  'innovation d(1:min(dim_obs_p,10))', resid_p(1:min(dim_obs_p,10))
             WRITE (*,*) '++ PDAF-debug PDAFobs_initialize:', debug, &
                  'MIN/MAX of innovation', MINVAL(resid_p), MAXVAL(resid_p)
          END IF

          CALL PDAF_timeit(51, 'new')
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
          ALLOCATE(HX_p(1,1))
          ALLOCATE(HXbar_p(1))
          ALLOCATE(obs_p(1))
!          ALLOCATE(resid_p(1))

          ! [Hx_1 ... Hx_N]
          obs_member = 0
          CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)

          DO member = 1, dim_ens
             ! Store member index to make it accessible with PDAF_get_obsmemberid
             obs_member = member

             ! [Hx_1 ... Hx_N]
             CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
          END DO

       END IF
    END IF haveobs
    CALL PDAF_timeit(12, 'old')

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFobs_initialize -- END'

  END SUBROUTINE PDAFobs_initialize

END MODULE PDAFobs
