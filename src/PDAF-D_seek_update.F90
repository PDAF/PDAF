! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_seek_update.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_update --- Control analysis update of the SEEK filter
!
! !INTERFACE:
SUBROUTINE  PDAF_seek_update(step, dim_p, dim_obs_p, dim_eof, state_p, &
     Uinv, eofV_p, eps, forget, irediag, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
     screen, subtype, incremental, type_forget, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the SEEK filter.
! 
! The analysis is performed by calling PDAF\_seek\_analysis
! and the rediagonalization is performed in PDAF\_seek\_rediag.
! In addition, the routine U\_prepoststep is called prior to
! the analysis and after the rediagonalization to allow the user
! to access the information of the modes and the state estimate.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_eof_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_eof     ! Number of EOFs
  REAL, INTENT(in)    :: forget      ! Forgetting factor
  REAL, INTENT(in)    :: eps         ! Epsilon for approximated TLM evolution
  INTEGER, INTENT(in) :: irediag     ! Interval to perform rediagonalization
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_eof,dim_eof) ! Inverse of matrix U
  REAL, INTENT(inout) :: eofV_p(dim_p,dim_eof) ! PE-local matrix V
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: subtype     ! Filter subtype
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prepoststep, &         ! User supplied pre/poststep routine
       U_prodRinvA              ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_seek
! Calls: U_prepoststep
! Calls: PDAF_seek_analysis
! Calls: PDAF_seek_rediag
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
!EOP

! *** local variables ***
  INTEGER :: i, col    ! Counters
  INTEGER :: minusStep ! Time step counter
  REAL :: epsinv       ! inverse of epsilon
  INTEGER :: countstep = 1   ! Internal step counter
  REAL, ALLOCATABLE :: Uinv_dyn(:,:) ! temporary matrix if Uinv is kept static



! **********************************************
! ***  Compute evolved basis of error space  ***
! **********************************************

  IF (subtype /= 2 .AND. subtype /= 3 .AND. subtype /= 5) THEN
     ! Do not do mode-ensemble handling for fixed-basis variants
     epsinv = 1.0 / eps
  
     DO  col = 1, dim_eof
        DO i = 1, dim_p
           eofV_p(i, col) = epsinv * (eofV_p(i, col) - state_p(i))
        END DO
     END DO
  END IF


! **********************
! ***  Update phase  ***
! **********************

! *** prestep for forecast modes and state ***
  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_eof, dim_eof_l, dim_obs_p, &
       state_p, Uinv, eofV_p, flag)
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')
  ! *** SEEK analysis with forgetting factor ***
  subt: IF (subtype /= 3) THEN
     CALL PDAF_seek_analysis(step, dim_p, dim_obs_p, dim_eof, state_p, &
          Uinv, eofV_p, forget, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prodRinvA, screen, incremental, type_forget, &
          flag)
  ELSE subt
     ! For fixed Uinv initialize dynamic Uinv which is hold only here
     ALLOCATE(Uinv_dyn(dim_eof, dim_eof))
     Uinv_dyn = Uinv

     ! Perform analysis
     CALL PDAF_seek_analysis(step, dim_p, dim_obs_p, dim_eof, state_p, &
          Uinv_dyn, eofV_p, forget, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prodRinvA, screen, incremental, type_forget, &
          flag)
  END IF subt
  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

! *** Re-orthogonalize the covariance modes ***
  re_diag: IF (irediag > 0) THEN
     re_diag2: IF ( MOD(countstep, irediag) == 0) THEN
        CALL PDAF_timeit(4, 'new')
        CALL PDAF_seek_rediag(dim_p, dim_eof, Uinv, eofV_p, subtype, &
             screen, flag)
        CALL PDAF_timeit(4, 'old')
        IF (mype == 0 .AND. screen > 1) THEN
           WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
                'PDAF', '--- re-diag duration:', PDAF_time_temp(4), 's'
        END IF
     END IF re_diag2
  END IF re_diag

  ! increment stepping counter
  countstep = countstep + 1
#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif


! *** poststep for analysis modes and state ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  IF (subtype /= 3) THEN
     CALL U_prepoststep(step, dim_p, dim_eof, dim_eof_l, dim_obs_p, &
          state_p, Uinv, eofV_p, flag)
  ELSE
     ! prepoststep with fixed Uinv - hand over Uinv as changed by analysis
     CALL U_prepoststep(step, dim_p, dim_eof, dim_eof_l, dim_obs_p, &
          state_p, Uinv_dyn, eofV_p, flag)
     DEALLOCATE(Uinv_dyn)
  END IF
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF

END SUBROUTINE PDAF_seek_update
