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
!$Id: PDAF-D_smoother_shift.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_smoother_shift --- Shift ensemble states in ensemble array for smoothing
!
! !INTERFACE:
SUBROUTINE PDAF_smoother_shift(dim_p, dim_ens, dim_lag, ens_p, sens_p, cnt_maxlag, screen)

! !DESCRIPTION:
! Routine to store a previous analysis ensemble for smoothing.
! The storage is performed by shifting all previous ensembles
! by one. Thus the previous ensembles are stored consecutively.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2012-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p         ! PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       ! Number of past time instances for smoother
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens, 1)   ! PE-local state ensemble
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)  ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count available number of time steps for smoothing
  INTEGER, INTENT(in) :: screen        ! Verbosity flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state
!EOP

! *** local variables ***
  INTEGER :: lag             ! Counters
  INTEGER :: n_lags          ! Available number of tiem instances for smoothing
  
  
! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoothe over full lag
     n_lags = cnt_maxlag
  END IF

  IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Store previous analysis for smoother'
  END IF


! ***********************
! *** Shift ensembles ***
! ***********************

  ! Shift past ensembles
  DO lag = n_lags-1, 1, -1

     IF (mype == 0 .AND. screen > 2) &
          write (*,*) 'PDAF: smoother: shift column', lag, 'to ',lag+1

     sens_p(:,:,lag+1) = sens_p(:,:,lag)

  END DO

  ! Store current ensemble
  IF (n_lags > 0) THEN

     IF (mype == 0 .AND. screen > 2) &
          write (*,*) 'PDAF: smoother: store current ensemble in smoother array'

     sens_p(:,:,1) = ens_p(:,:,1)
  END IF


! ********************
! *** Finishing up ***
! ********************

END SUBROUTINE PDAF_smoother_shift
