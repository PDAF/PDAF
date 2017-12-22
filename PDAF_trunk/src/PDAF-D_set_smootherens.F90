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
!$Id: PDAF-D_set_smootherens.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_set_smootherens --- Set pointer to smoother ensemble
SUBROUTINE PDAF_set_smootherens(sens_point, maxlag, status)

! !DESCRIPTION:
! Routine to set the pointer to the PDAF-internal smoother ensemble array
! and to set the value of the maximum lah to be considered.
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

  USE PDAF_mod_filter, &
       ONLY: sens, cnt_maxlag, dim_lag

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
  INTEGER, INTENT(in)        :: maxlag  ! Number of past timesteps in sens
  INTEGER, INTENT(out)       :: status  ! Status flag, 
                                        ! 0: no error, 1: maxlag too large

! !CALLING SEQUENCE:
! Called by: PDAF ensemble initialization routine
!EOP

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(sens)) THEN
     sens_point => sens

     status = 0
  END IF

  
! **************************************
! *** Set number of initialized lags ***
! **************************************

  IF (maxlag <= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     cnt_maxlag = maxlag
     
     status = 0
  ELSE
     ! Maxlag is larger than allocated smoother ensemble array
     cnt_maxlag = dim_lag

     status = 1
  END IF

END SUBROUTINE PDAF_set_smootherens
