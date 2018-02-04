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
!$Id: PDAF-D_get_smootherens.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_get_smootherens --- Set pointer to smoother ensemble
SUBROUTINE PDAF_get_smootherens(sens_point, maxlag, status)

! !DESCRIPTION:
! Routine to set the pointer to the PDAF-internal smoother ensemble array.
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
  INTEGER, INTENT(out)       :: maxlag  ! Number of past timesteps processed in sens
  INTEGER, INTENT(out)       :: status  ! Status flag 

! !CALLING SEQUENCE:
! Called by: U_prepoststep
!EOP

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(sens)) THEN
     sens_point => sens

     status = 0
  END IF

  ! Set number of initialized lags
  IF (cnt_maxlag-1 >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     maxlag = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     maxlag = cnt_maxlag-1
  END IF

END SUBROUTINE PDAF_get_smootherens
