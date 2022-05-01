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
! !ROUTINE: PDAF_get_ensstats --- Set pointer to ensemble statistics
SUBROUTINE PDAF_get_ensstats(skew_ptr, kurt_ptr, status)

! !DESCRIPTION:
! Routine to set the pointer to the PDAF-internal array of skewness and kurtosis.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2020-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: skewness, kurtosis

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, POINTER, INTENT(out) :: skew_ptr(:)  ! Pointer to skewness array
  REAL, POINTER, INTENT(out) :: kurt_ptr(:)  ! Pointer to kurtosis array
  INTEGER, INTENT(out)       :: status  ! Status flag 

! !CALLING SEQUENCE:
! Called by: U_prepoststep
!EOP

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(skewness)) THEN
     skew_ptr => skewness
     kurt_ptr => kurtosis

     status = 0
  END IF

END SUBROUTINE PDAF_get_ensstats
