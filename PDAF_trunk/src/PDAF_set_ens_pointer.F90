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
!$Id: PDAF_set_ens_pointer.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_set_ens_pointer --- Set pointer to ensemble array
SUBROUTINE PDAF_set_ens_pointer(ens_point, status)

! !DESCRIPTION:
! Routine to set the pointer to the PDAF-internal ensemble array.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2016-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: eofV

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, POINTER, INTENT(out) :: ens_point(:,:)  ! Pointer to smoother array
  INTEGER, INTENT(out)       :: status  ! Status flag

! !CALLING SEQUENCE:
! Called by: model code
!EOP

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(eofV)) THEN
     ens_point => eofV

     status = 0
  ELSE
     status = 1
  END IF
  
END SUBROUTINE PDAF_set_ens_pointer
