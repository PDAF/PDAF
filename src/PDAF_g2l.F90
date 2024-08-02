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
! !ROUTINE: PDAF_g2l - Project global to local vector according to index array
!
! !INTERFACE:
SUBROUTINE PDAF_g2l(dim_p, dim_l, idx_l_in_p, state_p, state_l)

! !DESCRIPTION:
! Project a global to a local state vector for the localized filters.
! The mapping is done using the provided index array.
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: idx_l_in_p(dim_l)     !< Index array for projection
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by user code
!EOP
  
! *** local variables ***
  INTEGER :: i                  ! Counter


! *************************************
! *** Initialize local state vector ***
! *************************************

  DO i = 1, dim_l
     ! Only use the index value if it is in a valid range
     IF (idx_l_in_p(i) > 0 .and. idx_l_in_p(i) <= dim_p) THEN
        state_l(i) = state_p(idx_l_in_p(i))
     END IF
  END DO
   
END SUBROUTINE PDAF_g2l
