! Copyright (c) 2004-2021 Lars Nerger
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
! !ROUTINE: PDAF_deallocate --- deallocate PDAF-internal arrays
!
! !INTERFACE:
SUBROUTINE PDAF_deallocate()

! !DESCRIPTION:
! Perform deallocation of PDAF-internal arrays
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2017-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE PDAF_mod_filter, &
       ONLY: dim_bias_p, state, state_inc, eofU, eofV, &
       sens, bias, dim_lag
  USE PDAF_mod_filtermpi, &
       ONLY: filterpe, COMM_couple

  IMPLICIT NONE

! !ARGUMENTS:

! !CALLING SEQUENCE:
! Called by: PDAF_deallocate
!EOP

! *** local variables ***

! ******************************
! *** Allocate filter fields ***
! ******************************

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     DEALLOCATE(state)

     IF (ALLOCATED(state_inc)) DEALLOCATE(state_inc)

     DEALLOCATE(eofU)

     ! Allocate full ensemble on filter-PEs
     DEALLOCATE(eofV)

     ! Allocate array for past ensembles for smoothing on filter-PEs
     IF (dim_lag > 0) THEN
        DEALLOCATE(sens)
     END IF

     IF (dim_bias_p > 0) THEN
        DEALLOCATE(bias)
     ENDIF
     
  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        DEALLOCATE(eofV)
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_deallocate
