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
! !ROUTINE: PDAF_reset_dim_p --- Reset state dimension and re-allocate state and ensemble
!
! !INTERFACE:
SUBROUTINE PDAF_reset_dim_p(dim_p_in, outflag)

! !DESCRIPTION:
! Reset state dimension and re-allocate the state vector
! and ensemble array.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2024-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mpi
!  USE PDAF_memcounting, &
!       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: screen, incremental, dim_ens, dim_p, &
       state, state_inc, eofV
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, filterpe, dim_ens_l, task_id, &
       COMM_couple

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p_in       ! Sub-type of filter
  INTEGER, INTENT(inout):: outflag      ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_alloc_filters
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER :: allocstat                  ! Status for allocate


! ************************************
! *** Reset state vector dimension ***
! ************************************

  dim_p = dim_p_in


! *********************************
! *** Re-allocate filter fields ***
! *********************************

  ! Initialize status flag
  outflag = 0
  
  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     IF (ALLOCATED(state)) DEALLOCATE(state)
     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in re-allocation of STATE'
        outflag = 20
     END IF

     IF (incremental == 1) THEN
        IF (ALLOCATED(state_inc)) DEALLOCATE(state_inc)
        ALLOCATE(state_inc(dim_p), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE_INC'
           outflag = 20
        END IF

        state_inc = 0.0
     END IF

     ! Allocate full ensemble on filter-PEs
     IF (ALLOCATED(eofV)) DEALLOCATE(eofV)
     ALLOCATE(eofV(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV'
        outflag = 20
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_p - re-allocate eofV of size ', &
          dim_ens, ' on pe(f) ', mype

  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        IF (ALLOCATED(eofV)) DEALLOCATE(eofV)
        ALLOCATE(eofV(dim_p, dim_ens_l), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV on model-pe'
           outflag = 20
        END IF

        IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_p - re-allocate eofV of size ', &
             dim_ens_l, ' on pe(m) ', mype_model, ' of model task ',task_id
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_reset_dim_p
