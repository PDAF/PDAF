! Copyright (c) 2014-2021 Paul Kirchgessner
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
!BOP
!
! !ROUTINE: PDAF_genobs_alloc --- PDAF-internal initialization of GENOBS
!
! !INTERFACE:
SUBROUTINE PDAF_genobs_alloc(subtype, outflag)

! !DESCRIPTION:
! Perform allocation of arrays for observation generation.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger
! Later revisions - see svn log
!
! !USES:
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: screen, dim_ens, dim_p, &
       state, eofU, eofV, sens
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, filterpe, dim_ens_l, task_id

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
  INTEGER, INTENT(out):: outflag        ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_alloc_filters
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER :: allocstat                  ! Status for allocate
  INTEGER :: subtype_dummy              ! Dummy variable to avoid compiler warning


! ******************************
! *** Allocate filter fields ***
! ******************************

  ! Initialize variable to prevent compiler warning
  subtype_dummy = subtype

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

    ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of STATE'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_p)

     ALLOCATE(eofU(1, 1), stat = allocstat)  ! This is not used for GENOBS
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of eofU'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_ens**2)

     ! Allocate full ensemble on filter-PEs
     ALLOCATE(eofV(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of eofV'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_ens)

     ALLOCATE(sens(1, 1, 1), stat = allocstat)  ! Not used in GENOBS
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of sens'
        outflag = 20
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: genobs_alloc - allocate ensmeble array of size ', &
          dim_ens, ' on pe(f) ', mype
     
  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble

     ! Allocate partial ensemble on model-only PEs
     ALLOCATE(eofV(dim_p, dim_ens_l), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of eofV on model-pe'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_ens_l)

     IF (screen > 2) WRITE (*,*) 'PDAF: genobs_alloc - allocate ensemble array of size ', &
          dim_ens_l, ' on pe(m) ', mype_model, ' of model task ',task_id

  END IF on_filterpe

END SUBROUTINE PDAF_genobs_alloc
