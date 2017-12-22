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
!$Id: PDAF-D_alloc_filters.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_alloc_filters --- internal interface to allocation routines
!
! !INTERFACE:
SUBROUTINE PDAF_alloc_filters(filterstr, subtype, flag)

! !DESCRIPTION:
! Interface routine to the filter-specific allocation
! routines.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-08 - Lars Nerger - Initial code for restructuring PDAF
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=10), INTENT(in) :: filterstr ! Name of filter algorithm
  INTEGER, INTENT(in) :: subtype             ! Sub-type of filter
  INTEGER, INTENT(inout)::flag               ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init
! Calls PDAF_seek_alloc
! Calls PDAF_seik_alloc
! Calls PDAF_enkf_alloc
! Calls PDAF_lseik_alloc
! Calls PDAF_etkf_alloc
! Calls PDAF_letkf_alloc
!EOP


! ***********************************************
! *** Call filter-specific allocation routine ***
! ***********************************************

  checkflag: IF (flag == 0) THEN
     IF (TRIM(filterstr) == 'SEEK') THEN
        CALL PDAF_seek_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'SEIK') THEN
        CALL PDAF_seik_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'ENKF') THEN
        CALL PDAF_enkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
        CALL PDAF_lseik_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'ETKF') THEN
        CALL PDAF_etkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'LETKF') THEN
        CALL PDAF_letkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
        CALL PDAF_estkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
        CALL PDAF_lestkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'LENKF') THEN
        CALL PDAF_lenkf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'NETF') THEN
        CALL PDAF_netf_alloc(subtype, flag)

     ELSE IF (TRIM(filterstr) == 'LNETF') THEN
        CALL PDAF_lnetf_alloc(subtype, flag)

     ENDIF
  END IF checkflag

END SUBROUTINE PDAF_alloc_filters
