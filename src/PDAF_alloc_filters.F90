! Copyright (c) 2004-2025 Lars Nerger
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
!
!> Interface routine to the filter-specific allocation routines.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2010-08 - Lars Nerger - Initial code from restructuring PDAF
!! Later revisions - see svn log
!!
SUBROUTINE PDAF_alloc_filters(filterstr, subtype, flag)

  USE PDAF_seek, ONLY: PDAF_seek_alloc
  USE PDAF_seik, ONLY: PDAF_seik_alloc
  USE PDAF_lseik, ONLY: PDAF_lseik_alloc
  USE PDAF_enkf, ONLY: PDAF_enkf_alloc
  USE PDAF_etkf, ONLY: PDAF_etkf_alloc
  USE PDAF_letkf, ONLY: PDAF_letkf_alloc
  USE PDAF_estkf, ONLY: PDAF_estkf_alloc
  USE PDAF_lestkf, ONLY: PDAF_lestkf_alloc
  USE PDAF_lenkf, ONLY: PDAF_lenkf_alloc
  USE PDAF_netf, ONLY: PDAF_netf_alloc
  USE PDAF_lnetf, ONLY: PDAF_lnetf_alloc
  USE PDAF_lknetf, ONLY: PDAF_lknetf_alloc
  USE PDAF_pf, ONLY: PDAF_pf_alloc
  USE PDAF_genobs, ONLY: PDAF_genobs_alloc
  USE PDAF_3dvar, ONLY: PDAF_3dvar_alloc

  IMPLICIT NONE

! *** Arguments ***
  CHARACTER(len=10), INTENT(in) :: filterstr ! Name of filter algorithm
  INTEGER, INTENT(in) :: subtype             ! Sub-type of filter
  INTEGER, INTENT(inout)::flag               ! Status flag


! ***********************************************
! *** Call filter-specific allocation routine ***
! ***********************************************

  checkflag: IF (flag == 0) THEN
     IF (TRIM(filterstr) == 'SEEK') THEN
        CALL PDAF_seek_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'SEIK') THEN
        CALL PDAF_seik_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'ENKF') THEN
        CALL PDAF_enkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
        CALL PDAF_lseik_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'ETKF') THEN
        CALL PDAF_etkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LETKF') THEN
        CALL PDAF_letkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
        CALL PDAF_estkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
        CALL PDAF_lestkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LENKF') THEN
        CALL PDAF_lenkf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'NETF') THEN
        CALL PDAF_netf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LNETF') THEN
        CALL PDAF_lnetf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
        CALL PDAF_lknetf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'PF') THEN
        CALL PDAF_pf_alloc(flag)

     ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
        CALL PDAF_genobs_alloc(flag)

     ELSE IF (TRIM(filterstr) == '3DVAR') THEN
        CALL PDAF_3dvar_alloc(subtype, flag)

     ENDIF
  END IF checkflag

END SUBROUTINE PDAF_alloc_filters
