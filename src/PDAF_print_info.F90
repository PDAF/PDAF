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
!> Print information for PDAF (timing and memory) to screen
!!
!! This routine displays the information from PDAF.
!! Possible are to display the timing information and
!! allocated memory.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2008-09 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_print_info(printtype)

  USE PDAF_mod_filter, ONLY: filterstr
  USE PDAF_seek, ONLY: PDAF_seek_memtime
  USE PDAF_seik, ONLY: PDAF_seik_memtime
  USE PDAF_lseik, ONLY: PDAF_lseik_memtime
  USE PDAF_enkf, ONLY: PDAF_enkf_memtime
  USE PDAF_lenkf, ONLY: PDAF_lenkf_memtime
  USE PDAF_estkf, ONLY: PDAF_estkf_memtime
  USE PDAF_lestkf, ONLY: PDAF_lestkf_memtime
  USE PDAF_etkf, ONLY: PDAF_etkf_memtime
  USE PDAF_letkf, ONLY: PDAF_letkf_memtime
  USE PDAF_netf, ONLY: PDAF_netf_memtime
  USE PDAF_lnetf, ONLY: PDAF_lnetf_memtime
  USE PDAF_lknetf, ONLY: PDAF_lknetf_memtime
  USE PDAF_pf, ONLY: PDAF_pf_memtime
  USE PDAF_3dvar, ONLY: PDAF_3dvar_memtime


  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: printtype       !< Type of screen output:  
                                         !< (1) timings, (2) memory


! ********************************
! *** Print screen information ***
! ********************************

  IF (TRIM(filterstr) == 'SEEK') THEN
     CALL PDAF_seek_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'SEIK') THEN
     CALL PDAF_seik_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'ENKF') THEN
     CALL PDAF_enkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_lseik_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'ETKF') THEN
     CALL PDAF_etkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_letkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_estkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_lestkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LENKF') THEN
     CALL PDAF_lenkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_netf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_lnetf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_lknetf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_pf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_3dvar_memtime(printtype)
  END IF

END SUBROUTINE PDAF_print_info
