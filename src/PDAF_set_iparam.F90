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
!> PDAF_set_iparam --- Set integer parameter for PDAF
!!
!! This routine lets the user set the value of a
!! method-specific integer parameter. This routine
!! builds the interface to the specific routine
!! provided by each DA method.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_set_iparam(id, value, flag)

  USE PDAF_mod_filter, &
       ONLY: filterstr
  USE PDAF_seik, &
       ONLY: PDAF_seik_set_iparam
  USE PDAF_enkf, &
       ONLY: PDAF_enkf_set_iparam
  USE PDAF_lseik, &
       ONLY: PDAF_lseik_set_iparam
  USE PDAF_etkf, &
       ONLY: PDAF_etkf_set_iparam
  USE PDAF_letkf, &
       ONLY: PDAF_letkf_set_iparam
  USE PDAF_estkf, &
       ONLY: PDAF_estkf_set_iparam
  USE PDAF_lestkf, &
       ONLY: PDAF_lestkf_set_iparam
  USE PDAF_lenkf, &
       ONLY: PDAF_lenkf_set_iparam
  USE PDAF_netf, &
       ONLY: PDAF_netf_set_iparam
  USE PDAF_lnetf, &
       ONLY: PDAF_lnetf_set_iparam
  USE PDAF_lknetf, &
       ONLY: PDAF_lknetf_set_iparam
  USE PDAF_pf, &
       ONLY: PDAF_pf_set_iparam
  USE PDAF_3dvar, &
       ONLY: PDAF_3dvar_set_iparam
  USE PDAF_genobs, &
       ONLY: PDAF_genobs_set_iparam

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT(in) :: id       !< Index of parameter
  INTEGER, INTENT(in) :: value    !< Parameter value
  INTEGER, INTENT(out) :: flag    !< Status flag: 0 for no error


! ********************************
! *** Print screen information ***
! ********************************
write (*,*) 'set_iparam, filterstr', filterstr
write (*,*) 'set_iparam: id', id,' value', value
  IF (TRIM(filterstr) == 'SEIK') THEN
     CALL PDAF_seik_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'ENKF') THEN
     CALL PDAF_enkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_lseik_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'ETKF') THEN
     CALL PDAF_etkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_letkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_estkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_lestkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LENKF') THEN
     CALL PDAF_lenkf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_netf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_lnetf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_lknetf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_pf_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_3dvar_set_iparam(id, value, flag)
  ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
     CALL PDAF_genobs_set_iparam(id, value, flag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: invalid DA method - likely PDAF is not yet initialized' 
  END IF

END SUBROUTINE PDAF_set_iparam
