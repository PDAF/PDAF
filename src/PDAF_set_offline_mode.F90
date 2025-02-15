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
!$Id$
!BOP
!
! !ROUTINE: PDAF_set_offline_mode --- Set offline mode of PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_set_offline_mode(screen)

! !DESCRIPTION:
! Helper routine for PDAF.
!
! This routine allows to activate the offline
! mode of PDAF. Thus, the functionality of
! PDAF to integrate an emsemble will be 
! deactivated.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
! 2023-08 - Lars Nerger - Initial code
! Other revisions - see repository log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: offline_mode, subtype_filter
  USE PDAF_utils_filters, &
       ONLY: PDAF_configinfo_filters

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(in) :: screen
!EOP

! *** Set offline mode ***

  offline_mode = .true.

  IF (screen > 0) THEN
     WRITE (*,'(/a,4x,a)') 'PDAF','Activate PDAF offline mode'

     ! Print configuration info
     CALL PDAF_configinfo_filters(subtype_filter, 1)
  END IF

END SUBROUTINE PDAF_set_offline_mode
