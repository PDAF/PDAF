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
! !ROUTINE: PDAF_reset_forget --- Manually reset forgetting factor
!
! !INTERFACE:
SUBROUTINE PDAF_reset_forget(forget_in)

! !DESCRIPTION:
! Helper routine for PDAF.
! The routine allows to manually set the forgetting
! factor to a new value. Usually this should be called
! in assimilate_pdaf before calling the analysis step
! routine.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: forget

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(in) :: forget_in    ! New value of forgetting factor
!EOP

! *** Set forgetting factor ***

  forget = forget_in

END SUBROUTINE PDAF_reset_forget
