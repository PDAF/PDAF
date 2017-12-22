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
!$Id: PDAF-D_incremental.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_incremental --- Interface for PDAF with incremental updating
!
! !INTERFACE:
SUBROUTINE PDAF_incremental(steps, U_dist_stateinc)

! !DESCRIPTION:
! Interface routine for PDAF with incremental
! analysis updating. It is called from the model 
! during the forecast. It has to provide the
! user-supplied routine U\_dist\_stateinc with
! the analysis increment from the analysis routine.
!
! !REVISION HISTORY:
! 2006-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: initevol, incremental, dim_p, state_inc
  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: steps ! Time steps over which increment is distributed

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the calling routine)
  EXTERNAL :: U_dist_stateinc  ! Add state increment during integration
  
! !CALLING SEQUENCE:
! Called by: model code
! Calls: U_dist_stateinc
!EOP


! *** INITIALIZATION ***
  IF (mype == 0 .AND. initevol > 0) THEN
     WRITE (*, '(a, 5x, a, i5, a)') &
          'PDAF', '--- Apply incremental analysis update over ', steps, ' time steps'
  END IF

! *** Call user routine to distribute and add increment ***

  IF (incremental == 1) THEN
     CALL U_dist_stateinc(dim_p, state_inc, initevol, REAL(steps))
  END IF

! *** Reset flag INITEVOL ***
  IF (initevol > 0) THEN
     initevol = 0
  END IF

END SUBROUTINE PDAF_incremental
