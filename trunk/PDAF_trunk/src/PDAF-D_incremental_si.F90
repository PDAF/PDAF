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
!$Id: PDAF-D_incremental_si.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_incremental_si --- Interface for PDAF with incremental updating
!
! !INTERFACE:
SUBROUTINE PDAF_incremental_si(steps)

! !DESCRIPTION:
! Interface routine for PDAF with incremental
! analysis updating. It is called from the model 
! during the forecast. It has to provide the
! user-supplied routine U\_dist\_stateinc with
! the analysis increment from the analysis routine.
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using a pre-defined routine name.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: steps ! Time steps over which increment is distributed

! ! Names of external subroutines 
  EXTERNAL :: distribute_stateinc_pdaf  ! Add state increment during integration
  
! !CALLING SEQUENCE:
! Called by: model code
! Calls: PDAF_incremental
!EOP


! *******************************************************
! *** Call the full routine for incremental updating  ***
! *******************************************************

  CALL PDAF_incremental(steps, distribute_stateinc_pdaf)

END SUBROUTINE PDAF_incremental_si
