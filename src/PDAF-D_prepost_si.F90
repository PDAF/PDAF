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
!$Id: PDAF-D_prepost_si.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_prepost_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_prepost_si(outflag)

! !DESCRIPTION:
! Interface routine called from the model during the 
! forecast of each ensemble state to transfer data
! from the model to PDAF and to perform the 
! pre/poststep. 
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using pre-defined routine names.
!
! Variant for domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector
       distribute_state_pdaf, &     ! Routine to distribute a state vector
       prepoststep_pdaf, &          ! User supplied pre/poststep routine
       next_observation_pdaf        ! Routine to provide time step, time and dimension
                                    !   of next observation

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_prepost
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_prepost(collect_state_pdaf, distribute_state_pdaf, &
       prepoststep_pdaf, next_observation_pdaf, outflag)

END SUBROUTINE PDAF_prepost_si
