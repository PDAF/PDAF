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
!$Id: PDAF-D_prepost.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_prepost --- Interface to PDAF for calling prepoststep
!
! !INTERFACE:
SUBROUTINE PDAF_prepost(U_collect_state, U_distribute_state, &
     U_prepoststep, U_next_observation, outflag)

! !DESCRIPTION:
! Interface routine called from the model at each time
! step during the forecast of each ensemble state. If
! the time of the next analysis step is reached the
! forecast state is transferred to PDAF and the 
! user-supplied routine U_prepoststep is executed
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
  USE PDAF_mod_filter, &
       ONLY: cnt_steps, nsteps
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world


  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(out) :: outflag  ! Status flag
  
! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
       U_prepoststep, &           ! User supplied pre/poststep routine
       U_next_observation, &      ! Routine to provide time step, time and dimension
                                  !   of next observation
       U_distribute_state         ! Routine to distribute a state vector

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_put_state_prepost
!EOP

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Increment time step counter
  cnt_steps = cnt_steps + 1


! ********************************
! *** At end of forecast phase ***
! ********************************

  IF (cnt_steps == nsteps) THEN

     IF (mype_world==0) WRITE(*,'(a,5x,a)') 'PDAF','Perform pre/poststep of PDAF'

     ! *** Call analysis step ***

     CALL PDAF_put_state_prepost(U_collect_state, U_prepoststep, outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
             U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     outflag = 0
  END IF

END SUBROUTINE PDAF_prepost
