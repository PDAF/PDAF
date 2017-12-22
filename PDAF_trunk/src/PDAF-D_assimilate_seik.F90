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
!$Id: PDAF-D_assimilate_seik.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_assimilate_seik --- Interface to PDAF for SEIK
!
! !INTERFACE:
SUBROUTINE PDAF_assimilate_seik(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
     U_init_obsvar, U_next_observation, outflag)

! !DESCRIPTION:
! Interface routine called from the model at each time
! step during the forecast of each ensemble state. If
! the time of the next analysis step is reached the
! forecast state is transferred to PDAF and the analysis
! is computed by calling PDAF_put_state_seik. Subsequently, 
! PDAF_get_state is called to initialize the next forecast
! phase. 
!
! The code is very generic. Basically the only
! filter-specific part are the calls to the
! routines PDAF\_put\_state\_X where the analysis
! is computed and PDAF\_get\_state to initialize the next
! forecast phase. The filter-specific call-back subroutines 
! are specified in the calls to the two core routines.
!
! Variant for SEIK with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code
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
  EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
       U_init_dim_obs, &      ! Initialize dimension of observation vector
       U_obs_op, &            ! Observation operator
       U_init_obsvar, &       ! Initialize mean observation error variance
       U_init_obs, &          ! Initialize observation vector
       U_prepoststep, &       ! User supplied pre/poststep routine
       U_prodRinvA, &         ! Provide product R^-1 A
       U_next_observation, &  ! Routine to provide time step, time and dimension
                              !   of next observation
       U_distribute_state     ! Routine to distribute a state vector

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_put_state_seik
! Calls: PDAF_get_state_seik
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

     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

     ! *** Call analysis step ***

     CALL PDAF_put_state_seik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
             U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     outflag = 0
  END IF

END SUBROUTINE PDAF_assimilate_seik
