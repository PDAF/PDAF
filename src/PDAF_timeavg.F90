! Copyright (c) 2004-2020 Lars Nerger
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
! The code is very generic. Basically the only
! filter-specific part if the call to the
! update-routine PDAF\_X\_update where the analysis
! is computed.  The filter-specific subroutines that
! are specified in the call to PDAF\_put\_state\_X
! are passed through to the update routine
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF-D_put_state_lestkf.F90 374 2020-02-26 12:49:56Z lnerger $
!BOP
!
! !ROUTINE: PDAF_timeavg --- Routine to increment state vector for time averaging
!
! !INTERFACE:
SUBROUTINE PDAF_timeavg(U_collect_state, avgsteps)

! !DESCRIPTION:
! Interface routine called from the model after the 
! forecast of each ensemble state to transfer data
! from the model to PDAF. For the parallelization 
! this involves transfer from model PEs to filter 
! PEs.\\
! During the forecast phase state vectors are 
! re-initialized from the forecast model fields
! by U\_collect\_state. 
! At the end of a forecast phase (i.e. when all 
! ensemble members have been integrated by the model)
! sub-ensembles are gathered from the model tasks.
! Subsequently the filter update is performed.
!
! Variant for LESTKF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filter, &
       ONLY: dim_p, nsteps, state, eofV, subtype_filter, &
       member, member_save, ensAvg
  USE PDAF_mod_filtermpi, &
       ONLY: modelpe, mype_model

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: avgsteps  ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state       ! Routine to collect a state vector

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: U_collect_state
! Calls: PDAF_gather_ens
! Calls: PDAF_lestkf_update
! Calls: PDAF_timeit
!EOP

! local variables
  INTEGER :: i   ! Counter


! **************************************************
! *** Save forecast state back to the ensemble   ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps > 0) THEN

     CALL PDAF_timeit(41, 'new')

     avgsteps = avgsteps + 1

     modelpes: IF (modelpe) THEN
        IF (mype_model == 0) &
             write (*,*) 'PDAF  timeavg: Inrement member', member, 'avgsteps', avgsteps
        
        ! Store member index for PDAF_get_memberid
        member_save = member

        IF (subtype_filter /= 2 .AND. subtype_filter /= 3) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, eofV(1 : dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1 : dim_p))
        END IF
        
        ! Increment array of time-averaged ensemble
        ensAvg(1: dim_p, member) = ensAvg(1: dim_p, member) + eofV(1 : dim_p, member)

     END IF modelpes

     CALL PDAF_timeit(41, 'old')
  END IF doevol

! ********************
! *** finishing up ***
! ********************

END SUBROUTINE PDAF_timeavg
