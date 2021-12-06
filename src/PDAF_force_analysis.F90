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
! !ROUTINE: PDAF_force_analysis --- Set ensemble index to force an analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_force_analysis()

! !DESCRIPTION:
! Helper routine for PDAF.
! The routine overwrite member index of the ensemble 
! state by local_dim_ens and the counter cnt_steps
! by nsteps-1. This forces that the analysis
! step is executed at the next call to PDAF_put_state
! or PDAF_assimilate.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: member, local_dim_ens, nsteps, cnt_steps, step_obs, &
       step, screen, use_PDAF_assim
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world

  IMPLICIT NONE
!EOP

! *** Set ensemble member ***

  member = local_dim_ens

  nsteps = cnt_steps + 1

  ! Only in case of using PDAF_assimilate, we need to reset the step counting
  IF (use_PDAF_assim) step_obs = step + nsteps - 1

  IF (screen>0 .AND. mype_world==0) THEN
     WRITE (*,'(a,5x,a,i8)') 'PDAF','!! Force analysis at step', step_obs
  END IF

END SUBROUTINE PDAF_force_analysis
