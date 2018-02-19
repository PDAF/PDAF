! Copyright (c) 2004-2018 Lars Nerger and Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_Btau --- Scalar to regulate the strength of the nudging.
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_Btau(t_now, t_last, t_obs, nudge, flag)


! !DESCRIPTION:
! Helper routine for proposal filter to regulate the strenght of the nudging
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-10 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!

  USE  PDAF_mod_ewpf, &
       ONLY: bt, &            ! Real to reduce the strenght of the nudging
       start_nudging          ! Relative timestep when the nudging should start

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: t_now     ! Current model timestep
  INTEGER, INTENT(in) :: t_last    ! Time step of previous assimilation
  INTEGER, INTENT(in) :: t_obs     ! Timestep of next observation
  REAL, INTENT(out) :: nudge       ! Forcing strength
  INTEGER, INTENT (in) :: flag     ! Flag to choose between different forcings
 
! Local variables
  REAL :: t_rel
  REAL :: a, b, c

! Options for the nudging:
! Flag    Forcingfunction
!  0        No forcing
!  1        Linear forcing starting previous_step till next_observations
!  2        No forcing till timestep= start_nudging, linear afterwards
!  3        No forcing till timestep=start_nudging, exponential afterwars

  t_rel = real( t_now - t_last) / real( t_obs - t_last )

  IF (flag == 0) THEN
     nudge = 0
  ELSEIF (flag == 1) THEN
     nudge = bt * t_rel
  ELSEIF (flag == 2) THEN
     IF ( t_rel <= start_nudging) THEN
        nudge = 0
     ELSE 
        nudge = bt* real( t_rel - start_nudging) / real( 1.0 - start_nudging )
     ENDIF
  ELSEIF (flag == 3) THEN
     IF ( t_rel <= start_nudging) THEN
        nudge = 0
     ELSE
        a = 10**(-3)
        b = LOG(real(10**3))
        c = real( t_rel - start_nudging) / real( 1.0 - start_nudging )
        nudge = bt*a*exp(b*c)
     ENDIF
  ENDIF

END SUBROUTINE PDAF_ewpf_Btau
