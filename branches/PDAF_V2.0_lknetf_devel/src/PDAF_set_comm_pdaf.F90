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
! !ROUTINE: PDAF_set_comm_pdaf --- Set MPI communicator for PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_set_comm_pdaf(in_COMM_pdaf)

! !DESCRIPTION:
! Helper routine for PDAF.
! This routine allows to set the overall (world)
! MPI communicator for PDAF. By default this is 
! MPI_COMM_WORLD. However, in the case that not
! all processes call PDAF or participate in the DA
! and forecasting - as for example in case of
! an IO server that uses separate MPI tasks - a
! separate communicator can be set.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filtermpi, &
       ONLY: isset_comm_pdaf, COMM_pdaf

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(in) :: in_COMM_pdaf    ! MPI communicator for PDAF
!EOP

! *** Set ensemble member ***

  COMM_pdaf = in_COMM_pdaf

  isset_comm_pdaf = .true.

END SUBROUTINE PDAF_set_comm_pdaf
