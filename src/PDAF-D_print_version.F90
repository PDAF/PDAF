! Copyright (c) 2004-2018 Lars Nerger
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
! !ROUTINE: PDAF_print_version --- Display version information for PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_print_version()

! !DESCRIPTION:
! This routine displays the information from PDAF.
! Possible are to display the timing information and
! allocated memory.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filtermpi, &
       ONLY: MPI_COMM_WORLD, mype_world, MPIerr

  IMPLICIT NONE


  ! Determine parallel rank of process
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)


! *********************************
! *** Print version information ***
! *********************************

  IF (mype_world==0) THEN
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++                      PDAF                       +++'
     WRITE(*, '(a)')  'PDAF    +++      Parallel Data Assimilation Framework       +++'
     WRITE(*, '(a)')  'PDAF    +++                                                 +++'     
     WRITE(*, '(a)')  'PDAF    +++            Version 1.13.2 - EGU 2019            +++' 
     WRITE(*, '(a)')  'PDAF    +++                                                 +++'     
     WRITE(*, '(a)')  'PDAF    +++ Note: This is a reduced version prepared for    +++'
     WRITE(*, '(a)')  'PDAF    +++ the short course on practical data assimilation +++'
     WRITE(*, '(a)')  'PDAF    +++ with PDAF given at the EGU Assembly 2019. For   +++'
     WRITE(*, '(a)')  'PDAF    +++ further use of PDAF, we recommend to download   +++'
     WRITE(*, '(a)')  'PDAF    +++ the full PDAF package at http://pdaf.awi.de     +++'
     WRITE(*, '(a)')  'PDAF    +++                                                 +++'     
     WRITE(*, '(a/)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  END IF

END SUBROUTINE PDAF_print_version
