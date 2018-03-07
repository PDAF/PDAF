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
!$Id: PDAF-D_print_version.F90 1872 2017-12-22 09:38:03Z lnerger $
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
     WRITE(*, '(/a)') 'PDAF    ++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++                      PDAF                      +++'
     WRITE(*, '(a)')  'PDAF    +++      Parallel Data Assimilation Framework      +++'
     WRITE(*, '(a)')  'PDAF    +++                                                +++'     
     WRITE(*, '(a)')  'PDAF    +++                Version 1.13.1                  +++' 
     WRITE(*, '(a/)') 'PDAF    ++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  END IF

END SUBROUTINE PDAF_print_version
