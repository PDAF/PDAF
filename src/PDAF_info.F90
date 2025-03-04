! Copyright (c) 2004-2025 Lars Nerger
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
!
!> Display information for PDAF
!!
MODULE PDAF_info

CONTAINS
!> Display version information for PDAF
!!
!! This routine displays the information from PDAF.
!! Possible are to display the timing information and
!! allocated memory.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_print_version()

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: COMM_pdaf, mype_world, MPIerr

  IMPLICIT NONE


  ! Determine parallel rank of process
  CALL MPI_Comm_rank(COMM_pdaf, mype_world, MPIerr)


! *********************************
! *** Print version information ***
! *********************************

  IF (mype_world==0) THEN
     WRITE(*, '(/a)') 'PDAF    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++                       PDAF                         +++'
     WRITE(*, '(a)')  'PDAF    +++       Parallel Data Assimilation Framework         +++'
     WRITE(*, '(a)')  'PDAF    +++                                                    +++'
     WRITE(*, '(a)')  'PDAF    +++                   Version 3.0                      +++'
     WRITE(*, '(a)')  'PDAF    +++                                                    +++'
     WRITE(*, '(a)')  'PDAF    +++                   Please cite                      +++'
     WRITE(*, '(a)')  'PDAF    +++ L. Nerger and W. Hiller, Computers and Geosciences +++'
     WRITE(*, '(a)')  'PDAF    +++ 2013, 55, 110-118, doi:10.1016/j.cageo.2012.03.026 +++'
     WRITE(*, '(a)')  'PDAF    +++   when publishing work resulting from using PDAF   +++'
     WRITE(*, '(a)')  'PDAF    +++                                                    +++'
     WRITE(*, '(a)')  'PDAF    +++          PDAF itself can also be cited as          +++'
     WRITE(*, '(a)')  'PDAF    +++  L. Nerger. Parallel Data Assimilation Framework   +++'
     WRITE(*, '(a)')  'PDAF    +++  (PDAF). Zenodo. 2024. doi:10.5281/zenodo.7861812  +++'
     WRITE(*, '(a/)') 'PDAF    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  END IF

END SUBROUTINE PDAF_print_version


!-------------------------------------------------------------------------------
!> Print information for PDAF (timing and memory) to screen
!!
!! This routine displays the information from PDAF.
!! Possible are to display the timing information and
!! allocated memory.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2008-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_print_info(printtype)

  USE PDAF_utils_filters, &
       ONLY: PDAF_print_info_filters

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: printtype       !< Type of screen output:  
                                         !< * (1) general timings
                                         !< * (3) timers focused on call-back routines (recommended)
                                         !< * (4,5) detailed and very detailed timers (to analyze filters)
                                         !< * (10) allocated memory of calling MPI task
                                         !< * (11) globally used memory (needs to e called by all processes)


! ********************************
! *** Print screen information ***
! ********************************

  CALL PDAF_print_info_filters(printtype)

END SUBROUTINE PDAF_print_info

END MODULE PDAF_info
