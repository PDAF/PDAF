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
!> Module providing timing functions
!!
!! This module provides methods to perform timings of 
!! parts of a program execution. It uses the MPI
!! function MPI_Wtime.
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2000-11 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_timer

  IMPLICIT NONE
  SAVE
  
  PUBLIC :: PDAF_timeit, PDAF_time_tot, PDAF_time_temp

  PRIVATE
  REAL, ALLOCATABLE :: t_start(:), t_end(:)
  REAL, ALLOCATABLE :: t_total(:), t_temp(:)

CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!> Initialize Counters and time regions
!!
!! Subroutine to initialize counters and to perform timing of a region
!! specified by timerID.
!! Usage:
!!  * CALL PDAF\_timeit(N,'ini') - Allocates and initializes N counters
!!  * CALL PDAF\_timeit(M,'new') - Start timing region for counter M
!!  * CALL PDAF\_timeit(M,'old') - End timing region for counter M
!!  * CALL PDAF\_timeit(M,'fin') - Finalized and deallocates all counters
!!
  SUBROUTINE PDAF_timeit(timerID, operation)

    USE mpi

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)          :: timerID   !< ID of timer
    CHARACTER(len=3), INTENT(in) :: operation !< Requested operation 

!$OMP MASTER
    ! Initialize timers
    IF (operation == 'ini') THEN
       IF ( .NOT. (ALLOCATED(t_start))) THEN
          ALLOCATE(t_start(timerID), t_end(timerID))
          ALLOCATE(t_total(timerID), t_temp(timerID))
       END IF
        
       t_total = 0.0
    END IF
    
    ! Begin timing region
    IF (operation == 'new') THEN
       t_start(timerID) = MPI_Wtime()
    END IF

    ! End timing region
    IF (operation == 'old') THEN
       t_end(timerID) = MPI_Wtime()
       t_temp(timerID) = t_end(timerID) - t_start(timerID)
       t_total(timerID) = t_total(timerID) + t_temp(timerID)
    END IF
    
    ! Finalize timers
    IF (operation == 'fin') THEN
       DEALLOCATE(t_start, t_end)
       DEALLOCATE(t_total, t_temp)
    END IF
!$OMP END MASTER
    
  END SUBROUTINE PDAF_timeit

!-------------------------------------------------------------------------------
!> Read out timers for last timing interval
!!
!! Read out the value of the timer in seconds for the last 
!! passage of the timing region defined by timerID.
!!
  REAL FUNCTION PDAF_time_temp(timerID)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: timerID             !< ID of timer

    PDAF_time_temp = t_temp(timerID)

  END FUNCTION PDAF_time_temp

!-------------------------------------------------------------------------------
!> Read out total time of a timing region.
!!
!! Read out the accumulated value of the timer in seconds
!! for the timing region define by timerID.
!!
    REAL FUNCTION PDAF_time_tot(timerID)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: timerID             !< ID of timer

    PDAF_time_tot = t_total(timerID)

  END FUNCTION PDAF_time_tot

END MODULE PDAF_timer
