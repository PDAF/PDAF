!$Id: timer.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !MODULE:
MODULE timer

! !DESCRIPTION: 
! This module provides methods to perform timings of 
! parts of a program execution. It uses the intrinsic 
! function SYSTEM\_CLOCK.
!
! !REVISION HISTORY:
! 2000-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
  
  PUBLIC :: timeit, time_tot, time_temp
!EOP

  PRIVATE
  INTEGER :: t_rate
  INTEGER, ALLOCATABLE :: t_start(:), t_end(:)
  REAL, ALLOCATABLE    :: t_total(:), t_temp(:)

CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: timeit - Initialize Counters and time regions
!
! !INTERFACE: timeit()
  SUBROUTINE timeit(timerID, operation)

! !DESCRIPTION:
! Subroutine to initialize counters and to perform timing of a region
! specified by timerID.
! Usage:\\
!   CALL PDAF\_timeit(N,'ini') - Allocates and initializes N counters\\
!   CALL PDAF\_timeit(M,'new') - Start timing region for counter M\\
!   CALL PDAF\_timeit(M,'old') - End timing region for counter M\\
!   CALL PDAF\_timeit(M,'fin') - Finalized and deallocates all counters\\

! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: timerID             ! ID of timer
    CHARACTER(len=3), INTENT(in) :: operation  ! Requested operation 
!EOP

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
       CALL SYSTEM_CLOCK(t_start(timerID))
    END IF

    ! End timing region
    IF (operation == 'old') THEN
       CALL SYSTEM_CLOCK(t_end(timerID), t_rate)
       t_temp(timerID) = REAL(t_end(timerID) - t_start(timerID)) &
            / REAL(t_rate)
       t_total(timerID) = t_total(timerID) + REAL(t_end(timerID) - &
            t_start(timerID)) / REAL(t_rate)
    END IF
    
    ! Finalize timers
    IF (operation == 'fin') THEN
       DEALLOCATE(t_start, t_end)
       DEALLOCATE(t_total, t_temp)
    END IF
    
  END SUBROUTINE timeit

!-------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: time_temp - Read out timers for last timing interval
!
! !INTERFACE: time_temp()
  REAL FUNCTION time_temp(timerID)

! !DESCRIPTION:
! Read out the value of the timer in seconds for the last 
! passage of the timing region defined by timerID.

! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: timerID             ! ID of timer
!EOP

    time_temp = t_temp(timerID)

  END FUNCTION time_temp

!-------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: PDAF_time_tot - Read out total time of a timing region.
!
! !INTERFACE: time_tot()
    REAL FUNCTION time_tot(timerID)

! !DESCRIPTION:
! Read out the accumulated value of the timer in seconds
! for the timing region define by timerID.

! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: timerID             ! ID of timer
!EOP

    time_tot = t_total(timerID)

  END FUNCTION time_tot

END MODULE timer
