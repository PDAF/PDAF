!$Id: main.F90 1410 2013-09-25 12:04:33Z lnerger $
!BOP
!
! !ROUTINE: main --- Main driver for PDAF testsuite
!
! !INTERFACE:
PROGRAM MAIN

! !DESCRIPTION:
! This is a simple model program to demonstrate the
! fully-parallel implementation of the online mode of PDAF. 
!
! The simple model has a 2-dimensional mesh. The initial state
! is read from a file. The time stepping consists in shifting
! the field vertically (in the direction of the first array index)
! by one grid point per time step. A period boundary condition is
! applied by inserting the field from the upper boundary into the
! lower one. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code based on dummy model example
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
!EOP


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  WRITE (*, '(/17x, a/)') '+++++ PDAF tutorial - online mode +++++'
  WRITE (*, '(16x, a)') 'Tutorial: 2D model without parallelization'
  WRITE (*, '(/)')
     
  ! Initialize model
  CALL initialize()


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration
  CALL integrate()

END PROGRAM MAIN
