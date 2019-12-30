!$Id$
!>  Main driver for PDAF tutorial (without assimilation)
!!
!! This is a simple model program to demonstrate the
!! fully-parallel implementation of the online mode of PDAF. 
!!
!! The simple model has a 2-dimensional grid. The initial state
!! is read from a file. The time stepping consists in shifting
!! the field vertically (in the direction of the first array index)
!! by one grid point per time step. A period boundary condition is
!! applied by inserting the field from the upper boundary into the
!! lower one. 
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code based on dummy model example
!! * Later revisions - see repository log
!!
PROGRAM MAIN

  IMPLICIT NONE

! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! *** Initial Screen output ***
  WRITE (*, '(/17x, a/)') '+++++ PDAF tutorial - online mode +++++'
  WRITE (*, '(16x, a)') 'Tutorial: 2D model without parallelization'
  WRITE (*, '(/)')
     
  ! *** Initialize model ***
  CALL initialize()


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration ***
  CALL integrate()

END PROGRAM MAIN
