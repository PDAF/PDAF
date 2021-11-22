!$Id: main.F90 332 2019-12-30 09:37:03Z lnerger $
!>  Main driver for PDAF tutorial
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

  USE mod_parallel_model, &      ! Model parallelization variables
       ONLY: mype_world, init_parallel, finalize_parallel

  IMPLICIT NONE

! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! Initialize parallelization ***
  CALL init_parallel()

  ! *** Initial Screen output ***
  IF (mype_world==0) THEN
     WRITE (*, '(/17x, a/)') '+++++ PDAF tutorial - online mode +++++'
     WRITE (*, '(17x, a)') 'Tutorial: 2D model with parallelization'
     WRITE (*, '(/)')
  END IF

  ! Initialize model ***
  CALL initialize()


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration ***
  CALL integrate()


! **************************
! ***      Clean up      ***
! **************************

  CALL finalize_parallel()

END PROGRAM MAIN
