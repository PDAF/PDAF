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
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
program main

  use model_init_mod, &              ! Model initialization
       only: initialize
  use model_step_mod, &              ! Model integration
       only: stepping
  use model_post_mod, &              ! Model post-processing
       only: postprocess

  implicit none

! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! *** Initial Screen output ***
  write (*, '(/10x, a/)') '+++++ PDAF tutorial - online mode +++++'
  write (*, '(10x, a)')   '        2D model with 2 fields'
  write (*, '(/)')

  ! *** Initialize model ***
  call initialize()


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration ***
  call stepping()


! **************************
! ***      Clean up      ***
! **************************

  ! Clean up model fields
  call postprocess()

end program main
