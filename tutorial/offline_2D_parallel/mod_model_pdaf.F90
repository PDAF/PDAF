!> Module to declare model variables for offline coupling
!!
!! This module declares model-related variables for offline
!! coupled DA. While for online coupling the variables
!! would be included by 'use' statements from model modules,
!! we have to declare the variables explicitly of the 
!! offline coupling.
!!
!! Implementation for the 2D offline example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code extracting from mod_assimilation
!! * Later revisions - see repository log
!!
MODULE mod_model_pdaf

  IMPLICIT NONE
  SAVE

! *** Variables specific for offline tutorial example ***

  INTEGER :: nx, ny                     !< Size of 2D grid
  INTEGER, ALLOCATABLE :: local_dims(:) !< Array for local state dimensions

END MODULE mod_model_pdaf
