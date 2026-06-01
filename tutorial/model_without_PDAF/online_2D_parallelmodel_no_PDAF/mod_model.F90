!$Id: mod_model.F90 332 2019-12-30 09:37:03Z lnerger $
!> Module for model-specific variables
!!
!! This module provides variables needed for the 
!! 2-dimensional tutorial model without parallelization.
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_model

  IMPLICIT NONE
  SAVE

! *** Variables specific for 2D tutorial model ***

  INTEGER :: nx                     !< Size of 2D grid in x-direction
  INTEGER :: ny                     !< Size of 2D grid in y-direction
  INTEGER :: total_steps            !< Total number of time steps
  REAL, ALLOCATABLE :: field_p(:,:) !< Process-local part of model field

  INTEGER :: nx_p                 !< Process-local size in x-direction

END MODULE mod_model
