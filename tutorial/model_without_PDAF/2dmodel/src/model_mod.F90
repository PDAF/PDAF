!> Module for 2D tutorial model
!!
!! This module provides variables for the 
!! 2-dimensional tutorial model with parallelization.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_mod

  implicit none
  save
  public

! *** Variables specific for 2D tutorial model ***

  integer :: nx                      !< Size of 2D grid in x-direction
  integer :: ny                      !< Size of 2D grid in y-direction
  integer :: n_dim                   !< Number of model dimensions
  integer :: total_steps             !< Total number of time steps
  real, allocatable :: fieldA(:,:)   !< Model field A
  real, allocatable :: fieldB(:,:)   !< Model field B
  real, allocatable :: coords_x(:)   !< Coordinates in x-direction
  real, allocatable :: coords_y(:)   !< Coordinates in y-direction

end module model_mod
