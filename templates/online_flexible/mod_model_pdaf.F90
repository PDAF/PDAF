!> Interface module between the model and PDAF
!!
!! This module includes the model module(s) to access
!! the model variables and make them accessible for
!! the PDAF user routines in a uniform way.
!!
!! __Revision history:__
!! * 2026-04 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_model_pdaf

  ! Include model variables (we include here all without 'only')
  USE mod_model

  IMPLICIT NONE
  SAVE

END MODULE mod_model_pdaf
