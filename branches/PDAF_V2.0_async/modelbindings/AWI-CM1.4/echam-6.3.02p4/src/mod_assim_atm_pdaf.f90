!$Id$
!>  Module holding variables for assimilation
!!
!! Module providing variables needed for the 
!! assimilation to be shared within the different
!! user-supplied routines for PDAF.
!! 
!! This module holds the variables that are specific
!! for the atmosphere-component of AWI-CM.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
MODULE mod_assim_atm_pdaf

  IMPLICIT NONE
  SAVE

! *** This module holds the variables specific for observations ***
! *** in the atmosphere. Their values are set in init_pdaf.     ***

! Settings for observations - available as namelist read-in
  INTEGER :: delt_obs_atm        ! time step interval between assimilation steps - Atmosphere
  INTEGER :: delt_obs_atm_offset ! time step offset until first analysis step

! Settings for precisions for ECHAM
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: wp = dp   

END MODULE mod_assim_atm_pdaf
