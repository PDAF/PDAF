!$Id: mod_assim_pdaf.f90 2179 2020-03-18 18:48:06Z lnerger $
!>  Module holding variables for assimilation
!!
!! Module providing variables needed for the 
!! assimilation to be shared within the different
!! user-supplied routines for PDAF.
!! 
!! This module holds the variables that are specific
!! for the ocean-component (FESOM) of AWI-CM.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
MODULE mod_assim_oce_pdaf

  IMPLICIT NONE
  SAVE

! *** This module holds the variables specific for observations ***
! *** in the ocean. Their values are set in init_pdaf.          ***

! Settings for observations - available as namelist read-in
  INTEGER :: delt_obs_ocn        ! time step interval between assimilation steps - Ocean
  INTEGER :: delt_obs_ocn_offset ! time step offset until first analysis step

! File output and input - available as as namelist read-in
  CHARACTER(len=100) :: path_obs_rawprof  = ''      ! Path to raw profile observation files
  CHARACTER(len=110) :: file_rawprof_prefix  = ''   ! file name prefix for profile observations 
  CHARACTER(len=110) :: file_rawprof_suffix  = '.nc'! file name suffix for profile observations 
  CHARACTER(len=110) :: file_syntobs = 'syntobs.nc' ! File name for synthetic observations

END MODULE mod_assim_oce_pdaf
