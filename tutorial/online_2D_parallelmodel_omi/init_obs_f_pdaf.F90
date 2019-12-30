!$Id$
!>  Initialize observation vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called in PDAF_X_update
!! before the loop over all local analysis domains
!! is entered. It has to provide the full observation 
!! vector according to current time step (where 'full' 
!! means 'all observations required for the localized 
!! analysis on the PE-local domain').  This routine 
!! is only used if a globally adaptive 
!! forgetting factor is requested, rather than an 
!! individual forgetting factor for each analysis 
!! domain. This routine has to be implemented 
!! consistently with the routines for the full 
!! observation dimension and the full observation 
!! operator. The forgetting factor will only be 
!! globally adaptive, if the full observation vector 
!! is the global observation vector.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! * Later revisions - see repository log
!!
SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_obs_f_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_obs_f_pdafomi(step, dim_obs_f, observation_f)

END SUBROUTINE init_obs_f_pdaf

