!$Id$
!>  Get mean observation error variance
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!!
!! This routine will only be called, if the adaptive
!! forgetting factor feature is used. Please note that
!! this is an experimental feature.
!!
!! The routine is called in global filters (like ESTKF)
!! during the analysis or in local filters (e.g. LESTKF)
!! before the loop over local analysis domains 
!! by the routine PDAF_set_forget that estimates an 
!! adaptive forgetting factor.  The routine has to 
!! initialize the mean observation error variance.  
!! For global filters this should be the global mean,
!! while for local filters it should be the mean for the
!! PE-local  sub-domain.  (See init_obsvar_l_pdaf()
!! for a localized variant for local filters.)
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! * Later revisions - see repository log
!!
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_obsvar_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
  REAL, INTENT(out)   :: meanvar       !< Mean observation error variance


! *****************************
! *** Compute mean variance ***
! *****************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)

END SUBROUTINE init_obsvar_pdaf
