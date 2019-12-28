!$Id$
!>  Get local mean observation error variance
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! This routine will only be called, if the 
!! local adaptive forgetting factor feature 
!! is used. Please note that this is an 
!! experimental feature.
!!
!! The routine is called in the loop over all
!! local analysis domains during each analysis
!! by the routine PDAF_set_forget_local that 
!! estimates a local adaptive forgetting factor.
!! The routine has to initialize the mean observation 
!! error variance for the current local analysis 
!! domain.  (See init_obsvar() for a global variant.)
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! \date Later revisions - see repository log
!!
SUBROUTINE init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, meanvar_l)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_obsvar_l_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)

END SUBROUTINE init_obsvar_l_pdaf
