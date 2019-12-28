!$Id$
!>  Set full dimension of observations
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to determine the dimension of the 
!! observation vector according to the current 
!! time step for all observations required for 
!! the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! \date Later revisions - see repository log
!!
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_dim_obs_f_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f  !< Dimension of full observation vector


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_dim_obs_f_pdafomi(step, dim_obs_f)

END SUBROUTINE init_dim_obs_f_pdaf

