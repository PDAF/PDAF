!$Id$
!>  Initialize local observation vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the analysis step
!! on each local analysis domain in 
!! PDAF_X_analysis.  It has to initialize 
!! the local vector of observations for the 
!! current local analysis domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! \date Later revisions - see repository log
!!
SUBROUTINE init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_obs_l_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

END SUBROUTINE init_obs_l_pdaf

