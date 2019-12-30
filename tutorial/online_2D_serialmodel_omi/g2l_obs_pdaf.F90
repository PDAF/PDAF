!$Id$
!>  Restrict an observation vector to local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the analysis step
!! on each of the local analysis domains.
!! It has to restrict the full vector of all 
!! observations required for the loop of localized 
!! analyses on the PE-local domain to the current 
!! local analysis domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code for PDAF-OMI
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: g2l_obs_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Current local analysis domain
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local obs. vector
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(in)    :: ostate_f(dim_obs_f)   !< Full PE-local obs. vector
  REAL, INTENT(out)   :: ostate_l(dim_obs_l)   !< Obs. vector on local domain


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)

END SUBROUTINE g2l_obs_pdaf
