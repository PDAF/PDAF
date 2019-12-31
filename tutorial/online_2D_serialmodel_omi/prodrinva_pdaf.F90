!$Id$
!>  Compute product of inverse of R with some matrix
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEEK/SEIK/ETKF/ESTKF
!!
!! The routine is called during the analysis step
!! of a global square-root filter.
!! It has to compute the product of the inverse of 
!! the observation error covariance matrix with
!! the matrix of observed ensemble perturbations.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code for PDAF-OMI
!! * Later revisions - see repository log
!!
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: prodRinvA_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p           !< PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank                !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    !< PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p,rank) !< Input matrix from PDAF analysis routine
  REAL, INTENT(out)   :: C_p(dim_obs_p,rank) !< Output matrix


! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL prodRinvA_pdafomi(step, dim_obs_p, rank, obs_p, A_p, C_p)

END SUBROUTINE prodRinvA_pdaf
