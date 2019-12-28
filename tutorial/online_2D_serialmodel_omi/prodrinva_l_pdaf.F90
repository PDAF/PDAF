!$Id$
!>  Compute product of inverse of R with some matrix
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the analysis step
!! on each local analysis domain. It has to 
!! compute the product of the inverse of the local
!! observation error covariance matrix with
!! the matrix of locally observed ensemble 
!! perturbations.
!! Next to computing the product, a localizing 
!! weighting can be applied to matrix A to perform
!! observation localization.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! \date Later revisions - see repository log
!!
SUBROUTINE prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  USE interface_pdafomi, &
       ONLY: prodRinvA_l_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p          !< Current local analysis domain
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
  
END SUBROUTINE prodRinvA_l_pdaf
