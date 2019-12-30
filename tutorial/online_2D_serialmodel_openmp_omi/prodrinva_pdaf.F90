!$Id$
!>  Compute product of inverse of R with some matrix
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEEK/SEIK/ETKF/ESTKF
!!
!! The routine is called during the analysis step.
!! It has to compute the product of the inverse of 
!! the observation error covariance matrix with
!! the matrix of observed ensemble perturbations
!! (SEIK/ETKF/ESTKF).
!!
!! Implementation for the 2D online example.
!!
!! __Revision history:__
!!  2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)

  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p           !< PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank                !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    !< PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p,rank) !< Input matrix from PDAF analysis routine
  REAL, INTENT(out)   :: C_p(dim_obs_p,rank) !< Output matrix

! *** local variables ***
  INTEGER :: i, j       ! index of observation component
  REAL :: ivariance_obs ! inverse of variance of the observations


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** initialize numbers
  ivariance_obs = 1.0 / rms_obs ** 2


! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

  DO j = 1, rank
     DO i = 1, dim_obs_p
        C_p(i, j) = ivariance_obs * A_p(i, j)
     END DO
  END DO

END SUBROUTINE prodRinvA_pdaf
