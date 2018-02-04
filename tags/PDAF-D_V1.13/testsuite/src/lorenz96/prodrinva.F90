!$Id: prodrinva.F90 814 2010-01-27 11:30:16Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA(step, dim_obs, rank, obs, A, C)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK, SEIK):
!
! The routine is called during the analysis step.
! It has to compute the product of the inverse of 
! the observation error covariance matrix with
! the matrix of observed EOF modes (SEEK) or 
! observed ensemble perturbations (SEIK).
!
! This variant is for the Lorenz96 model without
! parallelization. We assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step             ! Current time step
  INTEGER, INTENT(in) :: dim_obs          ! PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank             ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs(dim_obs)     ! PE-local vector of observations
  REAL, INTENT(in)    :: A(dim_obs, rank) ! Input matrix from SEEK_ANALYSIS
  REAL, INTENT(out)   :: C(dim_obs, rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis_newT   (as U_prodRinvA)
!EOP

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
     DO i = 1, dim_obs
        C(i, j) = ivariance_obs * A(i, j)
     END DO
  END DO

END SUBROUTINE prodRinvA
