!$Id: likelihood.F90 1673 2016-12-02 10:46:31Z lnerger $
!BOP
!
! !ROUTINE: likelihood --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood(step, dim_obs, obs, resid, likely)

! !DESCRIPTION:
! User-supplied routine for PDAF (NETF):
!
! The routine is called during the analysis step.
! It has to compute the likelihood of the
! ensemble according to the difference from the
! observation (residual) and the error distribution
! of the observations.
!
! In general this routine is similar to the routine
! prodRinvA used for ensemble square root Kalman
! filters. As an addition to this routine, we here have
! to evaluate the likelihood weight according the
! assumed observation error statistics.
!
! This variant is for the Lorenz96 model without
! parallelization. We assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs, obs_err_type

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step             ! Current time step
  INTEGER, INTENT(in) :: dim_obs          ! PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs(dim_obs)     ! PE-local vector of observations
  REAL, INTENT(in)    :: resid(dim_obs)   ! Input vector of residuum
  REAL, INTENT(out)   :: likely           ! Output vector - log likelihood

! !CALLING SEQUENCE:
! Called by: PDAF_netf_analysis        (as U_likelihood)
!EOP

! *** local variables ***
  INTEGER :: i          ! index of observation component
  REAL :: ivariance_obs ! inverse of variance of the observations
  REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** initialize numbers
  ivariance_obs = 1.0 / rms_obs ** 2


! ***************************************
! *** Before computing the likelihood ***
! *** scale by observation error      ***
! ***                   -1            ***
! ***      Rinvresid =  R  resid      ***
! ***                                 ***
! *** The inverse observation error   ***
! *** covariance matrix is not        ***
! *** computed explicitely.           ***
! ***************************************

  ALLOCATE(Rinvresid(dim_obs))

  DO i = 1, dim_obs
     Rinvresid(i) = ivariance_obs * resid(i)
  END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

  IF (obs_err_type==0) THEN

     ! Gaussian errors
     ! Calculate exp(-0.5*resid^T*R^-1*resid)
     CALL dgemv('t', dim_obs, 1, 0.5, resid, &
          dim_obs, Rinvresid, 1, 0.0, likely, 1)
     likely = EXP(-likely)

  ELSE

     ! Double-exponential errors
     ! Calculate exp(-SUM(ABS(resid)))
     likely = 0.0
     DO i = 1, dim_obs
        likely = likely + ABS(Rinvresid(i))
     END DO
     likely = EXP(-likely)

  END IF

! *** Clean up ***

  DEALLOCATE(Rinvresid)

END SUBROUTINE likelihood
