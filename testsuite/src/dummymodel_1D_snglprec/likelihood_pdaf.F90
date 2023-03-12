!$Id: likelihood_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: likelihood --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood_pdaf(step, dim_obs_p, obs_p, resid, likely)

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
! This routine is called by all filter processes.
!
! Implementation for the dummy model with domain
! decomposition. Here, we assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code based in prodRinvA_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs
  USE mod_parallel, &
       ONLY: npes_filter, COMM_filter, MPI_SUM, MPI_REAL, MPIerr

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
  REAL, INTENT(in)    :: resid(dim_obs_p)    ! Input vector of residuum
  REAL, INTENT(out)   :: likely              ! Output vector - log likelihood

! !CALLING SEQUENCE:
! Called by: PDAF_netf_analysis        (as U_likelihood)
!EOP

! *** local variables ***
  INTEGER :: i          ! index of observation component
  REAL :: ivariance_obs ! inverse of variance of the observations
  REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual
  REAL :: likely_p      ! PE-local part of likelihood


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

  ALLOCATE(Rinvresid(dim_obs_p))

  DO i = 1, dim_obs_p
     Rinvresid(i) = ivariance_obs * resid(i)
  END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

  ! Gaussian errors: compute exp(-0.5*resid^T*R^-1*resid)
  ! With MPI-parallelization this a only for the process-local part 
  CALL sgemv('t', dim_obs_p, 1, 0.5, resid, &
       dim_obs_p, Rinvresid, 1, 0.0, likely_p, 1)

  IF (npes_filter==1) THEN
     likely = EXP(-likely_p)
  ELSE
     ! With MPI-parallelization we need to sum over the partial values
     CALL MPI_allreduce(likely_p, likely, 1, &
       MPI_REAL, MPI_SUM, COMM_filter, MPIerr)

     likely = EXP(-likely)
  END IF


! *** Clean up ***

  DEALLOCATE(Rinvresid)

END SUBROUTINE likelihood_pdaf
