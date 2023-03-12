C     $Id$
#include "CPP_OPTIONS.h"
C!BOP
C
C     !ROUTINE: prodRinvA_pdaf --- Compute product of inverse of R with some matrix
C
C     !INTERFACE:
      SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank,
     &      obs_p, A_p, C_p)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: SEEK/SEIK/ETKF/ESTKF
C
C     The routine is called during the analysis step.
C     It has to compute the product of the inverse of 
C     the observation error covariance matrix with
C     the matrix of observed EOF modes (SEEK) or 
C     observed ensemble perturbations (SEIK/ETKF/ESTKF).
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &      ONLY: rms_obs

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step                ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
      INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
      _RL, INTENT(in)     :: obs_p(dim_obs_p)    ! PE-local vector of observations
      _RL, INTENT(in)     :: A_p(dim_obs_p,rank) ! Input matrix from SEEK_ANALYSIS
      _RL, INTENT(out)    :: C_p(dim_obs_p,rank) ! Output matrix

C     !CALLING SEQUENCE:
C     Called by: PDAF_etkf_analysis        (as U_prodRinvA)
C     Called by: PDAF_estkf_analysis       (as U_prodRinvA)
C!EOP

C     *** local variables ***
      INTEGER :: i, j       ! index of observation component


C     **********************
C     *** INITIALIZATION ***
C     **********************
  

C     *************************************
C     ***                -1             ***
C     ***           C = R   A           ***
C     ***                               ***
C     *** The inverse observation error ***
C     *** covariance matrix is not      ***
C     *** computed explicitely.         ***
C     *************************************

      C_p = A_p/rms_obs


      END SUBROUTINE prodRinvA_pdaf
