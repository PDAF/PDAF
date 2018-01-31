C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: prodRinvA_l_pdaf --- Compute product of inverse of R with some matrix
C
C     !INTERFACE:
      SUBROUTINE prodRinvA_l_pdaf(domain_p, step,
     &     dim_obs_l, rank, obs_l, A_l, C_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during the analysis step
C     on each local analysis domain. It has to 
C     compute the product of the inverse of the local
C     observation error covariance matrix with
C     the matrix of locally observed ensemble 
C     perturbations.
C     Next to computing the product,  a localizing 
C     weighting (similar to covariance localization 
C     often used in EnKF) can be applied to matrix A.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: local_range, locweight, srange, rms_obs,
     &     distance_l
      USE mod_parallel_pdaf,
     &      ONLY: mype_filter

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
      INTEGER, INTENT(in) :: step              ! Current time step
      INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
      INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
      _RL, INTENT(in)     :: obs_l(dim_obs_l)     ! Local vector of observations
      _RL, INTENT(inout)  :: A_l(dim_obs_l, rank) ! Input matrix
      _RL, INTENT(out)    :: C_l(dim_obs_l, rank) ! Output matrix

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_analysis   (as U_prodRinvA_l)
C     Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)

      
C     *** local variables ***
      INTEGER :: i, j           ! Index of observation component
      INTEGER :: verbose        ! Verbosity flag
      INTEGER :: verbose_w      ! Verbosity flag for weight computation
      INTEGER, SAVE :: domain_save = -1 ! Save previous domain index
      REAL*8  :: ivariance_obs  ! Inverse of variance of the observations
      INTEGER :: wtype          ! Type of weight function
      INTEGER :: rtype          ! Type of weight regulation
      REAL*8, ALLOCATABLE :: weight(:)   ! Localization weights
      REAL*8, ALLOCATABLE :: A_obs(:,:)  ! Array for a single row of A_l
      REAL*8  :: var_obs        ! Variance of observation error


C     This routine is generic and does not need any change is the
C     observation errors are constant. IF the observation errors
C     vary, one has to use a vector to store the observation errors
C     instead of the single value RMS_OBS used here.

      ivariance_obs = 1.0 / rms_obs**2
      var_obs = rms_obs**2

      
C **********************
C *** INITIALIZATION ***
C **********************

      IF ((domain_p <= domain_save .OR. domain_save < 0)
     &     .AND. mype_filter==0) THEN
         verbose = 1
      ELSE
         verbose = 0
      END IF
      domain_save = domain_p

C     Screen output
      IF (verbose == 1) THEN
         WRITE (*, '(a, 8x, a, f12.3)') 'PDAFuser:',
     &        '--- Use global rms for observations of ', rms_obs
         WRITE (*, '(a, 8x, a, 1x)') 'PDAFuser:',
     &        '--- Domain localization'
         WRITE (*, '(a, 12x, a, 1x, f12.2)') 'PDAFuser:',
     &        '--- Local influence radius', local_range
         
         IF (locweight > 0) THEN
            WRITE (*, '(a, 12x, a, a)') 'PDAFuser:',
     &           '--- Use distance-dependent weight ',
     &           'for observation errors'

            IF (locweight == 3) THEN
               write (*, '(a, 12x, a)') 'PDAFuser:',
     &              '--- Use regulated weight with mean error variance'
            ELSE IF (locweight == 4) THEN
               write (*, '(a, 12x, a, a)') 'PDAFuser:',
     &              '--- Use regulated weight with ',
     &              'single-point error variance'
            END IF
         END IF
      ENDIF


C     ********************************
C     *** Initialize weight array. ***
C     ********************************

C     Allocate weight array
      ALLOCATE(weight(dim_obs_l))

      if (locweight == 0) THEN
         ! Uniform (unit) weighting
         wtype = 0
         rtype = 0
      else if (locweight == 1) THEN
         ! Exponential weighting
         wtype = 1
         rtype = 0
      ELSE IF (locweight == 2 .OR. locweight
     &        == 3 .OR. locweight == 4) THEN
         ! 5th-order polynomial (Gaspari&Cohn, 1999)
         wtype = 2

         IF (locweight < 3) THEN
            ! No regulated weight
            rtype = 0
         ELSE
            ! Use regulated weight
            rtype = 1
         END IF

      END IF

      IF (locweight == 4) THEN
         ! Allocate array for single observation point
         ALLOCATE(A_obs(1, rank))
      END IF

      DO i=1, dim_obs_l

         ! Control verbosity of PDAF_local_weight
         IF (verbose==1 .AND. i==1) THEN
            verbose_w = 1
         ELSE
            verbose_w = 0
         END IF

         IF (locweight /= 4) THEN
            ! All localizations except regulated weight based on variance at 
            ! single observation point
            CALL PDAF_local_weight(wtype, rtype, 
     &           local_range, srange, distance_l(i),dim_obs_l, rank,
     &           A_l, var_obs, weight(i), verbose_w)
         ELSE
            ! Regulated weight using variance at single observation point
            A_obs(1,:) = A_l(i,:)
            CALL PDAF_local_weight(wtype, rtype,
     &           local_range, srange, distance_l(i),1, rank, A_obs,
     &           var_obs, weight(i), verbose_w)
         END IF
      END DO

      IF (locweight == 4) DEALLOCATE(A_obs)


C     ********************
C     *** Apply weight ***
C     ********************

      DO j = 1, rank
         DO i = 1, dim_obs_l
            C_l(i, j) = ivariance_obs * weight(i) * A_l(i, j)
         END DO
      END DO


C     *** Clean up ***

      DEALLOCATE(weight)
  
      END SUBROUTINE prodRinvA_l_pdaf
