!$Id: likelihood_local.F90 1673 2016-12-02 10:46:31Z lnerger $
!BOP
!
! !ROUTINE: likelihood_local --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood_local(domain, step, dim_obs_l, obs_l, resid_l, likely_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LNETF):
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the likelihood of the ensemble according
! to the difference from the observations (residual)
! and the error distribution of the observations.
! 
! In general this routine is similar to the routine
! prodRinvA_local used for ensemble square root Kalman
! filters. As an addition to this routine, we here have
! to evaluate the likelihood weight according the
! assumed observation error statistics.
!
! This variant is for the Lorenz96 model without
! parallelization. We assume uncorrelated observation
! errors. Thus, for Gaussian errors the observation 
! error covariance matrix is diagonal. As an alternative 
! the use of double-exponetial (Laplace) errors is
! implemented.
! In addition, a localizing weighting of matrix A or the 
! inverse of R by expotential decrease or a 5-th order 
! polynomial of compact support can be applied. This is 
! defined by the variables 'locweight', 'local_range, 
! 'local_range2' and 'srange' in the main program.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state
  USE mod_assimilation, &
       ONLY: local_range, local_range2, locweight, srange, rms_obs, &
       use_obs_mask, obsindx, obsindx_l, obs_err_type
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain             ! Current local analysis domain
  INTEGER, INTENT(in) :: step               ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l          ! Dimension of local observation vector
  REAL, INTENT(in)    :: obs_l(dim_obs_l)   ! Local vector of observations
  REAL, INTENT(inout) :: resid_l(dim_obs_l) ! Input matrix - residual
  REAL, INTENT(out)   :: likely_l           ! Output matrix - log likelihood

! !CALLING SEQUENCE:
! Called by: PDAF_lnetf_analysis    (as U_likelihood_l)
!EOP


! *** local variables ***
  INTEGER :: i             ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER, SAVE :: domain_save = -1   ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)      ! Localization weights
  REAL, ALLOCATABLE :: distance(:)    ! Localization distance
  REAL, ALLOCATABLE :: resid_obs(:)   ! Array for a single row of resid_l
  REAL, ALLOCATABLE :: Rinvresid_l(:) ! R^-1 times residual
  REAL    :: var_obs                  ! Variance of observation error


! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain < domain_save .OR. domain_save < 0) .AND. mype_filter == 0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF
  domain_save = domain

  ! Screen output
  IF (verbose == 1) THEN
     IF (obs_err_type==1) THEN
        WRITE (*, '(8x, a)') &
             '--- Assume double-exponential observation errors'
     ELSE
        WRITE (*, '(8x, a)') &
             '--- Assume Gaussian observation errors'
     END IF
     WRITE (*, '(8x, a, f12.3)') &
          '--- Use global rms for observations of ', rms_obs
     WRITE (*, '(8x, a, 1x)') &
          '--- Domain localization'
     WRITE (*, '(12x, a, 1x, f12.2)') &
          '--- Local influence radius', local_range
     IF (local_range /= local_range2) THEN
        WRITE (*, '(12x, a, f10.4)') &
             '--- Local influence radius on right hand side',local_range2
     END IF

     IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 5) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observed ensemble'
     ELSE IF (locweight == 3 .OR. locweight == 4 .OR. locweight == 6 &
          .OR. locweight == 7) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observation errors'

        IF (locweight == 6) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 7) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF
  
  ! *** initialize numbers
  ivariance_obs = 1.0 / rms_obs**2
  var_obs = rms_obs**2


! ***********************************************
! *** Before computing the likelihood, apply  ***
! *** the localization weight and scale by    ***
! *** observation error variance              ***
! ***                           -1            ***
! ***       Rinvresid = Weight R  resid       ***
! ***                                         ***
! *** Apply a weight matrix with correlations ***
! *** of compact support to residual or the   ***
! *** observation error variance.             ***
! ***********************************************

! *** Initialize array holding distance of an observation from 
! *** local analysis domain.

  ALLOCATE(distance(dim_obs_l))

  init_distance: DO i = 1, dim_obs_l
     ! distance between analysis point and current observation
     IF (.NOT. use_obs_mask) THEN
        distance(i) = ABS( REAL(local_range + 1 - i))
     ELSE
        distance(i) = ABS( REAL(domain - obsindx(obsindx_l(i))))
        if (distance(i) > local_range) distance(i) = ABS(dim_state - distance(i))
     END IF
  END DO init_distance

! *** Initialize weight array

  ! Allocate weight array
  ALLOCATE(weight(dim_obs_l))

  IF (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
  ELSE IF (locweight == 1 .OR. locweight == 3) THEN
     ! Exponential weighting
     wtype = 1
  ELSE IF (locweight == 2 .OR. locweight == 4 .OR. locweight == 5 &
       .OR. locweight == 6 .OR. locweight == 7) THEN
     ! 5th-order polynomial (Gaspari&Cohn, 1999)
     wtype = 2

     IF (locweight < 6) THEN
        ! No regulated weight
        rtype = 0
     ELSE IF (locweight == 6 .OR. locweight == 7) THEN
        ! Use regulated weight
        rtype = 1
     END IF

  END IF

  IF (locweight == 7) THEN
     ! Allocate array for single observation point
     ALLOCATE(resid_obs(1))
  END IF

  DO i=1, dim_obs_l

     ! Control verbosity of PDAF_local_weight
     IF (verbose==1 .AND. i==1) THEN
        verbose_w = 1
     ELSE
        verbose_w = 0
     END IF

     IF (locweight /= 7) THEN
        ! All localizations except regulated weight based on variance at 
        ! single observation point
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             dim_obs_l, 1, resid_l, var_obs, weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        resid_obs(1) = resid_l(i)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             1, 1, resid_obs, var_obs, weight(i), verbose_w)
     END IF
  END DO

  IF (locweight == 7) DEALLOCATE(resid_obs)


! *** Handling of special weighting types ***

  lw2: IF (locweight ==2) THEN
     ! Use square-root of 5th-order polynomial on A

     IF (verbose == 1) THEN
        WRITE (*, '(12x, a)') &
             '--- Use square-root of weight'
     END IF

     DO i = 1, dim_obs_l
        ! Check if weight >0 (Could be <0 due to numerical precision)
        IF (weight(i) > 0.0) THEN
           weight(i) = SQRT(weight(i))
        ELSE
           weight(i) = 0.0
        END IF
     END DO
  END IF lw2


! *** Apply weight

  ALLOCATE(Rinvresid_l(dim_obs_l))

  doweighting: IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 5) THEN

     ! *** Apply weight to matrix A
     DO i = 1, dim_obs_l
        resid_l(i) = weight(i) * resid_l(i)
     END DO

     ! ***       -1
     ! ***  C = R   A 
     DO i = 1, dim_obs_l
        Rinvresid_l(i) = ivariance_obs * resid_l(i)
     END DO
  
  ELSE doweighting

     ! *** Apply weight to matrix R only
     DO i = 1, dim_obs_l
        Rinvresid_l(i) = ivariance_obs * weight(i) * resid_l(i)
     END DO
     
  END IF doweighting


! ******************************
! *** Compute log likelihood ***
! ******************************

  IF (obs_err_type==0) THEN

     ! Gaussian errors
     ! Calculate exp(-0.5*resid^T*R^-1*resid)
     CALL dgemv('t', dim_obs_l, 1, 0.5, resid_l, &
          dim_obs_l, Rinvresid_l, 1, 0.0, likely_l, 1)
     likely_l = EXP(-likely_l)

  ELSE

     ! Double-exponential errors
     ! Calculate exp(-SUM(ABS(resid)))
     likely_l = 0.0
     DO i = 1, dim_obs_l
        likely_l = likely_l + ABS(Rinvresid_l(i))
     END DO
     likely_l = EXP(-likely_l)

  END IF

! *** Clean up ***

  DEALLOCATE(weight, distance, Rinvresid_l)

END SUBROUTINE likelihood_local
