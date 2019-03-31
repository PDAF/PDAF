!$Id: likelihood_l_pdaf.F90 1678 2016-12-11 12:32:53Z lnerger $
!BOP
!
! !ROUTINE: likelihood_l_pdaf --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l)

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
! This routine is called by all filter processes.
!
! Implementation for the dummy model with domain
! decomposition. Here, we assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance. In addition, a 
! localizing weighting of matrix A or the inverse of R
! by expotential decrease or a 5-th order polynomial 
! of compact support can be applied. This is defined 
! by the variables 'locweight', 'local_range, and 
! 'srange' in the main program.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code based on prodRinvA_l_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state, local_dims
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange, rms_obs
  USE mod_parallel, &
       ONLY: mype_filter
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p           ! Current local analysis domain
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
  INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
  INTEGER :: domain        ! Global domain index
  INTEGER, SAVE :: domain_save = -1   ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)      ! Localization weights
  REAL, ALLOCATABLE :: distance(:)    ! Localization distance
  REAL, ALLOCATABLE :: resid_obs(:)   ! Array for a single row of resid_l
  REAL, ALLOCATABLE :: Rinvresid_l(:) ! R^-1 times residual
  REAL    :: var_obs                  ! Variance of observation error
  INTEGER, SAVE :: mythread           ! Thread variable for OpenMP

! For OpenMP set the domain and the thread index to be thread private
!$OMP THREADPRIVATE(mythread, domain_save)


! **********************
! *** INITIALIZATION ***
! **********************

! For OpenMP parallelization, determine the thread index
#if defined (_OPENMP)
  mythread = omp_get_thread_num()
#else
  mythread = 0
#endif

  IF ((domain_p < domain_save .OR. domain_save < 0) .AND. mype_filter == 0) THEN
      verbose = 1

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) verbose = 0
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p

  ! Screen output
  IF (verbose == 1) THEN
     WRITE (*, '(8x, a, f12.3)') &
          '--- Use global rms for observations of ', rms_obs
     WRITE (*, '(8x, a, 1x)') &
          '--- Domain localization'
     WRITE (*, '(12x, a, 1x, f12.2)') &
          '--- Local influence radius', local_range

     IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 5) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observed ensemble'
     ELSE IF (locweight == 3 .OR. locweight == 4 .OR. locweight == 6 &
          .OR. locweight == 7) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observation errors'

        IF (locweight == 6) THEN
           WRITE (*, '(12x, a)') &
                '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 7) THEN
           WRITE (*, '(12x, a)') &
                '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF
  
  ! *** initialize numbers
  ivariance_obs = 1.0 / rms_obs**2
  var_obs = rms_obs**2


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

  ! Get domain index in global grid
  domain = domain_p
  DO i = 1, mype_filter
     domain = domain + local_dims(i)
  ENDDO

  ! Get grid index range for local observations
  IF (domain > INT(local_range)) THEN
     ilow = domain - INT(local_range)
  ELSE
     ilow = 1
  ENDIF
  IF (domain + INT(local_range) <= dim_state) THEN
     iup = domain + INT(local_range)
  ELSE
     iup = dim_state
  ENDIF

! *** Initialize array holding distance of an observation from 
! *** local analysis domain.

  ALLOCATE(distance(dim_obs_l))

  init_distance: DO i = ilow, iup
     ! distance between analysis point and current observation
     distance(i - ilow + 1) = ABS( REAL(domain - i))
  END DO init_distance


! *** Initialize weight array

  ! Allocate weight array
  ALLOCATE(weight(dim_obs_l))

  IF (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  ELSE IF (locweight == 1 .OR. locweight == 3) THEN
     ! Exponential weighting
     wtype = 1
     rtype = 0
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

  ! Gaussian errors: Calculate exp(-0.5*resid^T*R^-1*resid)
  CALL dgemv('t', dim_obs_l, 1, 0.5, resid_l, &
       dim_obs_l, Rinvresid_l, 1, 0.0, likely_l, 1)
  likely_l = EXP(-likely_l)


! *** Clean up ***

  DEALLOCATE(weight, distance, Rinvresid_l)
  
END SUBROUTINE likelihood_l_pdaf
