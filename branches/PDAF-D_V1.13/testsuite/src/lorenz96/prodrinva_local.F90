!$Id: prodrinva_local.F90 1633 2016-08-15 08:02:27Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_local --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_local(domain, step, dim_obs_l, rank, obs_l, A_l, C_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the product of the inverse of the local
! observation error covariance matrix with
! the matrix of locally observed ensemble 
! perturbations.
! Next to computing the product,  a localizing 
! weighting (similar to covariance localization 
! often used in EnKF) can be applied to matrix A
! or the inverse observation error covariance matrix.
!
! This variant is for the Lorenz96 model without
! parallelization. We assume a diagonal observation
! error covariance matrix with constant variances. 
! Thus, the product can be implemented efficiently 
! as a scaling of each element of the input matrix
! by the inverse variance. In addition, a 
! localizing weighting of matrix A or the inverse of R
! by expotential decrease or a 5-th order polynomial 
! of compact support can be applied. This is defined 
! by the variables 'locweight', 'local_range, 
! 'local_range2' and 'srange' in the main program.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state
  USE mod_assimilation, &
       ONLY: local_range, local_range2, locweight, srange, rms_obs, &
       use_obs_mask, obsindx, obsindx_l
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain            ! Current local analysis domain
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  ! Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis    (as U_prodRinvA_l)
! Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)
!EOP


! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)     ! Localization weights
  REAL, ALLOCATABLE :: distance(:)   ! Localization distance
  REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL    :: var_obs                 ! Variance of observation error


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
! ***                                 -1      ***
! ***     A = Weight*A or C = Weight R  A     ***
! ***                                         ***
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
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
     ALLOCATE(A_obs(1, rank))
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
             dim_obs_l, rank, A_l, var_obs, weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        A_obs(1,:) = A_l(i,:)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             1, rank, A_obs, var_obs, weight(i), verbose_w)
     END IF
  END DO

  IF (locweight == 7) DEALLOCATE(A_obs)


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

  doweighting: IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 5) THEN

     ! *** Apply weight to matrix A
     DO j = 1, rank
        DO i = 1, dim_obs_l
           A_l(i, j) = weight(i) * A_l(i, j)
        END DO
     END DO

     ! ***       -1
     ! ***  C = R   A 
     DO j = 1, rank
        DO i = 1, dim_obs_l
           C_l(i, j) = ivariance_obs * A_l(i, j)
        END DO
     END DO
  
  ELSE doweighting

     ! *** Apply weight to matrix R only
     DO j = 1, rank
        DO i = 1, dim_obs_l
           C_l(i, j) = ivariance_obs * weight(i) * A_l(i, j)
        END DO
     END DO
     
  END IF doweighting


! *** Clean up ***

  DEALLOCATE(weight, distance)

END SUBROUTINE prodRinvA_local
