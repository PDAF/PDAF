!$Id: prodrinva_l_pdaf.F90 1155 2011-09-14 08:22:04Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_l_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

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
! often used in EnKF) can be applied to matrix A.
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
! 2005-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state, local_dims
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange, rms_obs
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
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
  INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
  INTEGER :: domain        ! Global domain index
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)     ! Localization weights
  REAL, ALLOCATABLE :: distance(:)   ! Localization distance
  REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL    :: meanvar                 ! Mean variance in observation domain
  REAL    :: svarpovar               ! Mean state plus observation variance
  REAL    :: var_obs                 ! Variance of observation error


! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter == 0) THEN
     verbose = 1
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
  IF (domain > local_range) THEN
     ilow = domain - local_range
  ELSE
     ilow = 1
  ENDIF
  IF (domain + local_range <= dim_state) THEN
     iup = domain + local_range
  ELSE
     iup = dim_state
  ENDIF

! *** Initialize array holding distance of an observation from 
! *** local analysis domain.

  allocate(distance(dim_obs_l))

  init_distance: DO i = ilow, iup
     ! distance between analysis point and current observation
     distance(i - ilow + 1) = ABS( REAL(domain - i))
  END DO init_distance


! *** Initialize weight array

  ! Allocate weight array
  ALLOCATE(weight(dim_obs_l))

  if (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  else if (locweight == 1 .OR. locweight == 3) THEN
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

  end if

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
  
END SUBROUTINE prodRinvA_l_pdaf
