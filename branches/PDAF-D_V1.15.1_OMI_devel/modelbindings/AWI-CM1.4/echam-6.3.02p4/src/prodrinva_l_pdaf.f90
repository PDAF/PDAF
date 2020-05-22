!$Id: prodrinva_l_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_l_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the product of the inverse of the local
! observation error covariance matrix with
! the matrix of locally observed ensemble 
! perturbations.
! Next to computing the product, a localizing 
! weighting ("observation localization") can be
! applied to matrix A.
!
! This routine is called by all filter processes.
!
! This template is for a constant observation error
! given by the variable rms_obs. In this case, only
! the initialization of the distance has to be
! implemented.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE mod_assim_pdaf, &
       ONLY: local_range, locweight, srange, rms_obs
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
  REAL(dp), INTENT(in)    :: obs_l(dim_obs_l)  ! Local vector of observations
  REAL(dp), INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix
  REAL(dp), INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis    (as U_prodRinvA_l)
! Called by: PDAF_lestkf_analysis   (as U_prodRinvA_l)
! Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)
!EOP


! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
  INTEGER :: domain        ! Global domain index
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL(dp)    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL(dp), ALLOCATABLE :: weight(:)     ! Localization weights
  REAL(dp), ALLOCATABLE :: distance(:)   ! Localization distance
  REAL(dp), ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL(dp)    :: meanvar                 ! Mean variance in observation domain
  REAL(dp)    :: svarpovar               ! Mean state plus observation variance
  REAL(dp)    :: var_obs                 ! Variance of observation error


! **********************
! *** INITIALIZATION ***
! **********************
  if (mype_filter==0) write (*,*) 'ECHAM-PDAF TEMPLATE prodrinva_l_pdaf'

if (1==2) then     !Temporary deactivation
  IF (domain_p <= domain_save .OR. domain_save < 0) THEN
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
  
  ! *** initialize numbers (this is for constant observation errors)
  ivariance_obs = 1.0 / rms_obs**2
  var_obs = rms_obs**2


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

! *** Initialize array holding distance of an observation from 
! *** local analysis domain.

  allocate(distance(dim_obs_l))

!   DO i = ?, ?
     ! distance between analysis point and current observation
!      distance(i - ilow + 1) = ????
!   END DO

! *** NO CHANGES REQUIRED BELOW IF OBSERVATION ERRORS ARE CONSTANT ***

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

endif   !Temporary deactivation
  
END SUBROUTINE prodRinvA_l_pdaf
