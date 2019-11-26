!$Id: prodrinva_l_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
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
  USE mod_assim_pdaf, &
       ONLY: local_range, locweight, srange, rms_obs, local_obs_nod2d, &
       ocoord_n2d, bias_obs, loc_radius, ivariance_obs_l, distance
  USE mod_parallel_pdaf, &
       ONLY: mype_filter
  USE o_mesh, &
       ONLY: coord_nod2D
  USE o_param, &
       ONLY: r_earth
  USE g_rotate_grid

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
! Called by: PDAF_lestkf_analysis   (as U_prodRinvA_l)
! Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)
!EOP


! *** local variables ***
  INTEGER :: lnode         ! Counter for local node
  INTEGER :: i, j          ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
!  REAL :: ivariance_obs    ! Inverse of variance of the observations
  REAL, ALLOCATABLE :: var_obs_l(:)          ! Variance of observation error
  REAL :: wc_coord(2)      ! Coordinates of current water column
  REAL :: o_coord(2)       ! Coordinates of observation
  REAL :: dist2d(2)        ! Distance vector between water column and observation
  REAL, ALLOCATABLE :: weight(:)     ! Localization weights
  REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL, PARAMETER :: pi=3.14159265358979

! **********************
! *** INITIALIZATION ***
! **********************
  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p

  ! Screen output
  IF (verbose == 1) THEN
     WRITE (*, '(a, 5x, a, f12.3)') &
          'FESOM-PDAF', '--- Use global rms for observations of ', rms_obs
     IF (bias_obs /= 0.0) THEN
        WRITE (*, '(a, 5x, a, f12.3)') &
             'FESOM-PDAF', '--- Use global observation bias of ', bias_obs
     END IF
     WRITE (*, '(a, 5x, a, 1x)') &
          'FESOM-PDAF', '--- Domain localization'
     WRITE (*, '(a, 8x, a, 1x, es11.3)') &
          'FESOM-PDAF', '--- Local influence radius', local_range

     IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 5) THEN
        WRITE (*, '(a, 8x, a)') &
             'FESOM-PDAF', '--- Use distance-dependent weight for observed ensemble'
     ELSE IF (locweight == 3 .OR. locweight == 4 .OR. locweight == 6 &
          .OR. locweight == 7) THEN
        WRITE (*, '(a, 8x, a)') &
             'FESOM-PDAF', '--- Use distance-dependent weight for observation errors'

        IF (locweight == 6) THEN
           write (*, '(a, 8x, a)') &
                'FESOM-PDAF', '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 7) THEN
           write (*, '(a, 8x, a)') &
                'FESOM-PDAF', '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

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

  DO i = 1, dim_obs_l

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
             dim_obs_l, rank, A_l, var_obs_l(i), weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        A_obs(1,:) = A_l(i,:)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             1, rank, A_obs, var_obs_l(i), weight(i), verbose_w)
     END IF
  END DO
  
 IF (locweight == 7) DEALLOCATE(A_obs)

! *** Handling of special weighting types ***

  lw2: IF (locweight ==2) THEN
     ! Use square-root of 5th-order polynomial on A

     IF (verbose == 1) THEN
        WRITE (*, '(a, 12x, a)') &
             'FESOM-PDAF', '--- Use square-root of weight'
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
           C_l(i, j) = ivariance_obs_l(i) * A_l(i, j)
        END DO
     END DO
  
  ELSE doweighting

     ! *** Apply weight to matrix R only
     DO j = 1, rank
        DO i = 1, dim_obs_l
           C_l(i, j) = ivariance_obs_l(i) * weight(i) * A_l(i, j)
        END DO
     END DO
     
  END IF doweighting

! *** Clean up ***

  DEALLOCATE(weight)

END SUBROUTINE prodRinvA_l_pdaf
