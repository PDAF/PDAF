!$Id$
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
! weighting (similar to covariance localization 
! often used in EnKF) can be applied to matrix A.
!
! Implementation for the 2D online tutorial example.
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: local_range, locweight, srange
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, rms_obs_A, prodRinvA_l_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, rms_obs_B, prodRinvA_l_B
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_thread_num
#endif

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
  INTEGER :: verbose       ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  INTEGER, SAVE :: mythread          ! Thread variable for OpenMP

!$OMP THREADPRIVATE(mythread, domain_save)


! *** NO CHANGES REQUIRED BELOW IF THE OBSERVATION ERROR COVARIANCE MATRIX IS DIAGONAL ***


! **********************
! *** INITIALIZATION ***
! **********************

! For OpenMP parallelization, determine the thread index
#if defined (_OPENMP)
  mythread = omp_get_thread_num()
#else
  mythread = 0
#endif

  IF (domain_p <= domain_save .OR. domain_save < 0) THEN
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
          '--- Use global rms for observations A of ', rms_obs_A
     WRITE (*, '(8x, a, f12.3)') &
          '--- Use global rms for observations B of ', rms_obs_B
  ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

  IF (assim_A) CALL prodRinvA_l_A(verbose, dim_obs_l, rank, locweight, local_range, &
       srange, A_l, C_l)
  IF (assim_B) CALL prodRinvA_l_B(verbose, dim_obs_l, rank, locweight, local_range, &
       srange, A_l, C_l)
  
END SUBROUTINE prodRinvA_l_pdaf
