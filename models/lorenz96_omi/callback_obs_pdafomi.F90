!$Id: callback_obs_pdafomi.F90 451 2020-05-24 08:26:18Z lnerger $
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specifc routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this file. It is recommended to keep the order of the calls
!! consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs_f
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_f_TYPE.
!!
SUBROUTINE init_dim_obs_f_pdafomi(step, dim_obs_f)

  ! Include functions for different observations
  USE obs_gp_pdafomi, ONLY: assim_gp, init_dim_obs_f_gp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_gp ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_gp = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_gp) CALL init_dim_obs_f_gp(step, dim_obs_f_gp)

  dim_obs_f = dim_obs_f_gp

END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_f
!!
!! This routine calls the observation-specific
!! routines obs_op_f_TYPE.
!!
SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, ostate_f)

  ! Include functions for different observations
  USE obs_gp_pdafomi, ONLY: obs_op_f_gp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs_f            !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate_f(dim_obs_f)  !< PE-local full observed state

! *** local variables
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! Initialize offset
  offset_obs_f = 0

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_f_gp(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

END SUBROUTINE obs_op_f_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_f
!!
!! This routine calls the routine PDAFomi_init_obs_f
!! for each observation type
!!
SUBROUTINE init_obs_f_pdafomi(step, dim_obs_f, observation_f)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_f
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector

! *** local variables ***
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  ! Initialize offset (it will be incremented in PDAFomi_init_obs_f)
  offset_obs_f = 0

  ! The order of the calls has to be consistent with those in obs_op_f_pdafomi
  CALL PDAFomi_init_obs_f(gpobs, dim_obs_f, observation_f, offset_obs_f)

END SUBROUTINE init_obs_f_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar
!!
!! This routine calls the routine PDAFomi_init_obsvar_f
!! for each observation type
!!
SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_f
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
  REAL, INTENT(out)   :: meanvar       !< Mean observation error variance

! *** Local variables ***
  INTEGER :: cnt_obs_f


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_f = 0

  CALL PDAFomi_init_obsvar_f(gpobs, meanvar, cnt_obs_f)

END SUBROUTINE init_obsvar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l
  ! Include observation types
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  ! Include localization radius and local coordinates
  USE mod_assimilation, &   
       ONLY: local_range, coords_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: dim_obs_l_gp ! Dimension of observation type A
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL PDAFomi_init_dim_obs_l(gpobs_l, gpobs, coords_l, local_range, dim_obs_l_gp, &
       offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_gp

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for g2l_obs
!!
!! This routine calls the routine PDAFomi_g2l_obs
!! for each observation type
!!
SUBROUTINE g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_g2l_obs
  ! Include observation types
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local observation vector
  INTEGER, INTENT(in) :: dim_obs_l  !< Dimension of local observation vector
  REAL, INTENT(in)    :: ostate_f(dim_obs_f)   !< Full PE-local obs.ervation vector
  REAL, INTENT(out)   :: ostate_l(dim_obs_l)   !< Observation vector on local domain


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  CALL PDAFomi_g2l_obs(gpobs_l, gpobs, ostate_f, ostate_l)

END SUBROUTINE g2l_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_l
!!
!! This routine calls the routine PDAFomi_init_obs_l
!! for each observation type
!!
SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_l
  ! Include observation types
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  CALL PDAFomi_init_obs_l(gpobs_l, gpobs, observation_l)

END SUBROUTINE init_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar_l
!!
!! This routine calls the routine PDAFomi_init_obsvar_l
!! for each observation type
!!
SUBROUTINE init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_l
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance

! *** Local variables ***
  INTEGER :: cnt_obs_l


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_l = 0

  ! The order of the calls is not relevant
  CALL PDAFomi_init_obsvar_l(gpobs_l, gpobs, meanvar_l, cnt_obs_l)

END SUBROUTINE init_obsvar_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the routine PDAFomi_prodRinvA_l
!! for each observation type
!!
SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA_l
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  ! Include variables for localization
  USE mod_assimilation, ONLY: local_range, locweight, srange
  ! Include filter process rank
  USE mod_parallel, ONLY: mype_filter
#if defined (_OPENMP)
  ! Include OpenMP function to determine verbosity for OpenMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

! *** local variables ***
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  INTEGER, SAVE :: mythread          ! Thread variable for OpenMP

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

  ! Set verbosity flag (Screen output for first analysis domain)
  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) verbose = 0
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p


! *************************************
! *** Compute                       ***
! ***                  -1           ***
! ***           C = W R   A         ***
! ***                               ***
! *** where W are the localization  ***
! *** weights.                      ***
! *************************************

  CALL PDAFomi_prodRinvA_l(gpobs_l, gpobs, dim_obs_l, rank, &
       locweight, local_range, srange, A_l, C_l, verbose)
  
END SUBROUTINE prodRinvA_l_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for likelihood_l
!!
!! This routine calls the routine PDAFomi_likelihood_l
!! for each observation type
!!
SUBROUTINE likelihood_l_pdafomi(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_likelihood_l
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs, gpobs_l => thisobs_l

  ! Include variables for localization
  USE mod_assimilation, ONLY: local_range, locweight, srange
  ! Include filter process rank
  USE mod_parallel, ONLY: mype_filter
#if defined (_OPENMP)
  ! Include OpenMP function to determine verbosity for OpenMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p           ! Current local analysis domain
  INTEGER, INTENT(in) :: step               !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l          !< PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs_l(dim_obs_l)   !< PE-local vector of observations
  REAL, INTENT(inout) :: resid_l(dim_obs_l) !< Input vector of residuum
  REAL, INTENT(out)   :: lhood_l            !< Output vector - log likelihood

! *** local variables ***
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  INTEGER, SAVE :: mythread          ! Thread variable for OpenMP

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

  ! Set verbosity flag (Screen output for first analysis domain)
  IF ((domain_p < domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) verbose = 0
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p


! ********************************
! *** Compute local likelihood ***
! ********************************

  ! Initialize likelihood value before starting computation
  lhood_l = 0.0

  ! Increment likelihood
  CALL PDAFomi_likelihood_l(gpobs_l, gpobs, resid_l, locweight, &
       local_range, srange, lhood_l, verbose)

END SUBROUTINE likelihood_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA
!!
!! This routine calls the routine PDAFomi_prodRinvA
!! for each observation type
!!
SUBROUTINE prodRinvA_pdafomi(step, dim_obs_p, ncol, obs_p, A_p, C_p)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  INTEGER, INTENT(in) :: ncol              !< Number of columns in A_p and C_p
  REAL, INTENT(in)    :: obs_p(dim_obs_p)  !< PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p, ncol) !< Input matrix
  REAL, INTENT(out)   :: C_p(dim_obs_p, ncol) !< Output matrix


! *************************************
! *** Compute                       ***
! ***                -1             ***
! ***           C = R   A           ***
! *************************************

  CALL PDAFomi_prodRinvA(gpobs, ncol, A_p, C_p)
  
END SUBROUTINE prodRinvA_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for likelihood
!!
!! This routine calls the routine PDAFomi_likelihood
!! for each observation type
!!
SUBROUTINE likelihood_pdafomi(step, dim_obs, obs, resid, lhood)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_likelihood
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: dim_obs          !< PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs(dim_obs)     !< PE-local vector of observations
  REAL, INTENT(in)    :: resid(dim_obs)   !< Input vector of residuum
  REAL, INTENT(out)   :: lhood            !< Output vector - log likelihood


! **************************
! *** Compute likelihood ***
! **************************

  ! Initialize likelihood value before starting computation
  lhood = 0.0

  ! Increment likelihood
  CALL PDAFomi_likelihood(gpobs, dim_obs, obs, resid, lhood)

END SUBROUTINE likelihood_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for add_obs_error
!!
!! This routine calls the routine PDAFomi_add_obs_error
!! for each observation type
!!
SUBROUTINE add_obs_error_pdafomi(step, dim_obs_p, C_p)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_add_obs_error
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  REAL, INTENT(inout) :: C_p(dim_obs_p,dim_obs_p) ! Matrix to which R is added


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

  CALL PDAFomi_add_obs_error(gpobs, dim_obs_p, C_p)
  
END SUBROUTINE add_obs_error_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obscovar
!!
!! This routine calls the routine PDAFomi_init_obscovar
!! for each observation type
!!
SUBROUTINE init_obscovar_pdafomi(step, dim_obs, dim_obs_p, covar, m_state_p, &
     isdiag)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obscovar
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of observation vector
  INTEGER, INTENT(in) :: dim_obs_p            !< PE-local dimension of obs. vector
  REAL, INTENT(out) :: covar(dim_obs,dim_obs) !< Observation error covar. matrix
  REAL, INTENT(in) :: m_state_p(dim_obs_p)    !< Observation vector
  LOGICAL, INTENT(out) :: isdiag              !< Whether matrix R is diagonal


! *************************************
! ***   Initialize covariances      ***
! *************************************

  CALL PDAFomi_init_obscovar(gpobs, dim_obs_p, covar, isdiag)
  
END SUBROUTINE init_obscovar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type
!!
SUBROUTINE localize_covar_pdafomi(dim, dim_obs, HP, HPH)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_localize_covar
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  ! Include variables for localization
  USE mod_assimilation, ONLY: local_range, locweight, srange
  ! Include filter process rank
  USE mod_parallel, ONLY: mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim                   !< State dimension
  INTEGER, INTENT(in) :: dim_obs               !< number of observations
  REAL, INTENT(inout) :: HP(dim_obs, dim)      !< Matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

! *** local variables ***
  INTEGER :: i                       ! Counter
  INTEGER :: verbose                 ! Verbosity flag
  REAL, ALLOCATABLE :: coords(:,:)   ! Coordinates of sstate vector entries
  INTEGER :: offset_obs              ! Count offset of an observation type in full obs. vector


! **********************
! *** INITIALIZATION ***
! **********************

  ! Set verbosity flag (only first process writes)
  IF (mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF

  ! Initialize coordinate array

  ALLOCATE(coords(1, dim))

  DO i=1, dim
     coords(1, i) = REAL(i)
  END DO


! *************************************
! *** Apply covariance localization ***
! *************************************

  ! Initialize offset
  offset_obs = 0

  ! Apply localization
  CALL PDAFomi_localize_covar(gpobs, dim, locweight, local_range, srange, &
       coords, HP, HPH, offset_obs, verbose)

END SUBROUTINE localize_covar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obserr_f_pdaf
!!
!! This routine calls the routine PDAFomi_init_obserr_f
!! for each observation type
!!
SUBROUTINE init_obserr_f_pdafomi(step, dim_obs_f, obs_f, obserr_f)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obserr_f
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f           !< Full dimension of observation vector
  REAL, INTENT(in)    :: obs_f(dim_obs_f)    !< Full observation vector
  REAL, INTENT(out)   :: obserr_f(dim_obs_f) !< Full observation error stddev


! *****************************************************************************
! *** Initialize vector of observation errors for generating synthetic obs. ***
! *****************************************************************************

  CALL PDAFomi_init_obserr_f(gpobs, obserr_f)
  
END SUBROUTINE init_obserr_f_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for deallocate_obs
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!! for each observation type
!!
SUBROUTINE deallocate_obs_pdafomi(step)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_deallocate_obs
  ! Include observation types (rename generic name)
  USE obs_gp_pdafomi, ONLY: gpobs => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(gpobs)

END SUBROUTINE deallocate_obs_pdafomi
