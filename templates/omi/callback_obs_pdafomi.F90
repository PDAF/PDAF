!$Id$
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specifc routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! obs_TYPE_pdafomi (based on the template obs_TYEPE_pdafomi_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this module.
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
  USE obs_TYPE_pdafomi, &
       ONLY: assim_TYPE, init_dim_obs_f_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_TYPE ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_TYPE = 0

  ! Call observation specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_TYPE) CALL init_dim_obs_f_TYPE(step, dim_obs_f_TYPE)

  dim_obs_f = dim_obs_f_TYPE ! + dim_obs_f_TYPE2 ...

END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_f
!!
!! This routine calls the observation-specific
!! routines obs_op_f_TYPE.
!!
SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, ostate_f)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: obs_op_f_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs_f            !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate_f(dim_obs_f)  !< PE-local full observed state

! *** local variables
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! Initialize offset
  offset_obs_f = 0

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_f_TYPE(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

END SUBROUTINE obs_op_f_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for deallocate_obs
!!
!! This routine calls the observation-specific
!! routines deallocate_obs_TYPE.
!!
SUBROUTINE deallocate_obs_pdafomi(step)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: deallocate_obs_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL deallocate_obs_TYPE()

END SUBROUTINE deallocate_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_f
!!
!! This routine calls the observation-specific
!! routines init_obs_f_TYPE.
!!
SUBROUTINE init_obs_f_pdafomi(step, dim_obs_f, observation_f)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_obs_f_TYPE

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

  offset_obs_f = 0

  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL init_obs_f_TYPE(dim_obs_f, observation_f, offset_obs_f)

END SUBROUTINE init_obs_f_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar
!!
!! This routine calls the observation-specific
!! routines init_obsvar_TYPE.
!!
SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_obsvar_TYPE

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

  ! Initialize observation counter
  cnt_obs_f = 0

  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL init_obsvar_TYPE(meanvar, cnt_obs_f)

END SUBROUTINE init_obsvar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_l_TYPE.
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, &
       ONLY: init_dim_obs_l_TYPE

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
  INTEGER :: dim_obs_l_TYPE ! Dimension of observation type A
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets with zero
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL init_dim_obs_l_TYPE(coords_l, local_range, dim_obs_l_TYPE, &
       offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_TYPE

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_l
!!
!! This routine calls the observation-specific
!! routines init_obs_l_TYPE.
!!
SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_obs_l_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  CALL init_obs_l_TYPE(dim_obs_l, observation_l)

END SUBROUTINE init_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for g2l_obs
!!
!! This routine calls the observation-specific
!! routines g2l_obs_TYPE.
!!
SUBROUTINE g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: g2l_obs_TYPE

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

  CALL g2l_obs_TYPE(dim_obs_l, dim_obs_f, ostate_f, ostate_l)

END SUBROUTINE g2l_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the observation-specific
!! routines prodRinvA_l_TYPE.
!!
SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: prodRinvA_l_TYPE

  ! Include variables for localization
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange
  USE mod_parallel, &
       ONLY: mype_filter

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


! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
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

  CALL prodRinvA_l_TYPE(verbose, dim_obs_l, rank, locweight, local_range, &
       srange, A_l, C_l)
  
END SUBROUTINE prodRinvA_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar_l
!!
!! This routine calls the observation-specific
!! routines init_obsvar_l_TYPE.
!!
SUBROUTINE init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_obsvar_l_TYPE

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

  ! Initialize observation counter
  cnt_obs_l = 0

  CALL init_obsvar_l_TYPE(meanvar_l, cnt_obs_l)

END SUBROUTINE init_obsvar_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA
!!
!! This routine calls the observation-specific
!! routines prodRinvA_TYPE.
!!
SUBROUTINE prodRinvA_pdafomi(step, dim_obs_p, ncol, obs_p, A_p, C_p)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: prodRinvA_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  INTEGER, INTENT(in) :: ncol              !< Number of columns in A_p and C_p
  REAL, INTENT(in)    :: obs_p(dim_obs_p)  !< PE-local vector of observations
  REAL, INTENT(in) :: A_p(dim_obs_p, ncol) !< Input matrix
  REAL, INTENT(out)   :: C_p(dim_obs_p, ncol) !< Output matrix


! *************************************
! *** Compute                       ***
! ***                -1             ***
! ***           C = R   A           ***
! *************************************

  CALL prodRinvA_TYPE(ncol, A_p, C_p)
  
END SUBROUTINE prodRinvA_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for add_obs_error
!!
!! This routine calls the observation-specific
!! routines add_obs_err_TYPE.
!!
SUBROUTINE add_obs_error_pdafomi(step, dim_obs_p, C_p)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: add_obs_error_TYPE

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

  CALL add_obs_error_TYPE(dim_obs_p, C_p)
  
END SUBROUTINE add_obs_error_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obscovar
!!
!! This routine calls the observation-specific
!! routines init_obscovar_TYPE.
!!
SUBROUTINE init_obscovar_pdafomi(step, dim_obs, dim_obs_p, covar, m_state_p, &
     isdiag)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_obscovar_TYPE

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

  CALL init_obscovar_TYPE(dim_obs_p, covar, isdiag)
  
END SUBROUTINE init_obscovar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for likelihood
!!
!! This routine calls the observation-specific
!! routines likelihood_TYPE.
!!
SUBROUTINE likelihood_pdafomi(step, dim_obs, obs, resid, lhood)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: likelihood_TYPE

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
  CALL likelihood_TYPE(dim_obs, obs, resid, lhood)

END SUBROUTINE likelihood_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for likelihood_l
!!
!! This routine calls the observation-specific
!! routines likelihood_l_TYPE.
!!
SUBROUTINE likelihood_l_pdafomi(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: likelihood_l_TYPE

  ! Include variables for localization
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange
  USE mod_parallel, &
       ONLY: mype_filter

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


! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain_p < domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
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
  CALL likelihood_l_TYPE(verbose, dim_obs_l, obs_l, resid_l, &
       locweight, local_range, srange, lhood_l)

END SUBROUTINE likelihood_l_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the observation-specific
!! routines localize_covar_TYPE.
!!
SUBROUTINE localize_covar_pdafomi(dim, dim_obs, HP, HPH)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: localize_covar_TYPE

  ! Include variables for localization
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange
  USE mod_parallel, &
       ONLY: mype_filter

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

  IF (mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF

  ! Initialize coordinate array (an array of real)

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
  CALL localize_covar_TYPE(verbose, dim, dim_obs, &
       locweight, local_range, srange, coords, HP, HPH, offset_obs)

END SUBROUTINE localize_covar_pdafomi
