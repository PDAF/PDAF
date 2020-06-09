!$Id: callback_obs_pdafomi.F90 488 2020-06-07 11:00:55Z lnerger $
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
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)

  ! Include functions for different observations
  USE obs_gp_pdafomi, ONLY: init_dim_obs_l_gp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets
  INTEGER :: dim_obs_l_gp ! Dimension of observation type A


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL init_dim_obs_l_gp(domain_p, step, dim_obs_f, dim_obs_l_gp, offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_gp

END SUBROUTINE init_dim_obs_l_pdafomi



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



!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type
!!
SUBROUTINE localize_covar_pdafomi(dim, dim_obs, HP, HPH)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_localize_covar

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
  DO i=1, n_obstypes
     CALL PDAFomi_localize_covar(obs_f_all(i)%ptr, dim, locweight, local_range, srange, &
          coords, HP, HPH, offset_obs, verbose)
  END DO

END SUBROUTINE localize_covar_pdafomi
