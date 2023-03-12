!$Id$
!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  USE obs_OBSTYPE_pdafomi, ONLY: assim_OBSTYPE, init_dim_obs_OBSTYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_OBSTYPE ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_pdafomi: complete interface to observation modules'

  ! Initialize number of observations
  dim_obs_OBSTYPE = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_OBSTYPE) CALL init_dim_obs_OBSTYPE(step, dim_obs_OBSTYPE)

  dim_obs = dim_obs_OBSTYPE ! + dim_obs_OBSTYPE2 ...

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_OBSTYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_OBSTYPE_pdafomi, ONLY: obs_op_OBSTYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/obs_op_pdafomi: complete interface to observation modules'

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi
  CALL obs_op_OBSTYPE(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  USE obs_OBSTYPE_pdafomi, ONLY: init_dim_obs_l_OBSTYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/init_dim_obs_l_pdafomi: complete interface to observation modules'

  ! Call init_dim_obs_l specific for each observation
  CALL init_dim_obs_l_OBSTYPE(domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for localize_covar
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type to apply covariance
!! localization in the LEnKF.
!!
SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)

  ! Include functions for different observations
  USE obs_OBSTYPE_pdafomi, ONLY: localize_covar_OBSTYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs               !< number of observations
  REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

! *** local variables ***
  REAL, ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries


! **********************
! *** INITIALIZATION ***
! **********************

  WRITE (*, *) 'TEMPLATE callback_obs_pdafomi.F90/localize_covar_pdafomi: complete interface to observation modules'

  ! Initialize coordinate array

  ! One needs to provide the array COORDS_P holding the coordinates of each
  ! element of the process-local state vector. Each column of the array holds
  ! the information for one element. The array can be initialized here using
  ! information on the model grid.

  ! ALLOCATE(coords_p(NROWS, dim_p))

  ! coords_p = ...


! *************************************
! *** Apply covariance localization ***
! *************************************

  ! Call localize_covar specific for each observation
  CALL localize_covar_OBSTYPE(dim_p, dim_obs, HP_p, HPH, coords_p)


! ****************
! *** Clean up ***
! ****************

  ! DEALLOCATE(coords_p)

END SUBROUTINE localize_covar_pdafomi
