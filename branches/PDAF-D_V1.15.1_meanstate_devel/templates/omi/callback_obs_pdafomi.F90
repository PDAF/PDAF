!$Id$
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
  USE obs_TYPE_pdafomi, ONLY: assim_TYPE, init_dim_obs_f_TYPE

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

  ! Call observation-specific routines
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


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! Initialize offset
  offset_obs_f = 0

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_f_TYPE(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

END SUBROUTINE obs_op_f_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)

  ! Include functions for different observations
  USE obs_TYPE_pdafomi, ONLY: init_dim_obs_l_TYPE

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets
  INTEGER :: dim_obs_l_TYPE ! Dimension of observation type


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL init_dim_obs_l_TYPE(domain_p, step, dim_obs_f, dim_obs_l_TYPE, offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_TYPE ! + dim_obs_l_TYPE2 ...

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
  USE obs_TYPE_pdafomi, ONLY: obs_TYPE => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_TYPE)

END SUBROUTINE deallocate_obs_pdafomi
