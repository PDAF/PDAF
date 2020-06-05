!$Id: mod_obs_A_pdaf.F90 251 2019-11-19 08:43:39Z lnerger $
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
  USE obs_A_pdafomi, ONLY: assim_A, init_dim_obs_f_A
  USE obs_B_pdafomi, ONLY: assim_B, init_dim_obs_f_B
  USE obs_C_pdafomi, ONLY: assim_C, init_dim_obs_f_C

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_A ! Observation dimensions
  INTEGER :: dim_obs_f_B ! Observation dimensions
  INTEGER :: dim_obs_f_C ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_A = 0
  dim_obs_f_B = 0
  dim_obs_f_C = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
  IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)
  IF (assim_C) CALL init_dim_obs_f_C(step, dim_obs_f_C)

  dim_obs_f = dim_obs_f_A + dim_obs_f_B + dim_obs_f_C

END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_f
!!
!! This routine calls the observation-specific
!! routines obs_op_f_TYPE.
!!
SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, ostate_f)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: obs_op_f_A, obs_A => thisobs
  USE obs_B_pdafomi, ONLY: obs_op_f_B, obs_B => thisobs
  USE obs_C_pdafomi, ONLY: obs_op_f_C, obs_C => thisobs

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

  ! Initialize overall observation pointer array
  if (.not.allocated(obs_f_all)) then
     n_obstypes = 3
     allocate(obs_f_all(n_obstypes))
     obs_f_all(1)%ptr => obs_A
     obs_f_all(2)%ptr => obs_B
     obs_f_all(3)%ptr => obs_C
  endif

  ! Initialize offset
  offset_obs_f = 0

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_f_A(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_B(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_C(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

END SUBROUTINE obs_op_f_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs
  USE obs_B_pdafomi, ONLY: obs_B => thisobs
  USE obs_C_pdafomi, ONLY: obs_C => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_A)
  CALL PDAFomi_deallocate_obs(obs_B)
  CALL PDAFomi_deallocate_obs(obs_C)

END SUBROUTINE deallocate_obs_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l
  USE obs_C_pdafomi, ONLY: obs_C => thisobs, obs_C_l => thisobs_l

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
  INTEGER :: dim_obs_l_A ! Dimension of observation type A
  INTEGER :: dim_obs_l_B ! Dimension of observation type B
  INTEGER :: dim_obs_l_C ! Dimension of observation type C
  INTEGER :: offset_obs_l, offset_obs_f     ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL PDAFomi_init_dim_obs_l(obs_A_l, obs_A, coords_l, local_range, dim_obs_l_A, &
       offset_obs_l, offset_obs_f)
  CALL PDAFomi_init_dim_obs_l(obs_B_l, obs_B, coords_l, local_range, dim_obs_l_B, &
       offset_obs_l, offset_obs_f)
  CALL PDAFomi_init_dim_obs_l(obs_C_l, obs_C, coords_l, local_range, dim_obs_l_C, &
       offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_A + dim_obs_l_B + dim_obs_l_C

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the routine PDAFomi_prodRinvA_l
!! for each observation type
!!
SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA_l

  ! Include variables for localization
  USE mod_assimilation, ONLY: local_range, locweight, srange
  ! Include filter process rank
  USE mod_parallel_pdaf, ONLY: mype_filter

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
  INTEGER :: i                       ! Loop counter
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index


! **********************
! *** INITIALIZATION ***
! **********************

  ! Set verbosity flag (Screen output for first analysis domain)
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

  DO i=1, n_obstypes
     CALL PDAFomi_prodRinvA_l(obs_l_all(i)%ptr, obs_f_all(i)%ptr, dim_obs_l, rank, &
          locweight, local_range, srange, A_l, C_l, verbose)
  END DO
  
END SUBROUTINE prodRinvA_l_pdafomi
