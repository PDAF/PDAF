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

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_A ! Observation dimensions
  INTEGER :: dim_obs_f_B ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_A = 0
  dim_obs_f_B = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
  IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)

  dim_obs_f = dim_obs_f_A + dim_obs_f_B

END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_f
!!
!! This routine calls the observation-specific
!! routines obs_op_f_TYPE.
!!
SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, ostate_f)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: obs_op_f_A
  USE obs_B_pdafomi, ONLY: obs_op_f_B

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
  CALL obs_op_f_A(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_B(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

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

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_A)
  CALL PDAFomi_deallocate_obs(obs_B)

END SUBROUTINE deallocate_obs_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs
  USE obs_B_pdafomi, ONLY: obs_B => thisobs

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
  CALL PDAFomi_init_obs_f(obs_A, dim_obs_f, observation_f, offset_obs_f)
  CALL PDAFomi_init_obs_f(obs_B, dim_obs_f, observation_f, offset_obs_f)

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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs
  USE obs_B_pdafomi, ONLY: obs_B => thisobs

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

  CALL PDAFomi_init_obsvar_f(obs_A, meanvar, cnt_obs_f)
  CALL PDAFomi_init_obsvar_f(obs_B, meanvar, cnt_obs_f)

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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

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
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


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

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_A + dim_obs_l_B

END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_l
!!
!! This routine calls the routine PDAFomi_init_obs_l
!! for each observation type
!!
SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_l
  ! Include observation types (rename generic name)
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  CALL PDAFomi_init_obs_l(obs_A_l, obs_A, observation_l)
  CALL PDAFomi_init_obs_l(obs_B_l, obs_B, observation_l)

END SUBROUTINE init_obs_l_pdafomi



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
  ! Include observation types (rename generic name)
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

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

  CALL PDAFomi_g2l_obs(obs_A_l, obs_A, ostate_f, ostate_l)
  CALL PDAFomi_g2l_obs(obs_B_l, obs_B, ostate_f, ostate_l)

END SUBROUTINE g2l_obs_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

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

  CALL PDAFomi_prodRinvA_l(obs_A_l, obs_A, dim_obs_l, rank, &
       locweight, local_range, srange, A_l, C_l, verbose)
  CALL PDAFomi_prodRinvA_l(obs_B_l, obs_B, dim_obs_l, rank, &
       locweight, local_range, srange, A_l, C_l, verbose)
  
END SUBROUTINE prodRinvA_l_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

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

  CALL PDAFomi_init_obsvar_l(obs_A_l, obs_A, meanvar_l, cnt_obs_l)
  CALL PDAFomi_init_obsvar_l(obs_B_l, obs_B, meanvar_l, cnt_obs_l)

END SUBROUTINE init_obsvar_l_pdafomi



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
  USE obs_A_pdafomi, ONLY: obs_A => thisobs
  USE obs_B_pdafomi, ONLY: obs_B => thisobs

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

  CALL PDAFomi_prodRinvA(obs_A, ncol, A_p, C_p)
  CALL PDAFomi_prodRinvA(obs_B, ncol, A_p, C_p)
  
END SUBROUTINE prodRinvA_pdafomi

