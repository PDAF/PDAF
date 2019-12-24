!$Id$
!> \brief PDAF-OMI interface module
!!
!! \details This module provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specifc routines in this single file
!! and allows to make the call-back routines for the observations fully
!! generic.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! mod_obs_X_pdaf (based on the template mod_obs_pdaf_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this module.
!! 
!! \date 2019-12 - Lars Nerger - Initial code
!!
MODULE interface_pdafomi

  USE obs_A_pdafomi, &
       ONLY: assim_A, init_dim_obs_f_A, obs_op_f_A, deallocate_obs_A, &
       init_obs_f_A, init_obsvar_A, init_dim_obs_l_A, init_obs_l_A, &
       g2l_obs_A, prodRinvA_l_A, init_obsvar_l_A, rms_obs_A
  USE obs_B_pdafomi, &
       ONLY: assim_B, init_dim_obs_f_B, obs_op_f_B, deallocate_obs_B, &
       init_obs_f_B, init_obsvar_B, init_dim_obs_l_B, init_obs_l_B, &
       g2l_obs_B, prodRinvA_l_B, init_obsvar_l_B, rms_obs_B

  IMPLICIT NONE


!-------------------------------------------------------------------------------

CONTAINS

!> \brief Interface routine for init_dim_obs_f
!!
!! \details This routine calls the observation-specific
!! routines init_dim_obs_f_X.
!! It is called by the call-back routine for init_dim_obs_f.
!!
  SUBROUTINE init_dim_obs_f_pdafomi(step, dim_obs_f)

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

    ! Call observation specific routines
    ! The routines are independent, so it is not relevant
    ! in which order they are called
    IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
    IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)

    dim_obs_f = dim_obs_f_A + dim_obs_f_B

  END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!> \brief Interface routine for obs_op_f
!!
!! \details This routine calls the observation-specific
!! routines obs_op_f_X.
!! It is called by the call-back routine for obs_op_f.
!!
  SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, m_state_f)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step                 !< Current time step
    INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs_f            !< Dimension of full observed state
    REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
    REAL, INTENT(inout) :: m_state_f(dim_obs_f) !< PE-local full observed state

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
    IF (assim_A) CALL obs_op_f_A(dim_p, dim_obs_f, state_p, m_state_f, offset_obs_f)
    IF (assim_B) CALL obs_op_f_B(dim_p, dim_obs_f, state_p, m_state_f, offset_obs_f)

  END SUBROUTINE obs_op_f_pdafomi


!-------------------------------------------------------------------------------
!> \brief Interface routine for deallocate_obs
!!
!! \details This routine calls the observation-specific
!! routines deallocate_obs_X.
!! It is called by the call-back routine prepoststep_pdaf.
!!
  SUBROUTINE deallocate_obs_pdafomi(step)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

    IF (step > 0) THEN 
       CALL deallocate_obs_A()
       CALL deallocate_obs_B()
    END IF

  END SUBROUTINE deallocate_obs_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for init_obs_f
!!
!! \details This routine calls the observation-specific
!! routines init_obs_f_X.
!! It is called by the call-back routine for init_obs_f.
!!
  SUBROUTINE init_obs_f_pdafomi(step, dim_obs_f, observation_f)

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
    IF (assim_A) CALL init_obs_f_A(dim_obs_f, observation_f, offset_obs_f)
    IF (assim_B) CALL init_obs_f_B(dim_obs_f, observation_f, offset_obs_f)

  END SUBROUTINE init_obs_f_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for init_obsvar
!!
!! \details This routine calls the observation-specific
!! routines init_obsvar_X.
!! It is called by the call-back routine for init_obsvar.
!!
  SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)

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
    IF (assim_A) CALL init_obsvar_A(meanvar, cnt_obs_f)
    IF (assim_B) CALL init_obsvar_B(meanvar, cnt_obs_f)

  END SUBROUTINE init_obsvar_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for init_dim_obs_l
!!
!! \details This routine calls the observation-specific
!! routines init_dim_obs_l_X.
!! It is called by the call-back routine for init_dim_obs_l.
!!
  SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, coords_l, dim_obs_f, dim_obs_l)

    USE mod_assimilation, &   
         ONLY: local_range    !< localization radius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step       !< Current time step
    REAL, INTENT(in) :: coords_l(:)    !< Coordinates of local analysis domain
    INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
    INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
    INTEGER :: dim_obs_l_A ! Dimension of observation type A
    INTEGER :: dim_obs_l_B ! Dimension of observation type B
    INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Initialize offsets with zero
    offset_obs_l = 0
    offset_obs_f = 0

    ! Initialize local dimensions
    dim_obs_l_A = 0
    dim_obs_l_B = 0

    ! Call init_dim_obs_l specific for each observation
    ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
    IF (assim_A) CALL init_dim_obs_l_A(coords_l, local_range, dim_obs_l_A, &
         offset_obs_l, offset_obs_f)
    IF (assim_B) CALL init_dim_obs_l_B(coords_l, local_range, dim_obs_l_B, &
         offset_obs_l, offset_obs_f)

    ! Compute overall local observation dimension
    dim_obs_l = dim_obs_l_A + dim_obs_l_B

  END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for init_obs_l
!!
!! \details This routine calls the observation-specific
!! routines init_obs_l_X.
!! It is called by the call-back routine for init_obs_l.
!!
  SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
    INTEGER, INTENT(in) :: step       !< Current time step
    INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
    REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    IF (assim_A) CALL init_obs_l_A(dim_obs_l, observation_l)
    IF (assim_B) CALL init_obs_l_B(dim_obs_l, observation_l)

  END SUBROUTINE init_obs_l_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for g2l_obs
!!
!! \details This routine calls the observation-specific
!! routines g2l_obs_X.
!! It is called by the call-back routine for g2l_obs.
!!
  SUBROUTINE g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, &
       mstate_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain
    INTEGER, INTENT(in) :: step       !< Current time step
    INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local observation vector
    INTEGER, INTENT(in) :: dim_obs_l  !< Dimension of local observation vector
    REAL, INTENT(in)    :: mstate_f(dim_obs_f)   !< Full PE-local obs.ervation vector
    REAL, INTENT(out)   :: mstate_l(dim_obs_l)   !< Observation vector on local domain


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

    IF (assim_A) CALL g2l_obs_A(dim_obs_l, dim_obs_f, mstate_f, mstate_l)
    IF (assim_B) CALL g2l_obs_B(dim_obs_l, dim_obs_f, mstate_f, mstate_l)

  END SUBROUTINE g2l_obs_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for prodRinvA_l
!!
!! \details This routine calls the observation-specific
!! routines prodRinvA_l_X.
!! It is called by the call-back routine for prodRinvA_l.
!!
  SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

    USE mod_assimilation, &
         ONLY: local_range, locweight, srange

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
    INTEGER :: verbose       ! Verbosity flag
    INTEGER, SAVE :: domain_save = -1  ! Save previous domain index


! **********************
! *** INITIALIZATION ***
! **********************

    IF (domain_p <= domain_save .OR. domain_save < 0) THEN
       verbose = 1
    ELSE
       verbose = 0
    END IF
    domain_save = domain_p


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

    IF (assim_A) CALL prodRinvA_l_A(verbose, dim_obs_l, rank, locweight, local_range, &
         srange, A_l, C_l)
    IF (assim_B) CALL prodRinvA_l_B(verbose, dim_obs_l, rank, locweight, local_range, &
         srange, A_l, C_l)
  
  END SUBROUTINE prodRinvA_l_pdafomi



!-------------------------------------------------------------------------------
!> \brief Interface routine for init_obsvar_l
!!
!! \details This routine calls the observation-specific
!! routines init_obsvar_l_X.
!! It is called by the call-back routine for init_obsvar_l.
!!
  SUBROUTINE init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)

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

    IF (assim_A) CALL init_obsvar_l_A(meanvar_l, cnt_obs_l)
    IF (assim_B) CALL init_obsvar_l_B(meanvar_l, cnt_obs_l)

  END SUBROUTINE init_obsvar_l_pdafomi

END MODULE interface_pdafomi
