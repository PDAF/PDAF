!$Id: mod_obs_A_pdaf.F90 251 2019-11-19 08:43:39Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_interface_pdafomi
!
! !DESCRIPTION:
! This module provides interface routines between the call-back routines
! of PDAF and the observation-specific routines in PDAF-OMI. This structure
! collects all calls to observation-specifc routines in this single file
! and allows to make the call-back routines for the observations fully
! generic.
!
! When adding an observation type, one has to add one module
! mod_obs_X_pdaf (based on the template mod_obs_pdaf_TEMPLATE.F90).
! In addition one has to add a call to the different routines include
! in this module.
! 
! !REVISION HISTORY:
! 2019-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, init_dim_obs_f_A, obs_op_f_A, deallocate_obs_A, &
       init_obs_f_A, init_obsvar_A, init_dim_obs_l_A, init_obs_l_A, &
       g2l_obs_A, prodRinvA_l_A, init_obsvar_l_A, rms_obs_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, init_dim_obs_f_B, obs_op_f_B, deallocate_obs_B, &
       init_obs_f_B, init_obsvar_B, init_dim_obs_l_B, init_obs_l_B, &
       g2l_obs_B, prodRinvA_l_B, init_obsvar_l_B, rms_obs_B
  USE mod_obs_C_pdaf, &
       ONLY: assim_C, init_dim_obs_f_C, obs_op_f_C, deallocate_obs_C, &
       init_obs_f_C, init_obsvar_C, init_dim_obs_l_C, init_obs_l_C, &
       g2l_obs_C, prodRinvA_l_C, init_obsvar_l_C, rms_obs_C

  IMPLICIT NONE


!-------------------------------------------------------------------------------

CONTAINS
!BOP
!
! !ROUTINE: init_dim_obs_f_pdafomi --- Interface routine for init_dim_obs_f
!
! !INTERFACE:
  SUBROUTINE init_dim_obs_f_pdafomi(step, dim_obs_f)

! !DESCRIPTION:
! Interface routine for init_dim_obs_f.
! It is called by the call-back routine for init_dim_obs_f
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in)  :: step      ! Current time step
    INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector
!EOP

! *** Local variables
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

    ! Call observation specific routines
    ! The routines are independent, so it is not relevant
    ! in which order they are called
    IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
    IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)
    IF (assim_C) CALL init_dim_obs_f_C(step, dim_obs_f_C)

    dim_obs_f = dim_obs_f_A + dim_obs_f_B + dim_obs_f_C

  END SUBROUTINE init_dim_obs_f_pdafomi


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: obs_op_f_pdafomi --- Interface routine for obs_op_f
!
! !INTERFACE:
  SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! Interface routine for obs_op_f.
! It is called by the call-back routine for obs_op_f.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: step                 ! Current time step
    INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
    INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
    REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
    REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state
!EOP

! *** local variables ***
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
    IF (assim_C) CALL obs_op_f_C(dim_p, dim_obs_f, state_p, m_state_f, offset_obs_f)

  END SUBROUTINE obs_op_f_pdafomi


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: deallocate_obs_pdafomi --- Interface routine for deallocate_obs
!
! !INTERFACE:
  SUBROUTINE deallocate_obs_pdafomi(step)

! !DESCRIPTION:
! Interface routine for deallocate_obs.
! It is called by the call-back routine prepoststep_pdaf.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: step                 ! Current time step
!EOP


! *************************************
! *** Deallocate observation arrays ***
! *************************************

    IF (step > 0) THEN 
       CALL deallocate_obs_A()
       CALL deallocate_obs_B()
       CALL deallocate_obs_C()
    END IF

  END SUBROUTINE deallocate_obs_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obs_f_pdafomi --- Interface routine for init_obs_f
!
! !INTERFACE:
  SUBROUTINE init_obs_f_pdafomi(step, dim_obs_f, observation_f)

! !DESCRIPTION:
! Interface routine for init_obs_f.
! It is called by the call-back routine for init_obs_f.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: step        ! Current time step
    INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
    REAL, INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector
!EOP

! *** local variables ***
    INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    offset_obs_f = 0

    ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
    IF (assim_A) CALL init_obs_f_A(dim_obs_f, observation_f, offset_obs_f)
    IF (assim_B) CALL init_obs_f_B(dim_obs_f, observation_f, offset_obs_f)
    IF (assim_C) CALL init_obs_f_C(dim_obs_f, observation_f, offset_obs_f)

  END SUBROUTINE init_obs_f_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obsvar_pdafomi --- Interface routine for init_obsvar
!
! !INTERFACE:
  SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)

! !DESCRIPTION:
! Interface routine for init_obsvar.
! It is called by the call-back routine for init_obsvar.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: step          ! Current time step
    INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
    REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
    REAL, INTENT(out)   :: meanvar       ! Mean observation error variance
!EOP

! *** Local variables
    INTEGER :: cnt_obs_f


! *****************************
! *** Compute mean variance ***
! *****************************

    ! Initialize observation counter
    cnt_obs_f = 0

    ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
    IF (assim_A) CALL init_obsvar_A(meanvar, cnt_obs_f)
    IF (assim_B) CALL init_obsvar_B(meanvar, cnt_obs_f)
    IF (assim_C) CALL init_obsvar_C(meanvar, cnt_obs_f)

  END SUBROUTINE init_obsvar_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_dim_obs_l_pdafomi --- Interface routine for init_dim_obs_l
!
! !INTERFACE:
  SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, coords_l, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! Interface routine for init_dim_obs_l.
! It is called by the call-back routine for init_dim_obs_l.
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
    INTEGER, INTENT(in)  :: step       ! Current time step
    REAL, INTENT(in) :: coords_l(:)    ! Coordinates of local analysis domain
    INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
    INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector
!EOP


! *** local variables ***
    INTEGER :: dim_obs_l_A ! Dimension of observation type A
    INTEGER :: dim_obs_l_B ! Dimension of observation type B
    INTEGER :: dim_obs_l_C ! Dimension of observation type C
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
    dim_obs_l_C = 0

    ! Call init_dim_obs_l specific for each observation
    ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
    IF (assim_A) CALL init_dim_obs_l_A(coords_l, local_range, dim_obs_l_A, &
         offset_obs_l, offset_obs_f)
    IF (assim_B) CALL init_dim_obs_l_B(coords_l, local_range, dim_obs_l_B, &
         offset_obs_l, offset_obs_f)
    IF (assim_C) CALL init_dim_obs_l_C(coords_l, local_range, dim_obs_l_C, &
         offset_obs_l, offset_obs_f)

    ! Compute overall local observation dimension
    dim_obs_l = dim_obs_l_A + dim_obs_l_B + dim_obs_l_C

  END SUBROUTINE init_dim_obs_l_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obs_l_pdafomi --- Interface routine for init_obs_l
!
! !INTERFACE:
  SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! Interface routine for init_obs_l.
! It is called by the call-back routine for init_obs_l.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain index
    INTEGER, INTENT(in) :: step       ! Current time step
    INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
    REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector
!EOP


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    IF (assim_A) CALL init_obs_l_A(dim_obs_l, observation_l)
    IF (assim_B) CALL init_obs_l_B(dim_obs_l, observation_l)
    IF (assim_C) CALL init_obs_l_C(dim_obs_l, observation_l)

  END SUBROUTINE init_obs_l_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g2l_obs_pdafomi --- Interface routine for g2l_obs_pdafomi
!
! !INTERFACE:
  SUBROUTINE g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, &
       mstate_l)

! !DESCRIPTION:
! Interface routine for g2l_obs.
! It is called by the call-back routine for g2l_obs.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain
    INTEGER, INTENT(in) :: step       ! Current time step
    INTEGER, INTENT(in) :: dim_obs_f  ! Dimension of full PE-local obs. vector
    INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
    REAL, INTENT(in)    :: mstate_f(dim_obs_f)   ! Full PE-local obs. vector
    REAL, INTENT(out)   :: mstate_l(dim_obs_l)   ! Obs. vector on local domain
!EOP


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

    IF (assim_A) CALL g2l_obs_A(dim_obs_l, dim_obs_f, mstate_f, mstate_l)
    IF (assim_B) CALL g2l_obs_B(dim_obs_l, dim_obs_f, mstate_f, mstate_l)
    IF (assim_C) CALL g2l_obs_C(dim_obs_l, dim_obs_f, mstate_f, mstate_l)

  END SUBROUTINE g2l_obs_pdafomi



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: prodRinvA_l_pdafomi --- Interface routine for prodRinvA_l
!
! !INTERFACE:
  SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

! !DESCRIPTION:
! Interface routine for prodRinvA_l.
! It is called by the call-back routine for prodRinvA_l.
!
! !USES:
    USE mod_assimilation, &
         ONLY: local_range, locweight, srange

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
    INTEGER, INTENT(in) :: step              ! Current time step
    INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
    INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
    REAL, INTENT(in)    :: obs_l(dim_obs_l)  ! Local vector of observations
    REAL, INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix
    REAL, INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix
!EOP


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

    ! Screen output
    IF (verbose == 1) THEN
       IF (assim_A) WRITE (*, '(8x, a, f12.3)') &
            '--- Use global rms for observations A of ', rms_obs_A
       IF (assim_B) WRITE (*, '(8x, a, f12.3)') &
            '--- Use global rms for observations B of ', rms_obs_B
       IF (assim_C) WRITE (*, '(8x, a, f12.3)') &
            '--- Use global rms for observations C of ', rms_obs_C
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
    IF (assim_C) CALL prodRinvA_l_C(verbose, dim_obs_l, rank, locweight, local_range, &
         srange, A_l, C_l)
  
  END SUBROUTINE prodRinvA_l_pdafomi



!-------------------------------------------------------------------------------
!BOP
! !ROUTINE: init_obsvar_l_pdafomi --- Interface routine for init_obsvar_l
!
! !INTERFACE:
  SUBROUTINE init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)

! !DESCRIPTION:
! Interface routine for init_obsvar_l.
! It is called by the call-back routine for init_obsvar_l.
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: domain_p      ! Current local analysis domain
    INTEGER, INTENT(in) :: step          ! Current time step
    INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
    REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
    REAL, INTENT(out)   :: meanvar_l     ! Mean local observation error variance
!EOP

! *** Local variables
    INTEGER :: cnt_obs_l


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    ! Initialize observation counter
    cnt_obs_l = 0

    IF (assim_A) CALL init_obsvar_l_A(meanvar_l, cnt_obs_l)
    IF (assim_B) CALL init_obsvar_l_B(meanvar_l, cnt_obs_l)
    IF (assim_C) CALL init_obsvar_l_C(meanvar_l, cnt_obs_l)

  END SUBROUTINE init_obsvar_l_pdafomi

END MODULE mod_interface_pdafomi
