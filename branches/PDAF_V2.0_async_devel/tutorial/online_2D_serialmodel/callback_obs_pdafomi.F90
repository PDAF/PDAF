!$Id: mod_obs_A_pdaf.F90 251 2019-11-19 08:43:39Z lnerger $
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

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: assim_A, async_A, init_dim_obs_A, init_dim_obs_async_A
  USE obs_B_pdafomi, ONLY: assim_B, async_B, init_dim_obs_B, init_dim_obs_async_B
  USE obs_C_pdafomi, ONLY: assim_C, init_dim_obs_C

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_A ! Observation dimensions
  INTEGER :: dim_obs_B ! Observation dimensions
  INTEGER :: dim_obs_C ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_A = 0
  dim_obs_B = 0
  dim_obs_C = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A) THEN
     IF (async_A) THEN
        CALL init_dim_obs_async_A(step, dim_obs_A)
     ELSE
        CALL init_dim_obs_A(step, dim_obs_A)
     END IF
  END IF
  IF (assim_B) THEN
     IF (async_B) THEN
        CALL init_dim_obs_async_B(step, dim_obs_B)
     ELSE
        CALL init_dim_obs_B(step, dim_obs_B)
     END IF
  END IF

  IF (assim_C) CALL init_dim_obs_C(step, dim_obs_C)

  dim_obs = dim_obs_A + dim_obs_B + dim_obs_C

END SUBROUTINE init_dim_obs_pdafomi


!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: obs_op_A
  USE obs_B_pdafomi, ONLY: obs_op_B
  USE obs_C_pdafomi, ONLY: obs_op_C

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

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi
  CALL obs_op_A(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_B(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_C(dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: init_dim_obs_l_A
  USE obs_B_pdafomi, ONLY: init_dim_obs_l_B
  USE obs_C_pdafomi, ONLY: init_dim_obs_l_C
  
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Call init_dim_obs_l specific for each observation
  CALL init_dim_obs_l_A(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_B(domain_p, step, dim_obs, dim_obs_l)
  CALL init_dim_obs_l_C(domain_p, step, dim_obs, dim_obs_l)

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
  USE obs_A_pdafomi, ONLY: localize_covar_A
  USE obs_B_pdafomi, ONLY: localize_covar_B
  USE obs_C_pdafomi, ONLY: localize_covar_C

  ! Include information on model grid
  USE mod_model, &
       ONLY: nx, ny

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs               !< number of observations
  REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

! *** local variables ***
  INTEGER :: i, j, cnt               ! Counters
  REAL, ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries


! **********************
! *** INITIALIZATION ***
! **********************

  ! Initialize coordinate array

  ALLOCATE(coords_p(2, dim_p))

  cnt = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt = cnt + 1
        coords_p(1, cnt) = REAL(j)
        coords_p(2, cnt) = REAL(i)
     END DO
  END DO


! *************************************
! *** Apply covariance localization ***
! *************************************

  ! Call localize_covar specific for each observation
  CALL localize_covar_A(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_B(dim_p, dim_obs, HP_p, HPH, coords_p)
  CALL localize_covar_C(dim_p, dim_obs, HP_p, HPH, coords_p)


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(coords_p)

END SUBROUTINE localize_covar_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_async
!!
!! This routine is executed at for asynchronous
!! DA at the beginning of a forecast phase. 
!! The routine calls the observation-specific
!! routines init_dim_obs_TYPE to initialize
!! the observational information.
!!
SUBROUTINE init_dim_obs_async_pdafomi(step)

  USE mod_assimilation, ONLY: dim_obs_all

  ! Include functions for different observations
  USE obs_A_pdafomi, ONLY: assim_A, async_A, init_dim_obs_A
  USE obs_B_pdafomi, ONLY: assim_B, async_B, init_dim_obs_B
  USE obs_C_pdafomi, ONLY: assim_C, init_dim_obs_C

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step

! *** Local variables ***
  INTEGER :: dim_obs_A ! Observation dimensions
  INTEGER :: dim_obs_B ! Observation dimensions
  INTEGER :: dim_obs_C ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_A = 0
  dim_obs_B = 0
  dim_obs_C = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A .AND. async_A) CALL init_dim_obs_A(step, dim_obs_A)
  IF (assim_B .AND. async_B) CALL init_dim_obs_B(step, dim_obs_B)
!  IF (assim_C) CALL init_dim_obs_C(step, dim_obs_C)

  ! Store overall number of full observations
  dim_obs_all = dim_obs_A + dim_obs_B + dim_obs_C

END SUBROUTINE init_dim_obs_async_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_async
!!
!! For asynchronous DA this routine first counts
!! the observations valid at the current time
!! step. Subsequently, iy calls the asynchronous
!! observation operators obs_op_async_TYPE.
!!
SUBROUTINE obs_op_async_pdafomi(step)

  USE PDAF_interfaces_module, ONLY: PDAF_set_ens_pointer
  USE mod_assimilation, ONLY: dim_state_p
  USE mod_parallel_pdaf, ONLY: mype_world
  USE obs_A_pdafomi, ONLY: cnt_obs_async_A, obs_op_A_async
  USE obs_B_pdafomi, ONLY: cnt_obs_async_B, obs_op_B_async
  
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step

! *** Local variables
  INTEGER :: nobs                     ! Number of observations
  INTEGER :: status                   ! PDAF status flag
  REAL, POINTER :: ens_pointer(:,:)   ! Pointer to ensemble array


! *****************************************************************************
! *** Asynchronous application of observation operator H on a state vector ***
! *****************************************************************************

  ! *** Count current observations ***
  nobs = 0
  CALL cnt_obs_async_A(step, nobs)
  CALL cnt_obs_async_B(step, nobs)
  if (nobs>0 .and. mype_world==0) write (*,*) 'Number of obs at step ', step, ' =', nobs

  ! *** Apply asynchronous observation operator if we have observations ***
  IF (nobs>0) THEN

     ! Set pointer to ensemble array of PDAF and initialize state vector
     CALL PDAF_set_ens_pointer(ens_pointer, status)
     CALL collect_state_pdaf(dim_state_p, ens_pointer(:,1))

     ! Apply asynchronous process-local observation operator
     CALL obs_op_A_async(step, ens_pointer(:,1))
     CALL obs_op_B_async(step, ens_pointer(:,1))
  END IF

END SUBROUTINE obs_op_async_pdafomi



!-------------------------------------------------------------------------------
!> Gather observed ensemble at analysis time for asynchronous DA
!!
!! For asynchronous DA this routine gathers the 
!! observed ensemble on the filter-PEs.
!! In addition the observation operator is called 
!! by the processes which are not filter-PEs. 
!!
SUBROUTINE gather_obsens_async_pdafomi()
 
  USE obs_A_pdafomi, ONLY: gather_obsens_async_A
  USE obs_B_pdafomi, ONLY: gather_obsens_async_B

  IMPLICIT NONE


! ***********************************************************
! *** Asynchronous DA: Initialize observed state ensemble ***
! ***********************************************************

  CALL gather_obsens_async_A()
  CALL gather_obsens_async_B()

END SUBROUTINE gather_obsens_async_pdafomi
