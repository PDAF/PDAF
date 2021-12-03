C     $Id: callback_obs_pdafomi.F90 251 2019-11-19 08:43:39Z lnerger $
C     !> callback_obs_pdafomi

C     ! This file provides interface routines between the call-back routines
C     ! of PDAF and the observation-specific routines in PDAF-OMI. This structure
C     ! collects all calls to observation-specific routines in this single file
C     ! to make it easier to find the routines that need to be adapted.
C     !
C     ! The routines here are mainly pure pass-through routines. Thus they
C     ! simply call one of the routines from PDAF-OMI. Partly some addtional
C     ! variable is required, e.g. to specify the offset of an observation
C     ! in the observation vector containing all observation types. These
C     ! cases are described in the routines.
C     !
C     ! **Adding an observation type:**
C     ! When adding an observation type, one has to add one module
C     ! obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
C     ! In addition one has to add a call to the different routines include
C     ! in this file. It is recommended to keep the order of the calls
C     ! consistent over all files.
C     !
C     ! __Revision history:__
C     ! * 2021-03 - Nambir Mamnun - Initial code
C     ! * Later revisions - see repository log
C     !
C     !-------------------------------------------------------------------------------

C     !> Call-back routine for init_dim_obs
C     !
C     ! This routine calls the observation-specific
C     ! routines init_dim_obs_TYPE.
C     !
      SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

C     ! Include functions for different observations
      USE obs_etan_pdafomi, ONLY: assim_etan, init_dim_obs_etan

      IMPLICIT NONE

C     ! *** Arguments ***
      INTEGER, INTENT(in)  :: step     !< Current time step
      INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

C     ! *** Local variables ***
      INTEGER :: dim_obs_etan ! Observation dimensions


C     ! *********************************************
C     ! *** Initialize full observation dimension ***
C     ! *********************************************

C     ! Initialize number of observations

      dim_obs_etan = 0

C     ! Call observation-specific routines
C     ! The routines are independent, so it is not relevant
C     ! in which order they are called

      IF (assim_etan) CALL init_dim_obs_etan(step, dim_obs_etan)

      dim_obs = dim_obs_etan

      END


C     !-------------------------------------------------------------------------------
C     !> Call-back routine for obs_op
C     !
C     ! This routine calls the observation-specific
C     ! routines obs_op_TYPE.
C     !
      SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

C     ! Include functions for different observations
      USE obs_etan_pdafomi, ONLY: obs_op_etan

      IMPLICIT NONE

C! *** Arguments ***
      INTEGER, INTENT(in) :: step                    !< Current time step
      INTEGER, INTENT(in) :: dim_p                   !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs                 !< Dimension of full observed state
      REAL(kind=8), INTENT(in)    :: state_p(dim_p)  !< PE-local model state
      REAL(kind=8), INTENT(inout) :: ostate(dim_obs) !< PE-local full observed state


C     ! ******************************************************
C     ! *** Apply observation operator H on a state vector ***
C     ! ******************************************************

C     ! The order of these calls is not relevant as the setup
C     ! of the overall observation vector is defined by the
C     ! order of the calls in init_dim_obs_pdafomi
      CALL obs_op_etan(dim_p, dim_obs, state_p, ostate)
C     !  CALL obs_op_etan2(dim_p, dim_obs, state_p, ostate)
C     !  CALL obs_op_C(dim_p, dim_obs, state_p, ostate)

      END



C     !-------------------------------------------------------------------------------
C     !> Call-back routine for init_dim_obs_l
C     !
C     ! This routine calls the routine PDAFomi_init_dim_obs_l
C     ! for each observation type
C     !
      SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs,
     &    dim_obs_l)

C     ! Include functions for different observations
      USE obs_etan_pdafomi, ONLY: init_dim_obs_l_etan
      
      IMPLICIT NONE

C     ! *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step       !< Current time step
      INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
      INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


C     ! **********************************************
C     ! *** Initialize local observation dimension ***
C     ! **********************************************

C     ! Call init_dim_obs_l specific for each observation
      CALL init_dim_obs_l_etan(domain_p, step, dim_obs, dim_obs_l)

      END



C     !-------------------------------------------------------------------------------
C     !> Call-back routine for localize_covar
C     !
C     ! This routine calls the routine PDAFomi_localize_covar
C     ! for each observation type to apply covariance
C     ! localization in the LEnKF.
C     !
      SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)


C     ! Include functions for different observations
      USE obs_etan_pdafomi, ONLY: localize_covar_etan

      IMPLICIT NONE

#include "CPP_OPTIONS.h"
#include "SIZE.h"
#include "GRID.h"

C     ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                         !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs                       !< number of observations
      REAL(kind=8), INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL(kind=8), INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

C! *** local variables ***
      INTEGER :: i, j, cnt                       ! Counters
      REAL(kind=8), ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries


C     ! **********************
C     ! *** INITIALIZATION ***
C     ! **********************

C     ! Initialize coordinate array

      ALLOCATE(coords_p(2, dim_p))

      cnt = 0
      DO j = 1, sNx
        DO i= 1, sNy
          cnt = cnt + 1
          coords_p(1,cnt) = xC(j, i, 1, 1)
          coords_p(2,cnt) = yC(j, i, 1, 1)
        END DO
      END DO


C     ! *************************************
C     ! *** Apply covariance localization ***
C     ! *************************************

C     ! Call localize_covar specific for each observation
      CALL localize_covar_etan(dim_p, dim_obs, HP_p, HPH, coords_p)


C     ! ****************
C     ! *** Clean up ***
C     ! ****************

      DEALLOCATE(coords_p)

      END
