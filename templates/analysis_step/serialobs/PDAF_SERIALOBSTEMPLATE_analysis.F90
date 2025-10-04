!> Analysis routine for DA method SERIALOBSTEMPLATE
!!
!! This routine computes the analysis update of the DA method.
!! Thus, for ensemble DA, it transforms the forecast ensemble
!! into the analysis ensemble.
!!
!! ADAPTING THE TEMPLATE:
!! This template contains the typical steps of the EAKF.
!! On this basis one can implement another DA method. Below we 
!! describe the steps that are included in this code template.
!!
!! __Revision history:__
!! * 2025-10 - Lars Nerger - Initial code for template based on ENSRF
!! * Later revisions - see repository log
!!
MODULE PDAF_SERIALOBSTEMPLATE_analysis

CONTAINS
  SUBROUTINE PDAF_SERIALOBSTEMPLATE_ana(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, &
       HX_p, HXbar_p, obs_p, var_obs_p, &
       U_localize_covar_serial, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_timer, &                ! Routines for timings
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &          ! Routine for memory counting
         ONLY: PDAF_memcount
    USE PDAF_mod_parallel, &         ! Variables for parallelization
         ONLY: mype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step         !< Current time step
    INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_obs_p    !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p)           !< PE-local ensemble mean state
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
    REAL, INTENT(inout) :: HX_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
    REAL, INTENT(inout) :: HXbar_p(dim_obs_p)       !< PE-local observed state
    REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
    REAL, INTENT(in)    :: var_obs_p(dim_obs_p)     !< PE-local vector of observation eror variances
    INTEGER, INTENT(in) :: screen       !< Verbosity flag
    INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
    INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_localize_covar_serial !< Apply localization for single-observation vectors

! *** local variables ***
    INTEGER :: iobs, member             ! Counters
    REAL :: invdim_ensm1                ! inverse of ensemble size minus 1
    INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
    REAL :: HXbar_i                     ! mean observed ensemble for single observation
    REAL :: var_hx                      ! variance of observed ensemble for single observation 
    REAL :: var_ratio                   ! ratio of variances
    REAL, ALLOCATABLE :: HXpert_i(:)    ! observed ensemble perturbations for single observation
    REAL, ALLOCATABLE :: HXinc_i(:)     ! ensemble of observation increments for single observation
    REAL, ALLOCATABLE :: cov_xy_p(:)    ! covariances between state and single observation
    REAL, ALLOCATABLE :: cov_hxy_p(:)   ! covariances between full observed state and single observation
    REAL :: dummy                       ! dummy variable


! **********************
! *** INITIALIZATION ***
! **********************

    CALL PDAF_timeit(51, 'new')

    ! Dummy initialization to present compiler warning
    dummy = state_p(1)

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, i7, 3x, a)') &
            'PDAF ', step, 'SERIALOBSTEMPLATE with serial observation processing'
    END IF

    ! init numbers
    invdim_ensm1 = 1.0 / (REAL(dim_ens - 1))

    ! Allocate arrays
    ALLOCATE(HXpert_i(dim_ens))
    ALLOCATE(HXinc_i(dim_ens))
    ALLOCATE(cov_xy_p(dim_p))
    ALLOCATE(cov_hxy_p(dim_obs_p))
    IF (allocflag == 0) &
         CALL PDAF_memcount(3, 'r', 2*dim_p + 2*dim_ens)

    CALL PDAF_timeit(51, 'old')


! *************************************************************
! *** Loop over all single observations and update ensemble ***
! *************************************************************

    seqObs: DO iobs = 1, dim_obs_p

       ! ***********************************************
       ! *** Preparations                            ***
       ! ***********************************************

       CALL PDAF_timeit(10, 'new')
       CALL PDAF_timeit(51, 'new')
       CALL PDAF_timeit(30, 'new')

       ! Get mean of observed ensemble for single observation
       HXbar_i = HXbar_p(iobs)

       ! Get observed ensemble perturbations for single observation
       HXpert_i(:) = HX_p(iobs,:) - HXbar_i

       ! Compute variance of observed ensemble for single observation
       var_hx = 0.0
       DO member = 1, dim_ens
          var_hx = var_hx + HXpert_i(member)**2
       END DO
       var_hx = var_hx * invdim_ensm1

       ! Compute ratio of variances
       var_ratio = var_obs_p(iobs) / (var_hx + var_obs_p(iobs))

       CALL PDAF_timeit(30, 'old')
       CALL PDAF_timeit(31, 'new')

       ! Compute covariances between state ensemble and single observation
       cov_xy_p = 0.0
       DO member = 1, dim_ens
          cov_xy_p(:) = cov_xy_p(:) + ens_p(:, member) * HXpert_i(member)
       END DO
       cov_xy_p = cov_xy_p * invdim_ensm1

       ! Compute covariances between observed state ensemble and single observation
       cov_hxy_p = 0.0
       DO member = 1, dim_ens
          cov_hxy_p(:) = cov_hxy_p(:) + HX_p(:, member) * HXpert_i(member)
       END DO
       cov_hxy_p = cov_hxy_p * invdim_ensm1

       CALL PDAF_timeit(31, 'old')
       CALL PDAF_timeit(51, 'old')


       ! ****************************************************
       ! *** Apply localization to covariance arrays      ***
       ! ****************************************************

! +++ TEMPLATE: Serial observation processing methods usually apply localization
! +++ to the covariance matrices PH^T and HPH^T. The localization and handling
! +++ of HPH^T is needed for the parallelization.

       CALL PDAF_timeit(45, 'new')
       CALL U_localize_covar_serial(iobs, dim_p, dim_obs_p, cov_xy_p, cov_hxy_p)
       CALL PDAF_timeit(45, 'old')

       CALL PDAF_timeit(10, 'old')


       ! ********************************************************
       ! *** Update step 1: Compute the observation increment ***
       ! ********************************************************

       CALL PDAF_timeit(51, 'new')
       CALL PDAF_timeit(12, 'new')

       ! Update observed ensemble mean
       HXbar_i = var_ratio * HXbar_i + (1.0 - var_ratio) * obs_p(iobs)

       ! Update observed ensemble perturbations for iobs
       HXpert_i(:) = SQRT(var_ratio) * HXpert_i(:)

       ! Complete computation of ensemble of observation increments
       HXinc_i(:) = HXbar_i + HXpert_i(:) - HX_p(iobs, :)

       CALL PDAF_timeit(12, 'new')


       ! ****************************************************
       ! *** Step 2: Update state ensemble                ***
       ! ****************************************************

       CALL PDAF_timeit(14, 'new')

       ! store ratio of covariances in cov_xy_p
       cov_xy_p(:) = cov_xy_p(:) / var_hx

       ! Update ensemble members
       DO member = 1, dim_ens
          ens_p(:, member) = ens_p(:, member) + cov_xy_p(:) * HXinc_i(member) 
       END DO

       CALL PDAF_timeit(14, 'old')


       ! ******************************************************
       ! *** Step 2b: Update observed ensemble and its mean ***
       ! *** This step is required for parallelization      ***
       ! ******************************************************

       CALL PDAF_timeit(13, 'new')

       ! store ratio of covariances in cov_hxy_p
       cov_hxy_p(:) = cov_hxy_p(:) / var_hx

       ! Update observed ensemble members
       DO member = 1, dim_ens
          HX_p(:, member) = HX_p(:, member) + cov_hxy_p(:) * HXinc_i(member) 
       END DO

       ! Compute updated ensemble mean
       HXbar_p = 0.0
       DO member = 1, dim_ens
          HXbar_p(:) = HXbar_p(:) + HX_p(:, member)
       END DO
       HXbar_p = HXbar_p / REAL(dim_ens)

       CALL PDAF_timeit(13, 'old')
       CALL PDAF_timeit(51, 'old')

    END DO seqObs


! ********************
! *** Finishing up ***
! ********************

    ! Clean up
    DEALLOCATE(HXpert_i, HXinc_i)
    DEALLOCATE(cov_xy_p, cov_hxy_p)

! +++ TEMPLATE: Below is generic operation that is required
! +++ memory counting work

    IF (allocflag == 0) allocflag = 1

  END SUBROUTINE PDAF_SERIALOBSTEMPLATE_ana

END MODULE PDAF_SERIALOBSTEMPLATE_analysis
