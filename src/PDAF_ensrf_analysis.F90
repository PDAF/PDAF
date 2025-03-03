! Copyright (c) 2004-2025 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!
!> Analysis step for Kalman fitler with serial observation processing
!!
!! Analysis step of ensemble squre-root filter with 
!! serial observation processing. The two routines 
!! are for two subtypes:
!! * 0: ENSRF of Whitaber and Hamill (2002)
!! * 1: local least squares formulation of Anderson (2003)
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_ensrf_analysis

CONTAINS
!> Perform ENSRF analysis step
!!
!! Analysis step of ensemble square-root filter with 
!! serial observation processing following Whitaker
!! and Hamill (2002).  
!!
  SUBROUTINE PDAF_ensrf_ana(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, HX_p, HXbar_p, obs_p, var_obs_p, &
       U_localize_covar_serial, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &
         ONLY: mype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step          !< Current time step
    INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens       !< Size of state ensemble
    REAL, INTENT(inout) :: state_p(dim_p)           !< PE-local ensemble mean state
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
    REAL, INTENT(inout) :: HX_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
    REAL, INTENT(inout) :: HXbar_p(dim_obs_p)       !< PE-local observed state
    REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
    REAL, INTENT(in)    :: var_obs_p(dim_obs_p)     !< PE-local vector of observation eror variances
    INTEGER, INTENT(in) :: screen        !< Verbosity flag
    INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
    INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_localize_covar_serial  !< Apply localization for single-observation vectors

! *** local variables ***
    INTEGER :: iobs, member              ! Counters
    REAL :: invdim_ensm1                 ! inverse of ensemble size minus 1
    INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
    REAL :: alpha                        ! variance factor for ensemble perturbation update
    REAL :: HPH                          ! Value of HPH^T for single observation
    REAL :: HPHpR                        ! Value of HPH^T+R for single observation              
    REAL :: HXbar_i                      ! mean observed ensemble for single observation
    REAL :: innov_i                      ! innovation for single observation
    REAL, ALLOCATABLE :: Xpert_p(:,:)    ! array of ensemble perturbations
    REAL, ALLOCATABLE :: HXpert_i(:)     ! observed ensemble perturbations for single observation
    REAL, ALLOCATABLE :: HXpert_p(:,:)   ! observed ensemble perturbation for full observations
    REAL, ALLOCATABLE :: HP_p(:)         ! Matrix HP for PE-local state
    REAL, ALLOCATABLE :: HXY_p(:)        ! Matrix HX (HX)^T for all observations


! **********************
! *** INITIALIZATION ***
! **********************

    CALL PDAF_timeit(51, 'new')

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_ensrf_analysis -- START'

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, i7, 3x, a)') &
            'PDAF ', step, 'ENSRF analysis with serial observation processing'
    END IF

    ! init numbers
    invdim_ensm1 = 1.0 / (REAL(dim_ens - 1))

    ALLOCATE(Xpert_p(dim_p, dim_ens))
    ALLOCATE(HXpert_p(dim_obs_p, dim_ens))
    ALLOCATE(HP_p(dim_p))
    ALLOCATE(HXY_p(dim_obs_p))
    ALLOCATE(HXpert_i(dim_ens))
    IF (allocflag == 0) &
         CALL PDAF_memcount(3, 'r', (dim_p+dim_obs_p) * dim_ens + dim_p + dim_obs_p + dim_ens)


! *****************************************************
! *** Initial computation of ensemble perturbations ***
! *****************************************************

    CALL PDAF_timeit(10, 'new')
    CALL PDAF_timeit(30, 'new')

    ! Initialize ensemble perturbations
    ENSa: DO member = 1, dim_ens
       Xpert_p(:, member) = ens_p(:, member) - state_p(:)
    END DO ENSa

    ! Store observed ensemble perturbations
    DO member = 1, dim_ens
       HXpert_p(:, member) = HX_p(:, member) - HXbar_p(:)
    END DO

    CALL PDAF_timeit(30, 'old')
    CALL PDAF_timeit(10, 'old')
    CALL PDAF_timeit(51, 'old')


! *************************************************************
! *** Loop over all single observations and update ensemble ***
! *************************************************************

    seqObs: DO iobs = 1, dim_obs_p

       ! *******************************************************
       ! *** Compute the matrix HP and scalar HPH=HPH^T as   ***
       ! *** ensemble means from ensemble perturbations and  ***
       ! *** the perturbations of the observed ensemble.     ***
       ! *******************************************************

       CALL PDAF_timeit(10, 'new')
       CALL PDAF_timeit(51, 'new')

       ! *** Preparate ensemble perturbations and means ***

       CALL PDAF_timeit(30, 'new')

       ! Get mean of observed ensemble for current observation
       HXbar_i = HXbar_p(iobs)

       ! Get observed ensemble perturbations for current observation
       HXpert_i(:) = HXpert_p(iobs,:)

       CALL PDAF_timeit(30, 'old')


       ! *** Compute HP and HPH ***

       ! HP for model state

       CALL PDAF_timeit(31, 'new')

       HP_p = 0.0
       DO member = 1, dim_ens
          HP_p(:) = HP_p(:) + HXpert_i(member) * Xpert_p(:, member)
       END DO
       HP_p = HP_p * invdim_ensm1

       CALL PDAF_timeit(31, 'old')

       ! HP for observed model state

       CALL PDAF_timeit(32, 'new')

       HXY_p = 0.0
       DO member = 1, dim_ens
          HXY_p(:) = HXY_p(:) + HXpert_i(member) * HXpert_p(:, member)
       END DO
       HXY_p = HXY_p * invdim_ensm1

       CALL PDAF_timeit(32, 'old')

       ! Compute HPH^T

       CALL PDAF_timeit(34, 'new')

       HPH = 0.0
       DO member = 1, dim_ens
          HPH = HPH + HXpert_i(member)**2
       END DO
       HPH = HPH * invdim_ensm1

       CALL PDAF_timeit(34, 'old')


       ! *** Add observation error covariance ***

       HPHpR = HPH + var_obs_p(iobs)

       CALL PDAF_timeit(51, 'old')


       ! ********************************************
       ! *** Apply localization to HP_p and HXY_p ***
       ! ********************************************

       CALL PDAF_timeit(45, 'new')
       CALL U_localize_covar_serial(iobs, dim_p, dim_obs_p, HP_p, HXY_p)
       CALL PDAF_timeit(45, 'old')

       CALL PDAF_timeit(10, 'old')


       ! **************************************
       ! *** Compute innovation d = y - H x ***
       ! **************************************

       CALL PDAF_timeit(51, 'new')
       CALL PDAF_timeit(12, 'new')

       innov_i = obs_p(iobs) - HXbar_i


       ! *******************************************
       ! *** Compute Kalman gain HP (HPH)^-1     ***
       ! *******************************************

       ! For state
       HP_p = HP_p / HPHpR

       ! For observed state
       HXY_p = HXY_p / HPHpR


       ! **********************************************
       ! *** Compute scaling factor for Kalman gain ***
       ! *** (Whitaker/Hamill (2002), Eq. 13)       ***
       ! **********************************************

       alpha = 1.0 + SQRT(var_obs_p(iobs) / HPHpR)
       alpha = 1.0 / alpha

       CALL PDAF_timeit(12, 'old')


       ! ***********************************************
       ! *** Update ensemble mean and state ensemble ***
       ! ***********************************************

       CALL PDAF_timeit(14, 'new')

       ! Update mean state
       state_p(:) = state_p(:) + HP_p(:) * innov_i

       ! Update ensemble members
       DO member = 1, dim_ens
          Xpert_p(:, member) = Xpert_p(:, member) - alpha * HP_p(:) * HXpert_i(member) 
       END DO

       ! re-create full ensemble states
       DO member = 1, dim_ens
          ens_p(:, member) = Xpert_p(:, member) + state_p(:)
       END DO

       CALL PDAF_timeit(14, 'old')


       ! ***********************************************
       ! *** Update observed ensemble and its mean   ***
       ! ***********************************************

       CALL PDAF_timeit(13, 'new')

       ! Update mean observed state
       HXbar_p(:) = HXbar_p(:) + HXY_p(:) * innov_i

       ! Update observed ensemble members
       DO member = 1, dim_ens
          HXpert_p(:, member) = HXpert_p(:, member) - alpha * HXY_p(:) * HXpert_i(member)
       END DO

       ! re-create full observed ensemble states
       DO member = 1, dim_ens
          HX_p(:, member) = HXpert_p(:, member) + HXbar_p(:)
       END DO

       CALL PDAF_timeit(13, 'old')
       CALL PDAF_timeit(51, 'old')

    END DO seqObs


! ********************
! *** Finishing up ***
! ********************

    ! Clean up
    DEALLOCATE(Xpert_p)
    DEALLOCATE(HXpert_i)
    DEALLOCATE(HXpert_p)
    DEALLOCATE(HP_p)
    DEALLOCATE(HXY_p)

    IF (allocflag == 0) allocflag = 1

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_ensrf_analysis -- END'

  END SUBROUTINE PDAF_ensrf_ana


!-------------------------------------------------------------------------------
!> Serial observation processing filter with 2-step formulation 
!!
!! This variant of the serial processing filter uses the EAKF
!! local least squares formulation of Anderson (2003).
!! 1. The increment for the observed model state is computed.
!! 2. The increment is linearly regressed on the model state
!!
  SUBROUTINE PDAF_ensrf_ana_2step(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, HX_p, HXbar_p, obs_p, var_obs_p, &
       U_localize_covar_serial, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &
         ONLY: mype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step          !< Current time step
    INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens       !< Size of state ensemble
    REAL, INTENT(inout) :: state_p(dim_p)           !< PE-local ensemble mean state
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
    REAL, INTENT(inout) :: HX_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
    REAL, INTENT(inout) :: HXbar_p(dim_obs_p)       !< PE-local observed state
    REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
    REAL, INTENT(in)    :: var_obs_p(dim_obs_p)     !< PE-local vector of observation eror variances
    INTEGER, INTENT(in) :: screen        !< Verbosity flag
    INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
    INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_localize_covar_serial  !< Apply localization for single-observation vectors

! *** local variables ***
    INTEGER :: iobs, member              ! Counters
    REAL :: invdim_ensm1                 ! inverse of ensemble size minus 1
    INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
    REAL :: HXbar_i                      ! mean observed ensemble for single observation
    REAL :: var_hx                       ! variance of observed ensemble for single observation 
    REAL :: var_ratio                    ! ratio of variances
    REAL, ALLOCATABLE :: HXpert_i(:)     ! observed ensemble perturbations for single observation
    REAL, ALLOCATABLE :: HXinc_i(:)      ! ensemble of observation increments for single observation
    REAL, ALLOCATABLE :: cov_xy_p(:)     ! covariances between state and single observation
    REAL, ALLOCATABLE :: cov_hxy_p(:)    ! covariances between full observed state and single observation


! **********************
! *** INITIALIZATION ***
! **********************

    CALL PDAF_timeit(51, 'new')

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_ensrf_analysis -- START'

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, i7, 3x, a)') &
            'PDAF ', step, 'Local least squares EAKF with serial observation processing'
    END IF

    ! init numbers
    invdim_ensm1 = 1.0 / (REAL(dim_ens - 1))

    ! Allocate arrays
    ALLOCATE(HXpert_i(dim_ens))
    ALLOCATE(HXinc_i(dim_ens))
    ALLOCATE(cov_xy_p(dim_p))
    ALLOCATE(cov_hxy_p(dim_p))
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

       ! Get mean of observed ensemble for current observation
       HXbar_i = HXbar_p(iobs)

       ! Get observed ensemble perturbations for current observation
       HXpert_i(:) = HX_p(iobs,:) - HXbar_i

       ! Compute variance of observed ensemble for iobs
       var_hx = 0.0
       DO member = 1, dim_ens
          var_hx = var_hx + HXpert_i(member)**2
       END DO
       var_hx = var_hx * invdim_ensm1

       ! Compute ration of variance
       var_ratio = var_obs_p(iobs) / (var_hx + var_obs_p(iobs))

       CALL PDAF_timeit(30, 'old')
       CALL PDAF_timeit(31, 'new')

       ! Compute covariances between state ensemble and observation
       cov_xy_p = 0.0
       DO member = 1, dim_ens
          cov_xy_p(:) = cov_xy_p(:) + ens_p(:, member) * HXpert_i(member)
       END DO
       cov_xy_p = cov_xy_p * invdim_ensm1

       ! Compute covariances between observed state ensemble and observation
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

       CALL PDAF_timeit(45, 'new')
       CALL U_localize_covar_serial(iobs, dim_p, dim_obs_p, cov_xy_p, cov_hxy_p)
       CALL PDAF_timeit(45, 'old')

       CALL PDAF_timeit(10, 'old')


       ! ****************************************************
       ! *** Compute the observation increment            ***
       ! ****************************************************

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
       ! *** Update state ensemble                        ***
       ! ****************************************************

       CALL PDAF_timeit(14, 'new')

       ! store ratio of covariances in cov_xy_p
       cov_xy_p(:) = cov_xy_p(:) / var_hx

       ! Update ensemble members
       DO member = 1, dim_ens
          ens_p(:, member) = ens_p(:, member) + cov_xy_p(:) * HXinc_i(member) 
       END DO

       CALL PDAF_timeit(14, 'old')


       ! ****************************************************
       ! *** Update observed ensemble and its mean        ***
       ! ****************************************************

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

    IF (allocflag == 0) allocflag = 1

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_ensrf_analysis -- END'

  END SUBROUTINE PDAF_ensrf_ana_2step

END MODULE PDAF_ensrf_analysis
