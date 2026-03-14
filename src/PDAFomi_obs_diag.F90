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
!> PDAF-OMI observation diagnostics
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFomi_obs_diag

  USE PDAFomi_obs_f, &
       ONLY: obs_diag, obs_f_all, n_obstypes, have_obsmean_diag, &
       have_obsens_diag, rmsd, dim_obs_diag_p, obsstats, crps_allobs

  IMPLICIT NONE

  PRIVATE :: obs_diag, obs_f_all, n_obstypes, have_obsmean_diag, &
       have_obsens_diag, rmsd, dim_obs_diag_p, obsstats, crps_allobs

CONTAINS

!-------------------------------------------------------------------------------
!!> Set observation diagnostics flag
!!
!! This routine sets the flag for observation diagnostics 
!! in PDAF-OMI. One should set this flag before calling
!! PDAFomi_gather_obs.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_obs_diag(diag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: diag          !< Value for observation diagnostics mode
                                         !< >0 activates observation diagnostics

    ! Set value for observation diagnostics type
    obs_diag = diag

  END SUBROUTINE PDAFomi_set_obs_diag



!-------------------------------------------------------------------------------
!!> Return number of observation types in diagnostics
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_nobstypes(nobs)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs                   !< Number of observation types


! ***************************************
! *** Set number of observation types ***
! ***************************************

    ! Pre-initialize nobs
    nobs = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes > 0  .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN
       nobs = n_obstypes
    END IF

  END SUBROUTINE PDAFomi_diag_nobstypes



!-------------------------------------------------------------------------------
!!> Return dimension of observation vector for single obs. type
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_dimobs(dim_obs_ptr)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, POINTER, INTENT(inout) :: dim_obs_ptr(:)   !< Pointer to observation dimensions

! *** Local variables ***
    INTEGER :: iobs               ! Counters


! **************************************
! *** Set dimensions of observations ***
! **************************************

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0  .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       IF (ALLOCATED(dim_obs_diag_p)) DEALLOCATE(dim_obs_diag_p)
       ALLOCATE(dim_obs_diag_p(n_obstypes))

       ! Pre-initialize nobs
       dim_obs_diag_p = 0

       ! Initialize vector of observation dimensions
       DO iobs = 1, n_obstypes
          dim_obs_diag_p(iobs) = obs_f_all(iobs)%ptr%dim_obs_p
       END DO

       ! Set pointer
       dim_obs_ptr => dim_obs_diag_p

    ELSE

       IF (ALLOCATED(dim_obs_diag_p)) DEALLOCATE(dim_obs_diag_p)
       ALLOCATE(dim_obs_diag_p(1))
       dim_obs_diag_p = 0

       ! Set pointer
       dim_obs_ptr => dim_obs_diag_p

    END IF

  END SUBROUTINE PDAFomi_diag_dimobs


!-------------------------------------------------------------------------------
!!> Return pointers to observation vector and related coordinates
!!
!! This routine returns the pointer to the PE-local observation
!! and coordinate array.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_get_obs(id_obs, dim_obs_diag, ncoord, obs_p_ptr, ocoord_p_ptr)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: id_obs                    !< Index of observation type to return
    INTEGER, INTENT(out) :: dim_obs_diag             !< Observation dimension
    INTEGER, INTENT(out) :: ncoord                   !< Number of observation dimensions
    REAL, POINTER, INTENT(out) :: obs_p_ptr(:)       !< Pointer to observation vector
    REAL, POINTER, INTENT(out) :: ocoord_p_ptr(:,:)  !< Pointer to Coordinate array


! ***************************************************
! *** Set pointer to observations and coordinates ***
! ***************************************************

    ! Pre-initialize dimensions
    dim_obs_diag = 0
    ncoord = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0 .AND. n_obstypes >= id_obs &
         .AND. (have_obsmean_diag > 0 .OR. have_obsens_diag > 0)) THEN

       ! Set pointers
       obs_p_ptr => obs_f_all(id_obs)%ptr%obs_diag_p
       ocoord_p_ptr => obs_f_all(id_obs)%ptr%ocoord_diag_p

       IF (obs_f_all(id_obs)%ptr%dim_obs_p >0) THEN

          ! Set observation dimension and number of coordinates
          dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p
          ncoord = obs_f_all(id_obs)%ptr%ncoord

       END IF

    END IF

  END SUBROUTINE PDAFomi_diag_get_obs


!-------------------------------------------------------------------------------
!!> Return pointer to mean of observed ensemble
!!
!! This routine returns the pointer to the PE-local
!! observed ensemble mean.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_get_HXmean(id_obs, dim_obs_diag, HXmean_p_ptr)

    USE PDAF_diag, &
         ONLY: PDAF_diag_ensmean
    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: id_obs                    !< Index of observation type to return
    INTEGER, INTENT(out) :: dim_obs_diag             !< Observation dimension
    REAL, POINTER, INTENT(out) :: HXmean_p_ptr(:)    !< Pointer to observed ensemble mean

! *** Local variables ***
    INTEGER :: status                ! Status flag


! *********************************************
! *** Set pointer to observed ensemble mean ***
! *********************************************

    ! Pre-initialize observation dimension
    dim_obs_diag = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0 .AND. n_obstypes >= id_obs &
         .AND. (have_obsmean_diag > 0 .OR. have_obsens_diag > 0)) THEN

       IF (obs_f_all(id_obs)%ptr%dim_obs_p >0) THEN

          ! Set pointers
          IF (have_obsmean_diag>0) THEN

             ! *** Case in which the observed ensemble mean is directly initialized ***

             ! Set observation dimension
             dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p

          ELSE IF (have_obsens_diag>0) THEN

             ! *** Case in which only the observed ensemble is initialized ***
             ! *** here we need to compute the observed ensemble mean       ***

             ! Set observation dimension
             dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p

             ! Compute mean
             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, &
                  obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)

          END IF

       END IF

       ! Set pointer
       HXmean_p_ptr => obs_f_all(id_obs)%ptr%HXmean_diag_p

    END IF

  END SUBROUTINE PDAFomi_diag_get_HXmean


!-------------------------------------------------------------------------------
!!> Return pointer to observed ensemble
!!
!! This routine returns the pointer to the PE-local
!! observed ensemble.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_get_HX(id_obs, dim_obs_diag, HX_p_ptr)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: id_obs                    !< Index of observation type to return
    INTEGER, INTENT(out) :: dim_obs_diag             !< Observation dimension
    REAL, POINTER, INTENT(out) :: HX_p_ptr(:,:)      !< Pointer to observed ensemble mean


! ****************************************
! *** Set pointer to observed ensemble ***
! ****************************************

    ! Pre-initialize observation dimension
    dim_obs_diag = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0 .AND. n_obstypes >= id_obs &
         .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       IF (obs_f_all(id_obs)%ptr%dim_obs_p >0) THEN
          IF (have_obsens_diag>0) THEN

             ! Set observation dimension
             dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p

          END IF
       END IF
 
       ! Set pointer
       HX_p_ptr => obs_f_all(id_obs)%ptr%HX_diag_p

    END IF

  END SUBROUTINE PDAFomi_diag_get_HX


!-------------------------------------------------------------------------------
!!> Return pointer to observed ensemble
!!
!! This routine returns the pointer to the PE-local
!! observed ensemble.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_get_ivar(id_obs, dim_obs_diag, ivar_ptr)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: id_obs                    !< Index of observation type to return
    INTEGER, INTENT(out) :: dim_obs_diag             !< Observation dimension
    REAL, POINTER, INTENT(out) :: ivar_ptr(:)        !< Pointer to inverse observation error variances


! *******************************************************************
! *** Set pointer to vectorof inverse observation error variances ***
! *******************************************************************

    ! Pre-initialize observation dimension
    dim_obs_diag = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0 .AND. n_obstypes >= id_obs &
         .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       IF (obs_f_all(id_obs)%ptr%dim_obs_p >0) THEN
          ! Set observation dimension
          dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p
       END IF
 
       ! Set pointer
       ivar_ptr => obs_f_all(id_obs)%ptr%ivar_obs_diag_p

    END IF

  END SUBROUTINE PDAFomi_diag_get_ivar



!-------------------------------------------------------------------------------
!!> Compute RMS deviation beetween observation and observed ensemble mean
!!
!! This routine is an alias for PDAFomi_diag_rmsd
!!
!! __Revision history:__
!! * 2026-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_obs_rmsd(nobs, rmsd_pointer, verbose)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs                   !< Number of observation types
    REAL, POINTER, INTENT(inout) :: rmsd_pointer(:)  !< Vector of RMSD values
    INTEGER, INTENT(in) :: verbose                   !< Verbosity flag

    CALL PDAFomi_diag_rmsd(nobs, rmsd_pointer, verbose)
   
  END SUBROUTINE PDAFomi_diag_obs_rmsd



!-------------------------------------------------------------------------------
!!> Compute RMS deviation beetween observation and observed ensemble mean
!!
!! This routine computes the RMSD between the observation and the
!! observed ensemble mean.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_rmsd(nobs, rmsd_pointer, verbose)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE MPI
    USE PDAF_mod_parallel, &
         ONLY: COMM_filter
    USE PDAF_diag, &
         ONLY: PDAF_diag_ensmean
    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs                   !< Number of observation types
    REAL, POINTER, INTENT(inout) :: rmsd_pointer(:)  !< Vector of RMSD values
    INTEGER, INTENT(in) :: verbose                   !< Verbosity flag

! *** Local variables ***
    INTEGER :: i, id_obs             ! Counters
    INTEGER :: status                ! Status flag
    INTEGER :: dim_g                 ! Global number of observations of one obs. type
    REAL :: rmsd_p                   ! PE-local RMSd
    INTEGER :: MPIerr                ! MPI status flag


! ***********************
! *** Initialization  ***
! ***********************

    ! Pre-initialize nobs
    nobs = 0

    ! Allocate RMSD vector
    haveobs: IF (n_obstypes > 0  .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       ! Set number of obstypes
       nobs = n_obstypes


! ********************************
! *** Initialize RMSD pointer  ***
! ********************************

       IF (ALLOCATED(rmsd)) DEALLOCATE(rmsd)
       ALLOCATE(rmsd(n_obstypes))
       rmsd = 0.0

       ! Set pointer
       rmsd_pointer => rmsd


! ***********************************************
! *** Compute RMSD for each observation type  ***
! ***********************************************

       DO id_obs = 1, n_obstypes

          ! Get global state dimension
          CALL MPI_Allreduce(obs_f_all(id_obs)%ptr%dim_obs_p, dim_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

          IF (have_obsmean_diag==0 .AND. have_obsens_diag>0) THEN

             ! *** When only the observed ensemble is initialized ***
             ! *** we need to compute the observed ensemble mean   ***

             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)
          END IF

          ! Check that either the observed mean of ensmeble was initialized
          IF (have_obsmean_diag>0 .OR. have_obsens_diag>0) THEN

             rmsd_p = 0.0
             DO i = 1, obs_f_all(id_obs)%ptr%dim_obs_p
                rmsd_p = rmsd_p + (obs_f_all(id_obs)%ptr%obs_diag_p(i) - obs_f_all(id_obs)%ptr%HXmean_diag_p(i))**2
             END DO

             ! *** Complete computation of global RMSD ***

             ! normalize PE-local stddev
             rmsd_p = rmsd_p / REAL(dim_g)

             ! Get global stddev
             CALL MPI_Allreduce(rmsd_p, rmsd(id_obs), 1, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

             ! Complete computation of global rmsd
             rmsd(id_obs) = SQRT(rmsd(id_obs))

          END IF

       END DO

       IF (verbose>0) THEN
          WRITE(*,'(a, 5x, a)') 'PDAFomi', 'RMS deviations between observations and observed ensemble'
          WRITE(*,'(a, 3x, a, 4x, a)') 'PDAFomi', 'obs-ID', '  RMSD   '
          DO id_obs = 1, n_obstypes
             WRITE (*, '(a, 4x, i3, 3x, es12.3)') 'PDAFomi', id_obs, rmsd(id_obs)
          END DO
       END IF

    END IF haveobs

  END SUBROUTINE PDAFomi_diag_rmsd



!-------------------------------------------------------------------------------
!!> Compute statistics for difference beetween observation and observed ensemble mean
!!
!! This routine computes different statistics comparing
!! the observation vector and the observed ensemble mean.
!! The array obstats is organized so that the first
!! index is the type of statistics, while the second
!! index specified the ID of the observation type.
!!
!! The statistics are in the array obsstats as follows
!! * 1: correlation
!! * 2: centered RMS deviation
!! * 3: bias observation - observed ensemble mean
!! * 4: mean absolute deviation observation - observed ensemble mean
!! * 5: standard deviation of observations
!! * 6: standard deviation of observed ensemble mean
!! The statistics 1, 5, ans 6 are the usual values used
!! in Taylor diagrams.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_stats(nobs, obsstats_ptr, verbose)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE MPI
    USE PDAF_mod_parallel, &
         ONLY: COMM_filter
    USE PDAF_diag, &
         ONLY: PDAF_diag_ensmean, PDAF_diag_stats
    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs                     !< Number of observation types
    REAL, POINTER, INTENT(inout) :: obsstats_ptr(:,:)  !< Array of observation statistics
    INTEGER, INTENT(in) :: verbose                     !< Verbosity flag

! *** Local variables ***
    INTEGER :: i, id_obs             ! Counters
    INTEGER :: nstats                ! Number of observation statistics
    INTEGER :: status                ! Status flag
    INTEGER :: dim_g                 ! Global number of observations of one obs. type
    INTEGER :: MPIerr                ! MPI status flag
    REAL :: stats_p(6)               ! PE-local statistics array
    REAL :: stats_g(6)               ! Global statistics array
    REAL :: mean_obs                 ! mean of observation vector
    REAL :: mean_HXmean              ! mean of observed ensemble mean
    REAL :: means_p(2), means_g(2)   ! mean observations and obs. ensemble mean
    REAL :: mad_p                    ! PE-local mean absolute deviation
    REAL :: crmsd_p                  ! PE-local centered RMS difference
    REAL :: corr_p                   ! PE-local centered RMS difference
    REAL :: var_o_p                  ! PE-local centered RMS difference
    REAL :: var_Hx_p                 ! PE-local centered RMS difference


! ***********************
! *** Initialization  ***
! ***********************

    ! Pre-initialize nobs
    nobs = 0

    ! Number of observation statistics
    nstats = 6


    ! Allocate RMSD vector
    haveobs: IF (n_obstypes > 0  .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       ! Set number of obstypes
       nobs = n_obstypes


! ***********************************************
! *** Initialize pointer to statistics array  ***
! ***********************************************

       IF (ALLOCATED(obsstats)) DEALLOCATE(obsstats)
       ALLOCATE(obsstats(nstats, n_obstypes))
       obsstats = 0.0

       ! Set pointer
       obsstats_ptr => obsstats


! *****************************************************
! *** Compute statistics for each observation type  ***
! *****************************************************

       DO id_obs = 1, n_obstypes

          IF (have_obsens_diag>0) THEN
             ! *** When only the observed ensemble is initialized ***
             ! *** we need to compute the observed ensemble mean   ***

             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)
          END IF

          ! Call statistics routine of PDAF_diag
          CALL PDAF_diag_stats(obs_f_all(id_obs)%ptr%dim_obs_p, obs_f_all(id_obs)%ptr%obs_diag_p, &
               obs_f_all(id_obs)%ptr%HXmean_diag_p, obsstats(:,id_obs), 0)

       END DO

       IF (verbose>0) THEN
          WRITE(*,'(a, 5x, a)') 'PDAFomi', 'Statistics on deviations: observation y - observed ensemble mean Hx'
          WRITE(*,'(a, 3x, a, 1x, 5(a, 3x), a)') 'PDAFomi', 'obs-ID', '  corr   ', 'cRMSD  ', &
               'bias y-Hx', '   MAD   ', 'STDDEV(y)', 'STDDEV(Hx)'
          DO id_obs = 1, n_obstypes
             WRITE (*, '(a, 4x, i3, 3x, f7.3, 5es12.3)') &
                  'PDAFomi', id_obs, obsstats(1,id_obs), obsstats(2:6,id_obs)
          END DO
       END IF

    END IF haveobs

  END SUBROUTINE PDAFomi_diag_stats


!-------------------------------------------------------------------------------
!!> Compute CRPS between observation and observed ensemble
!!
!! This routine computes the CRPS between the observation and the
!! observed ensemble. It return a poiter to an array with values
!! for each observation type including the CRPS, reliability,
!! potential CRPS and uncertainty.
!!
!! __Revision history:__
!! * 2026-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_crps(nobs, crps_pointer, perturb, verbose)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_diag, &
         ONLY: PDAF_diag_crps
    USE PDAF_utils, &
         ONLY: PDAF_generate_rndvec
    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs                    !< Number of observation types
    REAL, POINTER, INTENT(inout) :: crps_pointer(:,:) !< Array of CRPS values
    INTEGER, INTENT(in) :: perturb                    !< 1 to add perturbations to observations
    INTEGER, INTENT(in) :: verbose                    !< Verbosity flag

! *** Local variables ***
    INTEGER :: id_obs                ! Counters
    INTEGER :: status                ! Status flag
    REAL, ALLOCATABLE :: rndvec(:)   ! Gaussian perturbations
    REAL, ALLOCATABLE :: obs_pert(:) ! Perturbed observations vector
    INTEGER, SAVE :: iseed(4)        ! Seed for random number generation
    LOGICAL, SAVE :: first=.true.    ! Flag for first call
    REAL :: CRPS                     ! CRPS
    REAL :: reli                     ! Reliability
    REAL :: pot_CRPS                 ! Potential CRPS
    REAL :: uncert                   ! Uncertainty


! ***********************
! *** Initialization  ***
! ***********************

    ! Pre-initialize nobs
    nobs = 0

    haveobs: IF (n_obstypes > 0  .AND. (have_obsens_diag>0)) THEN

       ! Set number of obstypes
       nobs = n_obstypes


! ********************************
! *** Initialize CRPS pointer  ***
! ********************************

       IF (ALLOCATED(crps_allobs)) DEALLOCATE(crps_allobs)
       ALLOCATE(crps_allobs(4, n_obstypes))
       crps_allobs = 0.0

       ! Set pointer
       crps_pointer => crps_allobs


! ***********************************************
! *** Compute CRPS for each observation type  ***
! ***********************************************

       pert: IF (perturb == 0) THEN

          ! Compute CRPS without perturbing the observations

          DO id_obs = 1, n_obstypes

             CALL PDAF_diag_crps(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, 0, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, obs_f_all(id_obs)%ptr%obs_diag_p, &
                  CRPS, reli, pot_CRPS, uncert, status)

             crps_allobs(1, id_obs) = CRPS
             crps_allobs(2, id_obs) = reli
             crps_allobs(3, id_obs) = pot_CRPS
             crps_allobs(4, id_obs) = uncert

          END DO

       ELSE pert

          ! Compute CRPS with perturbed observations

          DO id_obs = 1, n_obstypes

             ALLOCATE(obs_pert(obs_f_all(id_obs)%ptr%dim_obs_p))
             ALLOCATE(rndvec(obs_f_all(id_obs)%ptr%dim_obs_p))

             IF (first) THEN
                iseed(1) = 1053
                iseed(2) = 153
                iseed(3) = 53
                iseed(4) = 2053
                first = .false.
             END IF

             ! Generate normally distribute random vector 
             CALL PDAF_generate_rndvec(obs_f_all(id_obs)%ptr%dim_obs_p, rndvec, 1.0, 1, iseed)

             ! Create perturbed observation vector
             obs_pert = rndvec/SQRT(obs_f_all(id_obs)%ptr%ivar_obs_diag_p(:)) + obs_f_all(id_obs)%ptr%obs_diag_p(:)

             ! Compute CPRS and its decomposition
             CALL PDAF_diag_crps(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, 0, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, obs_pert, &
                  CRPS, reli, pot_CRPS, uncert, status)

             crps_allobs(1, id_obs) = CRPS
             crps_allobs(2, id_obs) = reli
             crps_allobs(3, id_obs) = pot_CRPS
             crps_allobs(4, id_obs) = uncert

             DEALLOCATE(obs_pert, rndvec)

          END DO

       END IF pert

       IF (verbose>0) THEN
          WRITE(*,'(a, 5x, a)') 'PDAFomi', 'Decomposed CRPS for each observation type'
          WRITE(*,'(a, 3x, a, 1x, a,2x, 4(a, 1x))') 'PDAFomi', 'obs-ID', &
               '   CRPS    ', 'reliability', ' pot. CRPS ', 'uncertainty'
          DO id_obs = 1, n_obstypes
             WRITE (*, '(a, 4x, i3, 1x, 4es12.3)') &
                  'PDAFomi', id_obs, crps_allobs(1:4, id_obs)
          END DO
       END IF

    END IF haveobs

  END SUBROUTINE PDAFomi_diag_crps


!-------------------------------------------------------------------------------
!!> Set omitted observation by high observation error for diagnistics only
!!
!! This routine checks the difference between observations and observed
!! ensemble mean and sets the inverse observation error to a very 
!! small value is the difference is too large.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_omit_by_inno()

    USE PDAF_diag, &
         ONLY: PDAF_diag_ensmean
    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Local variables ***
    INTEGER :: i, id_obs, cnt       ! Counters
    INTEGER :: status               ! Status flag
    REAL :: inno2                   ! Squared innovation
    REAL :: limit2                  ! Squared limit


    allobs: DO id_obs = 1, n_obstypes

       omit: IF (obs_f_all(id_obs)%ptr%inno_omit > 0.0) THEN

          IF (have_obsmean_diag==0 .AND. have_obsens_diag>0) THEN

             ! *** When only the observed ensemble is initialized ***
             ! *** we need to compute the observed ensemble mean   ***

             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)
          END IF
          
          haveobs: IF (have_obsmean_diag>0 .OR. have_obsens_diag>0) THEN

             ! Squared limit factor
             limit2 = obs_f_all(id_obs)%ptr%inno_omit * obs_f_all(id_obs)%ptr%inno_omit

             ! Check for observations to be excluded
             cnt = 0
             DO i = 1, obs_f_all(id_obs)%ptr%dim_obs_p

                ! Squared innovation
                inno2 = (obs_f_all(id_obs)%ptr%obs_diag_p(i) - obs_f_all(id_obs)%ptr%HXmean_diag_p(i))**2

                IF (inno2 > limit2 / obs_f_all(id_obs)%ptr%ivar_obs_diag_p(i)) THEN

                   ! Exclude observation by increased its observation error
                   obs_f_all(id_obs)%ptr%ivar_obs_diag_p(i) = obs_f_all(id_obs)%ptr%inno_omit_ivar

                   ! Count excluded obs
                   cnt = cnt + 1
                END IF
             ENDDO

          END IF haveobs

       END IF omit

    END DO allobs

  END SUBROUTINE PDAFomi_diag_omit_by_inno

END MODULE PDAFomi_obs_diag
