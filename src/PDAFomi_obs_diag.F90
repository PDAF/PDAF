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
       have_obsens_diag, rmsd, dim_obs_diag_p, obsstats

  IMPLICIT NONE

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

    ! Set value for observation diagnostics type
    obs_diag = diag

  END SUBROUTINE PDAFomi_set_obs_diag



!-------------------------------------------------------------------------------
!!> Return number of observation type in diagnostics
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_diag_nobs(nobs)

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

  END SUBROUTINE PDAFomi_diag_nobs



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

    ! Pre-initialize observation dimension
    dim_obs_diag = 0

    ! Check whether OMI is used and observation diagnostics are used
    IF (n_obstypes >= 0 .AND. n_obstypes >= id_obs &
         .AND. (have_obsmean_diag > 0 .OR. have_obsens_diag > 0)) THEN

       IF (obs_f_all(id_obs)%ptr%dim_obs_p >0) THEN

          ! Set observation dimension and number of coordinates
          dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p
          ncoord = obs_f_all(id_obs)%ptr%ncoord

          ! Set pointers
          obs_p_ptr => obs_f_all(id_obs)%ptr%obs_diag_p
          ocoord_p_ptr => obs_f_all(id_obs)%ptr%ocoord_diag_p

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

             ! Set pointer
             HXmean_p_ptr => obs_f_all(id_obs)%ptr%HXmean_diag_p

          ELSE IF (have_obsens_diag>0) THEN

             ! *** Case in which only the observed ensemble is initialized ***
             ! *** here we need to compue the observed ensemble mean       ***

             ! Set observation dimension
             dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p

             ! Compute mean
             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, &
                  obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)
 
             ! Set pointer
             HXmean_p_ptr => obs_f_all(id_obs)%ptr%HXmean_diag_p

          END IF

       END IF

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

          ! Set pointer
          IF (have_obsens_diag>0) THEN

             ! *** Case in which only the observed ensemble is initialized ***
             ! *** here we need to compue the observed ensemble mean       ***

             ! Set observation dimension
             dim_obs_diag = obs_f_all(id_obs)%ptr%dim_obs_p
 
             ! Set pointer
             HX_p_ptr => obs_f_all(id_obs)%ptr%HX_diag_p

          END IF

       END IF

    END IF

  END SUBROUTINE PDAFomi_diag_get_HX



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
  SUBROUTINE PDAFomi_diag_obs_rmsd(nobs, rmsd_pointer)

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

          IF (have_obsens_diag>0) THEN

             ! *** When only the observed ensemble is initialized ***
             ! *** we need to compue the observed ensemble mean   ***

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

    END IF haveobs

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
  SUBROUTINE PDAFomi_diag_obs_stats(nobs, obsstats_ptr)

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
    INTEGER, INTENT(inout) :: nobs                     !< Number of observation types
    REAL, POINTER, INTENT(inout) :: obsstats_ptr(:,:)  !< Array of observation statistics

! *** Local variables ***
    INTEGER :: i, id_obs             ! Counters
    INTEGER :: nstats                ! Number of observation statistics
    INTEGER :: status                ! Status flag
    INTEGER :: dim_g                 ! Global number of observations of one obs. type
    REAL :: rmsd_p                   ! PE-local RMSd
    REAL :: stat_g                   ! Global statistic
    INTEGER :: MPIerr                ! MPI status flag


! ***********************
! *** Initialization  ***
! ***********************

    ! Pre-initialize nobs
    nobs = 0

    ! Number of observation statistics
    nstats = 1


    ! Allocate RMSD vector
    habeobs: IF (n_obstypes > 0  .AND. (have_obsmean_diag>0 .OR. have_obsens_diag>0)) THEN

       ! Set number of obstypes
       nobs = n_obstypes


! ***********************************************
! *** Initialize pointer to statistics array  ***
! ***********************************************

       IF (ALLOCATED(obsstats)) DEALLOCATE(obsstats)
       ALLOCATE(obsstats(n_obstypes, nstats))
       obsstats = 0.0

       ! Set pointer
       obsstats_ptr => obsstats


! ***********************************************
! *** Compute RMSD for each observation type  ***
! ***********************************************

       DO id_obs = 1, n_obstypes

          ! Get global state dimension
          CALL MPI_Allreduce(obs_f_all(id_obs)%ptr%dim_obs_p, dim_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

          IF (have_obsens_diag>0) THEN

             ! *** When only the observed ensemble is initialized ***
             ! *** we need to compue the observed ensemble mean   ***

             CALL PDAF_diag_ensmean(obs_f_all(id_obs)%ptr%dim_obs_p, dim_ens, obs_f_all(id_obs)%ptr%HXmean_diag_p, &
                  obs_f_all(id_obs)%ptr%HX_diag_p, status)
          END IF

          IF (have_obsmean_diag>0 .OR. have_obsens_diag>0) THEN

          ! Check that either the observed mean of ensmeble was initialized
             rmsd_p = 0.0
             DO i = 1, obs_f_all(id_obs)%ptr%dim_obs_p
                rmsd_p = rmsd_p + (obs_f_all(id_obs)%ptr%obs_diag_p(i) - obs_f_all(id_obs)%ptr%HXmean_diag_p(i))**2
             END DO

             ! *** Complete computation of global RMSD ***

             ! normalize PE-local stddev
             rmsd_p = rmsd_p / REAL(dim_g)

             ! Get global stddev
             CALL MPI_Allreduce(rmsd_p, stat_g, 1, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

             ! Complete computation of global rmsd
             obsstats(id_obs, 1) = SQRT(stat_g)

          END IF

       END DO

    END IF habeobs

  END SUBROUTINE PDAFomi_diag_obs_stats

END MODULE PDAFomi_obs_diag
