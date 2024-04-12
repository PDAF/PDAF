! Copyright (c) 2004-2024 Lars Nerger
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
!$Id: PDAFomi_obs_f.F90 1147 2023-03-12 16:14:34Z lnerger $

!> Utility routines for analysis step
!!
MODULE PDAF_analysis_utils

  ! Array for observation statistics in local analysis
  INTEGER :: obsstats(4)           ! PE-local statistics
  ! obsstats(1): Local domains with observations
  ! obsstats(2): Local domains without observations
  ! obsstats(3): Sum of all available observations for all domains
  ! obsstats(4): Maximum number of observations over all domains


CONTAINS

!-------------------------------------------------------------------------------

!> Print statistics on local analysis domains
!!
  SUBROUTINE PDAF_print_domain_stats(n_domains_p)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, COMM_filter, MPIerr

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: n_domains_p   ! Number of PE-local analysis domains

! *** Local Variables ***
    INTEGER :: n_domains_stats(4)        ! Gobal statistics for number of analysis domains


    IF (npes_filter>1) THEN
       CALL MPI_Reduce(n_domains_p, n_domains_stats(1), 1, MPI_INTEGER, MPI_MIN, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(n_domains_p, n_domains_stats(2), 1, MPI_INTEGER, MPI_MAX, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(n_domains_p, n_domains_stats(3), 1, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
       IF (mype == 0) THEN
          WRITE (*, '(a, 5x, a, i9, 1x, i9, 1x, f11.1)') &
               'PDAF', '--- local analysis domains (min/max/avg):', n_domains_stats(1:2), &
               REAL(n_domains_stats(3)) / REAL(npes_filter)
       END IF
    ELSE
       ! This is a work around for working with nullmpi.F90
       IF (mype == 0) THEN
          WRITE (*, '(a, 5x, a, i9)') &
               'PDAF', '--- local analysis domains:', n_domains_p
       END IF
    END IF

  END SUBROUTINE PDAF_print_domain_stats


!-------------------------------------------------------------------------------

!> Initialize local observation statistics
!!
  SUBROUTINE PDAF_init_local_obsstats()

    IMPLICIT NONE

! *** Set obsstats to zero ***

    obsstats = 0

  END SUBROUTINE PDAF_init_local_obsstats


!> Update local observation statistics
!!
  SUBROUTINE PDAF_incr_local_obsstats(dim_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l   ! Number of locally assimilated observations

! *** Update observation statistics ***

!$OMP CRITICAL
    IF (dim_obs_l > obsstats(4)) obsstats(4) = dim_obs_l
    IF (dim_obs_l > 0) THEN
       obsstats(3) = obsstats(3) + dim_obs_l
       obsstats(1) = obsstats(1) + 1
    ELSE
       obsstats(2) = obsstats(2) + 1
    END IF
!$OMP END CRITICAL

  END SUBROUTINE PDAF_incr_local_obsstats



!-------------------------------------------------------------------------------

!> Print observation statistics
!!
  SUBROUTINE PDAF_print_local_obsstats(screen, n_domains_with_obs)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, COMM_filter, MPIerr
    USE PDAFomi, &
         ONLY: omi_n_obstypes => n_obstypes, PDAFomi_obsstats_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: screen      ! Verbosity flag
    INTEGER, OPTIONAL, INTENT(out) :: n_domains_with_obs

! *** Local variables ***
    INTEGER :: obsstats_g(4)           ! Global statistics

    ! *** Print statistics for local analysis to the screen ***
    IF ( npes_filter>1) THEN
       CALL MPI_Reduce(obsstats, obsstats_g, 3, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(obsstats(4), obsstats_g(4), 1, MPI_INTEGER, MPI_MAX, &
            0, COMM_filter, MPIerr)
    ELSE
       ! This is a work around for working with nullmpi.F90
       obsstats_g = obsstats
    END IF

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- Global statistics for local analysis:'
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Local domains with observations:       ', obsstats_g(1)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Local domains without observations:    ', obsstats_g(2)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Maximum local observation dimension:   ', obsstats_g(4)
       IF (obsstats_g(2) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Total avg. local observation dimension:', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1) + obsstats_g(2))
       END IF
       IF (obsstats_g(2) > 0 .AND. obsstats_g(1) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Avg. for domains with observations:    ', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1))
       END IF
    END IF

    IF (omi_n_obstypes > 0) CALL PDAFomi_obsstats_l(screen)

    IF (PRESENT(n_domains_with_obs)) n_domains_with_obs = obsstats_g(1)

  END SUBROUTINE PDAF_print_local_obsstats



!-------------------------------------------------------------------------------

!> Compute innovation and call PDAFomi_omit obs
!!
!! This routine is used by some of the global filters
!! (EnKF, LEnKF, PF, NETF) with OMI to omit observations
!! if the innovation is too large. These filters do not
!! compute the innovation for the ensemble mean, but this
!! is needed to omit observations. Thus, this step is 
!! done here and then the omission routine of OMI is called.
!! 
  SUBROUTINE PDAF_omit_obs_omi(dim_p, dim_obs_p, dim_ens, state_p, ens_p, &
       obs_p, U_init_obs, U_obs_op, compute_mean, screen)
  
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_mod_filter, &
         ONLY: obs_member, debug

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p          !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_obs_p      !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens        !< Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p) !< on exit: PE-local forecast mean state
    REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
    REAL, INTENT(inout) :: obs_p(dim_obs_p)   !< PE-local observation vector
    INTEGER, INTENT(in) :: compute_mean   !< (1) compute mean; (0) state_p holds mean
    INTEGER, INTENT(in) :: screen         !< Verbosity flag

! External subroutines 
! (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_obs, &             ! Initialize observation vector
         U_obs_op                         ! Observation operator


! *** local variables ***
    INTEGER :: row, member                ! Counters
    REAL    :: invdimens                  ! Inverse global ensemble size
    REAL, ALLOCATABLE :: resid_p(:)       ! PE-local observation residual


! ******************************************
! *** Compute residual for ensmeble mean ***
! ******************************************

    IF (compute_mean==1) THEN
       CALL PDAF_timeit(51, 'new')

    ! *** Compute mean forecast state
       state_p = 0.0
       invdimens = 1.0 / REAL(dim_ens)
       DO member = 1, dim_ens
          DO row = 1, dim_p
             state_p(row) = state_p(row) + invdimens * ens_p(row, member)
          END DO
       END DO

       CALL PDAF_timeit(51, 'old')
    END IF

    ALLOCATE(resid_p(dim_obs_p))

    ! Apply observation operator to ensemble mean state
    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_omit_obs_omi -- call obs_op'

    CALL PDAF_timeit(44, 'new')
    obs_member = 0
    CALL U_obs_op(step, dim_p, dim_obs_p, state_p, resid_p)
    CALL PDAF_timeit(44, 'old')

    ! Initialize vector of observations
    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_omit_obs_omi -- call init_obs'

    CALL PDAF_timeit(50, 'new')
    CALL U_init_obs(step, dim_obs_p, obs_p)
    CALL PDAF_timeit(50, 'old')

    ! Compute residual
    CALL PDAF_timeit(51, 'new')
    resid_p = obs_p - resid_p


! **************************************************
! *** Omit observations with too high innovation ***
! **************************************************

    CALL PDAFomi_omit_by_inno_cb(dim_obs_p, resid_p, obs_p)

    CALL PDAF_timeit(51, 'old')

    DEALLOCATE(resid_p)

  END SUBROUTINE PDAF_omit_obs_omi

END MODULE PDAF_analysis_utils
