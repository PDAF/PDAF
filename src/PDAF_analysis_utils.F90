MODULE PDAF_analysis_utils

  ! Array for observation statistics in local analysis
  INTEGER :: obsstats(4)           ! PE-local statistics
  ! obsstats(1): Local domains with observations
  ! obsstats(2): Local domains without observations
  ! obsstats(3): Sum of all available observations for all domains
  ! obsstats(4): Maximum number of observations over all domains

CONTAINS

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
          WRITE (*, '(a, 5x, a, i7, 1x, i7, 1x, f9.1)') &
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


!> Print observation statistics
!!
  SUBROUTINE PDAF_print_local_obsstats(screen)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, COMM_filter, MPIerr
    USE PDAFomi, &
         ONLY: omi_n_obstypes => n_obstypes, PDAFomi_obsstats

    IMPLICIT NONE

    INTEGER, INTENT(in) :: screen      ! Verbosity flag

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
       IF (obsstats_g(2) > 0 .AND. obsstats_g(1) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Total avg. local observation dimension:', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1) + obsstats_g(2))
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Avg. for domains with observations:    ', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1))
       END IF
    END IF

    if (omi_n_obstypes > 0) CALL PDAFomi_obsstats(obsstats_g, screen)

  END SUBROUTINE PDAF_print_local_obsstats

END MODULE PDAF_analysis_utils
