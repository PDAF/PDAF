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
!> Module for ensmeble communication in COMM_couple
!!
!! This modules provides the routines to perform the ensemble
!! communication between the filter processes and the model processes
!! within the communicator COMM_couple.
!!
!! Before the integration the ensemble is scattered using
!! PDAF_scatter_ens, while after the model intregration the
!! ensemble is gathered on the filter processes using
!! PDAF_gather_ens. 
!!
!! The implementation collects the ensemble on rank=0 of each
!! COMM_couple. This assumes usaully that the analysis is computed
!! on one model task. For other cases, e.g. for distributing the
!! analysis step over all model tasks one would need to change
!! the communication scheme.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code from restructuring
!! * Other revisions - see repository log
!!
MODULE PDAF_communicate_ens

CONTAINS

!-------------------------------------------------------------------------------
!> Gather distributed ensemble on filter PEs
!!
!! If the ensemble integration is distributed over multiple
!! model tasks, this routine collects the distributed
!! ensmble information onto the processes that perform
!! the analysis step (filterpe==.true.).
!!
!! __Revision history:__
!! * 2011-12 - Lars Nerger - Initial code extracted from PDAF_put_state_seik
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_ens(dim_p, dim_ens_p, ens, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_mod_parallel, &
         ONLY: mype_filter, mype_couple, npes_couple, filterpe, &
         all_dim_ens_l, all_dis_ens_l, COMM_couple, MPIerr, &
         filter_no_model
    USE PDAF_iau, &
         ONLY: PDAF_iau_update_ens

    IMPLICIT NONE
  
! *** Arguments *** 
    INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens_p   !< Size of ensemble
    REAL, INTENT(inout) :: ens(:, :)   !< PE-local state ensemble
    INTEGER, INTENT(in) :: screen      !< Verbosity flag

! local variables
    INTEGER :: pe_rank, col_frst, col_last  ! Counters
    INTEGER, ALLOCATABLE :: MPIreqs(:)      ! Array of MPI requests
    INTEGER, ALLOCATABLE :: MPIstats(:,:)   ! Array of MPI statuses


! ********************************************
! *** Store local ensembles in case of IAU ***
! ********************************************

    CALL PDAF_iau_update_ens(ens)


! **********************************************
! *** Gather forecast ensemble on filter PEs ***
! **********************************************

    IF (filterpe .AND. mype_filter == 0 .AND. screen > 0 .AND. npes_couple > 1) &
         WRITE (*, '(a, 5x, a)') 'PDAF', '--- Gather sub-ensembles on filter task'

    ! *** call timer
    CALL PDAF_timeit(4, 'new')

    ! *** Send from model PEs that are not filter PEs ***
    subensS: IF (.NOT.filterpe .AND. npes_couple > 1) THEN

        ! Send sub-ensembles to couple PEs with rank 0
       CALL MPI_SEND(ens, dim_p * dim_ens_p, MPI_REALTYPE, 0, mype_couple, &
            COMM_couple, MPIerr)

       IF ((screen>2)) WRITE (*,*) 'PDAF: put_state - send subens of size ', &
            dim_ens_p,' from rank(couple) ',mype_couple, &
            ' in couple task ', mype_filter+1

    END IF subensS

    ! *** Receive on filter PEs ***
    subensR: IF (filterpe .AND. npes_couple > 1) THEN

       ALLOCATE(MPIreqs(npes_couple-1))
       ALLOCATE(MPIstats(MPI_STATUS_SIZE, npes_couple-1))

       ! Receive sub-ensembles on filter PEs
       FnM: IF (filter_no_model) THEN
          taskloopB: DO pe_rank = 1, npes_couple - 1
             col_frst = all_dis_ens_l(pe_rank) + 1
             col_last = col_frst + all_dim_ens_l(pe_rank) - 1 

#ifdef BLOCKING_MPI_EXCHANGE
             CALL MPI_Recv(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIstatus, MPIerr)
#else
             CALL MPI_Irecv(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIreqs(pe_rank), MPIerr)
#endif

             IF (screen > 2) &
                  WRITE (*,*) 'PDAF: put_state - recv subens members ', &
                  col_frst,' to ', col_last,' from rank(couple): ',pe_rank, &
                  ' in couple task ', mype_filter+1
          END DO taskloopB
       ELSE FnM
          taskloopC: DO pe_rank = 1, npes_couple - 1
             col_frst = all_dis_ens_l(pe_rank + 1) + 1
             col_last = col_frst + all_dim_ens_l(pe_rank + 1) - 1 

#ifdef BLOCKING_MPI_EXCHANGE
             call MPI_recv(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank + 1), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIstatus, MPIerr)
#else
             CALL MPI_Irecv(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank + 1), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIreqs(pe_rank), MPIerr)
#endif
             IF (screen > 2) &
                  WRITE (*,*) 'PDAF: put_state - recv subens members ', &
                  col_frst,' to ', col_last,' from rank(couple): ',pe_rank, &
                  ' in couple task ', mype_filter+1
          END DO taskloopC
       END IF FnM

#ifndef BLOCKING_MPI_EXCHANGE
       ! Check for completion of receives
       CALL MPI_Waitall(npes_couple-1, MPIreqs, MPIstats, MPIerr)
#endif

       CALL PDAF_timeit(4, 'old')

       DEALLOCATE(MPIreqs, MPIstats)
     
       IF (screen > 2) &
            WRITE (*,*) 'PDAF: put_state - recv in couple task ', mype_filter+1, ' completed'
    END IF subensR

  END SUBROUTINE PDAF_gather_ens

!-------------------------------------------------------------------------------
!> Gather ensemble to model PEs
!!
!! If the ensemble integration is distributed over multiple
!! model tasks, this routine distributes the ensemble from
!! the processes that perform the analysis step (filterpe==.true.)
!! to all model processes.
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code extracted from PDAF_get_state
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_scatter_ens(dim_p, dim_ens_p, ens, state, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_mod_parallel, &
         ONLY: mype_filter, mype_couple, npes_couple, filterpe, &
         all_dim_ens_l, all_dis_ens_l, COMM_couple, MPIerr, &
         filter_no_model, MPIstatus, statetask
    USE PDAF_mod_core, &
         ONLY: ensemblefilter
    USE PDAF_iau, &
         ONLY: PDAF_iau_update_inc

    IMPLICIT NONE
  
! *** Arguments *** 
    INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens_p   !< Size of ensemble
    REAL, INTENT(inout) :: ens(:, :)   !< PE-local state ensemble
    REAL, INTENT(inout) :: state(:)    !< PE-local state vector (for SEEK)
    INTEGER, INTENT(in) :: screen      !< Verbosity flag

! local variables
    INTEGER :: pe_rank, col_frst, col_last  ! Counters
    INTEGER, ALLOCATABLE :: MPIreqs(:)      ! Array of MPI requests
    INTEGER, ALLOCATABLE :: MPIstats(:,:)   ! Array of MPI statuses
    INTEGER :: dummy                        ! Dummy variable


! *************************************************
! *** Scatter forecast ensemble from filter PEs ***
! *************************************************

    ! Dummy initialization preventing compiler warning
    dummy = dim_ens_p

    CALL PDAF_timeit(4, 'new')

    ! *** Send from filter PEs ***
    subensS: IF (filterpe .AND. npes_couple > 1) THEN

       IF (mype_filter == 0 .AND. screen > 0) &
            WRITE (*, '(a, 5x, a)') 'PDAF', '--- Distribute sub-ensembles'

       ALLOCATE(MPIreqs(npes_couple-1))
       ALLOCATE(MPIstats(MPI_STATUS_SIZE, npes_couple-1))

       ! Send sub-ensembles to each model PE within coupling communicator
       FnM: IF (filter_no_model) THEN
          taskloopB: DO pe_rank = 1, npes_couple - 1
             col_frst = all_dis_ens_l(pe_rank) + 1
             col_last = col_frst + all_dim_ens_l(pe_rank) - 1 

#ifdef BLOCKING_MPI_EXCHANGE
             CALL MPI_Send(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIerr)
#else
             CALL MPI_Isend(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIreqs(pe_rank), MPIerr)
#endif

             IF (screen > 2) &
                  WRITE (*,*) 'PDAF: get_state - send subens members ', &
                  col_frst,' to ', col_last,' to rank(couple): ',pe_rank, &
                  ' in couple task ', mype_filter+1
          END DO taskloopB

          ! SEEK: Send central state to STATETASK
          ifSEEK1: IF (.NOT.ensemblefilter) THEN

             CALL MPI_SEND(state, dim_p, MPI_REALTYPE, &
                  statetask-1, statetask - 1, COMM_couple, MPIerr)

             IF (screen > 2) WRITE (*,*) &
                  'PDAF: get_state - send state to statetask ',statetask, &
                  ' in couple task ', mype_filter + 1
          END IF ifSEEK1
       ELSE
          taskloopC: DO pe_rank = 1, npes_couple - 1
             col_frst = all_dis_ens_l(pe_rank + 1) + 1
             col_last = col_frst + all_dim_ens_l(pe_rank + 1) - 1 

#ifdef BLOCKING_MPI_EXCHANGE
             CALL MPI_Send(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank + 1), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIerr)
#else
             CALL MPI_Isend(ens(1, col_frst), &
                  dim_p * all_dim_ens_l(pe_rank + 1), MPI_REALTYPE, &
                  pe_rank, pe_rank, COMM_couple, MPIreqs(pe_rank), MPIerr)
#endif

             IF (screen > 2) &
                  WRITE (*,*) 'PDAF: get_state - send subens members ', &
                  col_frst,' to ', col_last,' to rank(couple): ',pe_rank, &
                  ' in couple task ', mype_filter+1
          END DO taskloopC

          ! SEEK: Send central state to STATETASK
          ifSEEK3: IF ((.NOT.ensemblefilter) .AND. statetask > 1) THEN

             CALL MPI_SEND(state, dim_p, MPI_REALTYPE, &
                  statetask - 1, statetask - 1, COMM_couple, MPIerr)

             IF (screen > 2) WRITE (*,*) &
                  'PDAF: get_state - send state to statetask ',statetask, &
                  ' in couple task ', mype_filter + 1
          END IF ifSEEK3
       END IF FnM

#ifndef BLOCKING_MPI_EXCHANGE
       ! Check for completion of sends
       CALL MPI_Waitall(npes_couple-1, MPIreqs, MPIstats, MPIerr)
#endif

       DEALLOCATE(MPIreqs, MPIstats)

       IF (screen > 2) &
            WRITE (*,*) 'PDAF: get_state - send in couple task ', mype_filter+1, ' completed'

    END IF subensS

    ! *** Receive on model PEs that are not filter PEs ***
    subensRA: IF (.NOT.filterpe .AND. npes_couple > 1) THEN
       FnMA: IF (filter_no_model) THEN

          ! Receive sub-ensemble on each model PE 0
          CALL MPI_RECV(ens, dim_p * all_dim_ens_l(mype_couple), &
               MPI_REALTYPE, 0, mype_couple, COMM_couple, MPIstatus, MPIerr)
          IF (screen > 2) &
               WRITE (*,*) 'PDAF: get_state - recv subens of size ', &
               all_dim_ens_l(mype_couple),' on rank(couple) ',mype_couple, &
               ' in couple task ', mype_filter+1

          ! SEEK: Receive central state on model PE 0 of STATETASK
          ifSEEK4: IF ((.NOT.ensemblefilter) .AND. mype_couple == statetask) THEN
                 
             CALL MPI_RECV(state, dim_p, MPI_REALTYPE, &
                  0, mype_couple, COMM_couple, MPIstatus, MPIerr)
             IF (screen > 2) WRITE (*,*) &
                  'PDAF: get_state - recv state on statetask ', &
                  statetask, ' in couple task ', mype_filter + 1
          END IF ifSEEK4

       ELSE

          ! Receive sub-ensemble on each model PE 0
          CALL MPI_RECV(ens, dim_p * all_dim_ens_l(mype_couple + 1), &
               MPI_REALTYPE, 0, mype_couple, COMM_couple, MPIstatus, MPIerr)
          IF (screen > 2) &
               WRITE (*,*) 'PDAF: get_state - recv subens of size ', &
               all_dim_ens_l(mype_couple+1),' on rank(couple) ',mype_couple, &
               ' in couple task ', mype_filter+1

          ! SEEK: Receive central state on model PE 0 of STATETASK
          ifSEEK2: IF ((.NOT.ensemblefilter) .AND. mype_couple+1 == statetask) THEN
                 
             CALL MPI_RECV(state, dim_p, MPI_REALTYPE, &
                  0, mype_couple, COMM_couple, MPIstatus, MPIerr)
             IF (screen > 2) WRITE (*,*) &
                  'PDAF: get_state - recv state on statetask ', &
                  statetask, ' in couple task ', mype_filter + 1
          END IF ifSEEK2

       END IF FnMA
    END IF subensRA

    CALL PDAF_timeit(4, 'old')


! *****************************************
! *** Update ensemble increment for IAU ***
! *****************************************

    CALL PDAF_iau_update_inc(ens)

  END SUBROUTINE PDAF_scatter_ens

END MODULE PDAF_communicate_ens
