! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_seik_analysis.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seik_analysis --- Perform SEIK analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_seik_analysis(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Uinv, ens_p, state_inc_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, flag)

! !DESCRIPTION:
! Analysis step of the SEIK filter
! with adaptive forgetting factor.
!
! Variant for domain decomposed states. 
! Old formulation regarding application of T.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: filterstr, obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local model state
  REAL, INTENT(inout) :: Uinv(rank, rank)      ! Inverse of eigenvalue matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obsvar, &         ! Initialize mean observation error variance
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA              ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_seik_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_set_forget
! Calls: PDAF_seik_matrixT
! Calls: PDAF_seik_Uinv
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP
       
! *** local variables ***
  INTEGER :: member, row              ! counters
  REAL :: invdimens                   ! Inverse global ensemble size
  REAL :: forget_ana                  ! forgetting factor used for analysis
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HL_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_p(:,:)    ! local Uinv
  REAL, ALLOCATABLE :: Uinv_inc(:,:)  ! increment for Uinv
  REAL, ALLOCATABLE :: resid_p(:)     ! observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! observation vector
  REAL, ALLOCATABLE :: m_state_p(:)   ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHLd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)     ! local RiHLd

  ! *** Variables for LU solver GESV
  REAL, ALLOCATABLE :: temp_Uinv(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)     ! vector of pivot indices for GESV
  INTEGER :: gesv_info                ! control flag for GESV


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Assimilating observations - SEIK old formulation'
  END IF


! ***********************************
! *** Compute mean forecast state ***
! ***********************************
  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')
  

! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


  ! ************************
  ! *** Compute residual ***
  ! ***   d = y - H x    ***
  ! ************************

  CALL PDAF_timeit(12, 'new')

  haveobsB: IF (dim_obs_p > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ALLOCATE(resid_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(m_state_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_obs_p)
     
     ! Project state onto observation space
     obs_member = 0 ! Store member index (0 for central state)
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, m_state_p)

     ! get observation vector
     CALL U_init_obs(step, dim_obs_p, obs_p)

     ! get residual as difference of observation and
     ! projected state
     resid_p = obs_p - m_state_p

  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! **********************************************
! ***   Compute analyzed matrix Uinv         ***
! ***                                        ***
! ***  -1                T        T  -1      ***
! *** U  = forget*(r+1) T T + (HL)  R  (HL)  ***
! ***  i                          i  i     i ***
! **********************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! HL = [Hx_1 ... Hx_(r+1)] T 
     ALLOCATE(HL_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     ENS: DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        ! [Hx_1 ... Hx_(r+1)]
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HL_p(:, member))
     END DO ENS

     ! Set forgetting factor
     forget_ana = forget
     IF (type_forget == 1) THEN
        CALL PDAF_set_forget(step, filterstr, dim_obs_p, dim_ens, HL_p, &
             m_state_p, obs_p, U_init_obsvar, forget, forget_ana)
     ENDIF
     DEALLOCATE(m_state_p)

     ! HL = [Hx_1 ... Hx_(r+1)] T
     CALL PDAF_seik_matrixT(dim_obs_p, dim_ens, HL_p)

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     ALLOCATE(RiHL_p(dim_obs_p, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * rank)

     CALL U_prodRinvA(step, dim_obs_p, rank, obs_p, HL_p, RiHL_p)
     DEALLOCATE(obs_p)

     ! *** Initialize Uinv = (r+1) T^T T ***
     CALL PDAF_seik_Uinv(rank, Uinv)


     ! *** Finish computation of Uinv  ***
     ! ***   -1          -1    T       ***
     ! ***  U  = forget U  + HL RiHL   ***
     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(Uinv_inc(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_p, &
          1.0, HL_p, dim_obs_p, RiHL_p, dim_obs_p, &
          0.0, Uinv_p, rank)

     DEALLOCATE(HL_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Uinv  ***
 
     CALL PDAF_timeit(31, 'new')

     ! Set forgetting factor
     forget_ana = forget

     ! Initialize Uinv = N T^T T 
     CALL PDAF_seik_Uinv(rank, Uinv)

     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(Uinv_inc(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     ! No observation-contribution to Uinv from this domain
     Uinv_p = 0.0
  END IF haveobsA

  ! get total sum on all filter PEs
  CALL MPI_allreduce(Uinv_p, Uinv_inc, rank**2, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  Uinv = forget_ana * Uinv + Uinv_inc

  DEALLOCATE(Uinv_p, Uinv_inc)

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ************************************
! ***      update model state      ***
! ***                              ***
! ***  a   f          T         f  ***
! *** x = x + L U RiHV  (y - H x ) ***
! ***                              ***
! ************************************


  ! ************************************************
  ! *** Compute matrix L = [x_1, ..., x_(r+1)] T ***
  ! ************************************************
  CALL PDAF_timeit(16, 'new')
  CALL PDAF_seik_matrixT(dim_p, dim_ens, ens_p)
  CALL PDAF_timeit(16, 'old')

  CALL PDAF_timeit(13, 'new')
  ! ************************
  ! *** RiHLd = RiHV^T d ***
  ! ************************
  ALLOCATE(RiHLd_p(rank))
  ALLOCATE(RiHLd(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, rank, 1.0, RiHL_p, &
          dim_obs_p, resid_p, 1, 0.0, RiHLd_p, 1)
     DEALLOCATE(RiHL_p, resid_p)
  ELSE haveobsC
     RiHLd_p = 0.0
  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHLd_p, RiHLd, rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHLd_p)

 
  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                   ***
  ! ****************************************

  ALLOCATE(temp_Uinv(rank, rank))
  ALLOCATE(ipiv(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

  ! save matrix Uinv
  temp_Uinv = Uinv

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(rank, 1, temp_Uinv, rank, ipiv, &
       RiHLd, rank, gesv_info)
  DEALLOCATE(temp_Uinv, ipiv)

  CALL PDAF_timeit(13, 'old')

  ! *** check if solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x,a/)') 'PDAF-ERROR(1): Problem in solve for state analysis !!!'
     flag = 1
  ELSE

     CALL PDAF_timeit(14, 'new')

     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + L RiHLd  ***
     ! **************************

     CALL gemvTYPE('n', dim_p, rank, 1.0, ens_p, &
          dim_p, RiHLd, 1, 0.0, state_inc_p, 1)
     DEALLOCATE(RiHLd)
    
     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_p = state_p + state_inc_p
     END IF

     CALL PDAF_timeit(14, 'old')
    
  END IF update


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_seik_analysis
