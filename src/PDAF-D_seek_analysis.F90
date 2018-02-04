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
!$Id: PDAF-D_seek_analysis.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_analysis --- Perform SEEK analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_seek_analysis(step, dim_p, dim_obs_p, dim_eof, state_p, &
     Uinv, V_p, forget, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, screen, incremental, type_forget, &
     flag)

! !DESCRIPTION:
! Analysis step of the SEEK filter
! with forgetting factor.
!
! Variant for domain decomposed states.
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
       ONLY: state_inc, obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_eof      ! Number of EOFs
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  ! ! *** The covariance P is decomposed as P = V U V^T ***
  REAL, INTENT(inout) :: Uinv(dim_eof,dim_eof) ! Inverse of matrix U
  REAL, INTENT(inout) :: V_p(dim_p,dim_eof)    ! PE-local matrix V
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA              ! Provide product Rinv A for SEEK analysis

! !CALLING SEQUENCE:
! Called by: PDAF_seek_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: j                        ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HV_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHV_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_p(:,:)    ! local Uinv
  REAL, ALLOCATABLE :: Uinv_inc(:,:)  ! increment for Uinv
  REAL, ALLOCATABLE :: resid_p(:)     ! observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! observation vector
  REAL, ALLOCATABLE :: m_state_p(:)   ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHVd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHVd_p(:)     ! local RiHVd
  REAL, ALLOCATABLE :: temp_Uinv(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)     ! vector of pivot indices for GESV
  INTEGER :: gesv_info              ! Control flag for GESV

  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') 'PDAF ', step, 'Assimilating observations - SEEK'
  END IF


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL U_init_dim_obs(step, dim_obs_p)

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')

  haveobsB: IF (dim_obs_p > 0) THEN
     ! The residual only exists for domains with observations

     ALLOCATE(resid_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(m_state_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_p)

     ! *** Project state onto observation space and    ***
     ! *** compute observation residual (innovation) d ***

     ! Project state onto observation space
     obs_member = 0 ! Store member index (0 for central state)
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, m_state_p)

     ! get observation vector
     CALL U_init_obs(step, dim_obs_p, obs_p)

     ! get residual as difference of observation an projected state
     resid_p = obs_p - m_state_p

     DEALLOCATE(m_state_p)
  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! *****************************************
! ***   Compute analyzed matrix Uinv    ***
! ***                                   ***
! ***  -1          -1    T  T  -1       ***
! *** U  = forget*U   + V  H  R   H  V  ***
! ***  i           i-1   i  i  i   i  i ***
! *****************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! HV = H V
     ALLOCATE(HV_p(dim_obs_p, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_eof)

     DO j = 1, dim_eof
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = j

        ! Call observation operator
        CALL U_obs_op(step, dim_p, dim_obs_p, V_p(:, j), HV_p(:, j))
     ENDDO

     CALL PDAF_timeit(30,'old')

     ! ***                RiHV = Rinv HV                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     ALLOCATE(RiHV_p(dim_obs_p, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_eof)

     CALL U_prodRinvA(step, dim_obs_p, dim_eof, obs_p, HV_p, RiHV_p)

     ! *** Finish computation of Uinv  ***
     ! ***   -1          -1    T       ***
     ! ***  U  = forget U  + HV  RiHV  ***
     ALLOCATE(Uinv_p(dim_eof, dim_eof))
     ALLOCATE(Uinv_inc(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof * dim_eof)

     CALL PDAF_timeit(31, 'new')
     ! partial increment
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_obs_p, &
          1.0, HV_p, dim_obs_p, RiHV_p, dim_obs_p, &
          0.0, Uinv_p, dim_eof)

     DEALLOCATE(HV_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Uinv  ***
 
     CALL PDAF_timeit(31, 'new')

     ALLOCATE(Uinv_p(dim_eof, dim_eof))
     ALLOCATE(Uinv_inc(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof * dim_eof)

     ! No observation-contribution to Uinv from this domain
     Uinv_p = 0.0
  END IF haveobsA

  ! get total increment on all filter PEs
  CALL MPI_allreduce(Uinv_p, Uinv_inc, dim_eof * dim_eof, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  Uinv = forget * Uinv + Uinv_inc

  DEALLOCATE(Uinv_p, Uinv_inc)

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! *******************************************************
! ***             update model state                  ***
! ***                                                 ***
! ***  a   f            f     f          T        f   ***
! *** x = x + K (y - H x ) = x + V U RiHV (y - H x )  ***
! ***                                                 ***
! *******************************************************


  ! ************************
  ! *** RiHVd = RiHV^T d ***
  ! ************************

  CALL PDAF_timeit(13, 'new')

  ALLOCATE(RiHVd(dim_eof))
  ALLOCATE(RiHVd_p(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, dim_eof, 1.0, RiHV_p, &
          dim_obs_p, resid_p, 1, 0.0, RiHVd_p, 1)

     DEALLOCATE(RiHV_p, resid_p)
  ELSE haveobsC
     RiHVd_p = 0.0
  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHVd_p, RiHVd, dim_eof, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHVd_p)


  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                    ***
  ! ****************************************

  ALLOCATE(temp_Uinv(dim_eof, dim_eof))
  ALLOCATE(ipiv(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_eof * dim_eof)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', dim_eof)

  ! save matrix Uinv
  temp_Uinv = Uinv

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(dim_eof, 1, temp_Uinv, dim_eof, ipiv, RiHVd, dim_eof, gesv_info)
  DEALLOCATE(temp_Uinv, ipiv)

  CALL PDAF_timeit(13, 'old')

  ! *** check whether solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in solve for state analysis !!!'
     flag = 1
  ELSE

     CALL PDAF_timeit(14, 'new')

     ! **************************
     ! *** Update model state ***
     ! ***     a   f          ***
     ! ***   x = x + V RiHVd  ***
     ! **************************

     IF (incremental == 0) THEN
        ! Allocate only if no incremental updating is used. 
        ! With incremental STATE_INC is allocated in filter_init.
        ALLOCATE(state_inc(dim_p))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_p)
     END IF

     CALL gemvTYPE('n', dim_p, dim_eof, 1.0, V_p, &
          dim_p, RiHVd, 1, 0.0, state_inc, 1)
     DEALLOCATE(RiHVd)

     IF (incremental == 0) THEN
        ! update state only if incremental updating is not used
        state_p = state_p + state_inc
        DEALLOCATE(state_inc)
     END IF

     CALL PDAF_timeit(14, 'old')
        
  END IF update

    
! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_seek_analysis
