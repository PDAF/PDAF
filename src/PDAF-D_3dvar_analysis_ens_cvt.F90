! Copyright (c) 2004-2020 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_3dvar_analysis_ens_cvt --- ensemble 3DVAR with CVT
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_analysis_ens_cvt(step, dim_p, dim_obs_p, dim_ens, &
     dim_cvec_ens, state_p, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
     U_cvtmat_ens, screen, incremental, flag)

! !DESCRIPTION:
! Analysis step of incremental ensemble 3DVAR with control
! variable transformation.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
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
  USE PDAF_mod_filtermpi, &
       ONLY: mype !, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE
  USE PDAF_mod_filter, &
       ONLY: obs_member, type_opt

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: dim_cvec_ens ! Size of control vector (ensemble part)
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)      ! PE-local state analysis increment
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA, &           ! Provide product R^-1 A
       U_cvtmat_ens             ! Initialize CVT transform matrix in obs. space

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_ETKF_Tright
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
!EOP

! *** local variables ***
  INTEGER :: iter, member, col, row    ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: invdimens                    ! Inverse global ensemble size
  REAL, ALLOCATABLE :: HV_p(:,:)       ! observed ensemble perturbations
  REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
  REAL, ALLOCATABLE :: dy_p(:)         ! PE-local observation background residual
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local analysis increment vector
  INTEGER :: incremental_dummy         ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)         ! Dummy variable to avoid compiler warning
  REAL :: fact                         ! Scaling factor for transforming from v_p to x_p


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)


  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - ensemble 3DVAR incremental, transformed'
     IF (type_opt==0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: LBFGS' 
     ELSEIF (type_opt==1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: CG+' 
     ELSEIF (type_opt==2) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: plain CG' 
     END IF
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
  CALL PDAF_timeit(51, 'new')


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')
  
  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


  haveobsB: IF (dim_obs_p > 0) THEN

! ********************************************************
! *** Background innovation and ensemble perturbations ***
! ***          d = y - H xb,  H (X - meanX) **         ***
! ********************************************************

     CALL PDAF_timeit(12, 'new')
  
     ! *** Observation background innovation ***

     ! get observation vector
     ALLOCATE(obs_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! Get observed state estimate
     ALLOCATE(dy_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     obs_member = 0 ! Store member index (0 for central state)
     CALL PDAF_timeit(44, 'new')
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, dy_p)
     CALL PDAF_timeit(44, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     dy_p = obs_p - dy_p
     CALL PDAF_timeit(51, 'old')


     ! *** Observed localized/scaled background ensemble perturbations ***

     ! Get observed ensemble perturbations
     ALLOCATE(HV_p(dim_obs_p, dim_cvec_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     ! Get ensemble perturbation matrix X'=X-xmean subtract forecast state
     DO col = 1, dim_ens
        DO row = 1, dim_p
           ens_p(row, col) = ens_p(row, col) - state_p(row)
        END DO
     END DO

     CALL PDAF_timeit(44, 'new')
     CALL U_cvtmat_ens(step, dim_p, dim_ens, dim_obs_p, dim_cvec_ens, &
          ens_p(:, :), HV_p(:, :))
     CALL PDAF_timeit(44, 'old')

     CALL PDAF_timeit(12, 'old')


! ****************************
! ***   Iterative solving  ***
! ****************************
     
     ! Prepare control vector for optimization
     ALLOCATE(v_p(dim_cvec_ens))
     v_p = 0.0


     opt: IF (type_opt==0) THEN

        ! LBFGS solver
        CALL PDAF_3dvar_optim_lbfgs(step, dim_cvec_ens, dim_obs_p, &
             obs_p, dy_p, HV_p, v_p, U_prodRinvA, screen)

     ELSEIF (type_opt==1) THEN

        ! CG+ solver
        CALL PDAF_3dvar_optim_cgplus(step, dim_cvec_ens, dim_obs_p, &
             obs_p, dy_p, HV_p, v_p, U_prodRinvA, screen)

     ELSEIF (type_opt==2) THEN

        ! CG solver
        CALL PDAF_3dvar_optim_cg(step, dim_cvec_ens, dim_obs_p, &
             obs_p, dy_p, HV_p, v_p, U_prodRinvA, screen)

     ELSE
        ! Further solvers - not implemented
        WRITE (*,*) 'PDAF', 'INVALID SOLVER!'
        flag=10
     END IF opt


! ***************************************************
! ***   Solving completed: Update state estimate  ***
! ***************************************************

     CALL PDAF_timeit(13, 'new')

     CALL cvec2state_ens_pdaf(step, dim_p, dim_ens, dim_cvec_ens, &
          ens_p(:, :), v_p, state_p)

     ! Add analysis state to ensemble perturbations
     DO col = 1, dim_ens
        DO row = 1, dim_p
           ens_p(row, col) = ens_p(row, col) + state_p(row)
        END DO
     END DO

     CALL PDAF_timeit(13, 'old')

  END IF haveobsB


! ********************
! *** Finishing up ***
! ********************

  IF (dim_obs_p > 0) THEN
     DEALLOCATE(obs_p, dy_p, HV_p)
     DEALLOCATE(v_p)
  END IF

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_analysis_ens_cvt
