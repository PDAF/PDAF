! Copyright (c) 2004-2021 Lars Nerger
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
! !ROUTINE: PDAF_3dvar_analysis_cvt --- 3DVAR with CVT
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_cvec, &
     state_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     screen, incremental, type_opt, flag)

! !DESCRIPTION:
! Analysis step of incremental 3DVAR with control
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
       ONLY: mype, comm_filter, mpierr
  USE PDAF_mod_filter, &
       ONLY: obs_member

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_cvec   ! Size of control vector
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(inout) :: state_inc_p(dim_p)      ! PE-local state analysis increment
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_opt     ! Type of minimizer for 3DVar
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA, &           ! Provide product R^-1 A
       U_cvt, &                 ! Apply control vector transform matrix to control vector
       U_cvt_adj, &             ! Apply adjoint control vector transform matrix
       U_obs_op_lin, &          ! Linearized observation operator
       U_obs_op_adj             ! Adjoint observation operator

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
  REAL, ALLOCATABLE :: dy_p(:)         ! PE-local observation background residual
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local analysis increment vector
  INTEGER :: opt_parallel              ! Whether to run solver with decomposed control vector


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     IF (type_opt==1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: LBFGS' 
     ELSEIF (type_opt==2) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: CG+' 
     ELSEIF (type_opt==3) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: plain CG' 
     ELSEIF (type_opt==12) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: CG+ parallelized' 
     ELSEIF (type_opt==13) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: plain CG parallelized' 
     END IF
  END IF


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(43, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(43, 'old')
  
  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


  haveobsB: IF (dim_obs_p > 0) THEN

! *******************************
! *** Background innovation   ***
! ***      d = y - H xb       ***
! *******************************

     CALL PDAF_timeit(12, 'new')
  
     ! *** Observation background innovation ***

     ! Get observed state estimate
     ALLOCATE(dy_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     obs_member = 0 ! Store member index (0 for central state)
     CALL PDAF_timeit(44, 'new')
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, dy_p)
     CALL PDAF_timeit(44, 'old')

     ! get observation vector
     ALLOCATE(obs_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     dy_p = obs_p - dy_p
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(12, 'old')


! ****************************
! ***   Iterative solving  ***
! ****************************

     CALL PDAF_timeit(52, 'new')

     opt_parallel = 0

     ! Prepare control vector for optimization
     ALLOCATE(v_p(dim_cvec))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_cvec)
     v_p = 0.0

     ! Choose solver
     opt: IF (type_opt==1) THEN

        ! LBFGS solver
        CALL PDAF_3dvar_optim_lbfgs(step, dim_p, dim_cvec, dim_obs_p, &
             obs_p, dy_p, v_p, U_prodRinvA, &
             U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
             opt_parallel, screen)

     ELSEIF (type_opt==2 .OR. type_opt==12) THEN

        IF (type_opt==12) opt_parallel = 1

        ! CG+ solver
        CALL PDAF_3dvar_optim_cgplus(step, dim_p, dim_cvec, dim_obs_p, &
             obs_p, dy_p, v_p, U_prodRinvA, &
             U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
             opt_parallel, screen)

     ELSEIF (type_opt==3 .OR. type_opt==13) THEN

        IF (type_opt==13) opt_parallel = 1

        ! CG solver
        CALL PDAF_3dvar_optim_cg(step, dim_p, dim_cvec, dim_obs_p, &
             obs_p, dy_p, v_p, U_prodRinvA, &
             U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
             opt_parallel, screen)

     ELSE
        ! Further solvers - not implemented
        WRITE (*,'(a)') 'PDAF-ERROR: Invalid solver chosen!'
        flag=10
     END IF opt

     CALL PDAF_timeit(52, 'old')


! ****************************************************
! ***   Solving completed: Compute state increment ***
! ****************************************************

     CALL PDAF_timeit(14, 'new')

     ! State increment: Apply V to control vector v_p
     CALL PDAF_timeit(49, 'new')
     CALL U_cvt(0, dim_p, dim_cvec, v_p, state_inc_p)
     CALL PDAF_timeit(49, 'old')

     CALL PDAF_timeit(51, 'new')
     IF (incremental<1) THEN
        ! Add analysis increment to state vector
        state_p = state_p + state_inc_p
     END IF
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(14, 'old')

  END IF haveobsB


! ********************
! *** Finishing up ***
! ********************

  IF (dim_obs_p > 0) THEN
     DEALLOCATE(obs_p, dy_p)
     DEALLOCATE(v_p)
  END IF

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_analysis_cvt
