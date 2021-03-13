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
! !ROUTINE: PDAF_3dvar_analysis_cvt --- incr. 3DVAR with variable transformation
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens, &
     state_p, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
     screen, incremental, flag)

! !DESCRIPTION:
! Analysis step of incremental 3DVAR with control variable transformation.
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
       U_prodRinvA              ! Provide product R^-1 A

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
  INTEGER :: converged                 ! Flag whether the optimization converged
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: maxiter                   ! Maximum number of minimization iterations
  REAL :: eps_min                      ! Limit of change of J to stop iterations
  REAL :: invdimens                    ! Inverse global ensemble size
  REAL, ALLOCATABLE :: HV_p(:,:)       ! observed ensemble perturbations
  REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
  REAL, ALLOCATABLE :: deltay_p(:)     ! PE-local observation background residual
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local analysis increment vector
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local produce HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL :: J_tot, J_old                 ! Cost function terms
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  INTEGER :: incremental_dummy         ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)         ! Dummy variable to avoid compiler warning
  REAL :: fact                         ! Scaling factor for transforming from v_p to x_p


! **********************
! *** INITIALIZATION ***
! **********************

  ! Settings for steepest descent
  maxiter = 200
  eps_min = 1.0e-7
  converged = 0

  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)


  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - 3DVAR incremental, transformed'
     IF (type_opt==0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: LBFGS' 
     ELSEIF (type_opt==1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: CG+' 
     ELSEIF (type_opt==2) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: plain CG' 
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- solver: steepest descent' 
        WRITE (*, '(a, 8x, a, es10.2)') 'PDAF', '--- Convergence limit:', eps_min 
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
     ALLOCATE(deltay_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     obs_member = 0 ! Store member index (0 for central state)
     CALL PDAF_timeit(44, 'new')
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, deltay_p)
     CALL PDAF_timeit(44, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     deltay_p = obs_p - deltay_p
     CALL PDAF_timeit(51, 'old')


     ! *** Observated background ensemble perturbations ***

     ! Get observed ensemble
     ALLOCATE(HV_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     CALL PDAF_timeit(44, 'new')
     ENS1: DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        ! [Hx_1 ... Hx_N]
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HV_p(:, member))
     END DO ENS1
     CALL PDAF_timeit(44, 'old')

     HV_p = HV_p / SQRT(REAL(dim_ens-1))

     ! Compute observed ensemble perturbations
     ! Subtract ensemble mean: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_etkf_Tright(dim_obs_p, dim_ens, HV_p)
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(12, 'old')


! ****************************
! ***   Iterative solving  ***
! ****************************
     
     ! Prepare control vector for optimization
     ALLOCATE(v_p(dim_ens))
     v_p = 0.0


     opt: IF (type_opt==0) THEN

        ! LBFGS solver
        CALL PDAF_3dvar_optim_lbfgs(step, dim_ens, dim_obs_p, &
             obs_p, deltay_p, HV_p, v_p, U_prodRinvA, screen)

     ELSEIF (type_opt==1) THEN

        ! CG+ solver
        CALL PDAF_3dvar_optim_cgplus(step, dim_ens, dim_obs_p, &
             obs_p, deltay_p, HV_p, v_p, U_prodRinvA, screen)

     ELSEIF (type_opt==2) THEN

        ! CG solver
        CALL PDAF_3dvar_optim_cg(step, dim_ens, dim_obs_p, &
             obs_p, deltay_p, HV_p, v_p, U_prodRinvA, screen)

     ELSE

        ! *** Simple steepest descent

        ! Allocate arrays for iterations
        ALLOCATE(HVv_p(dim_obs_p))
        ALLOCATE(RiHVv_p(dim_obs_p, 1))
        ALLOCATE(gradJ_p(dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_ens + 2*dim_obs_p)


        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(a, 5x, a)') 'PDAF', '--- Optimizing ...' 

        minloop: DO iter = 1, maxiter

           ! *** Evaluate cost function ***

           IF (iter>1) J_old = J_tot

           CALL PDAF_3dvar_costf_cvt(step, dim_ens, dim_obs_p, &
                obs_p, deltay_p, HV_p, v_p, J_tot, gradJ_p, &
                U_prodRinvA, screen)

           IF (mype==0 .AND. screen > 2) &
                WRITE (*,'(a, 8x, a, i5, es12.4)') 'PDAF', '--- iter, J: ', iter, J_tot


           ! *** Update control vector ***

           v_p = v_p - 0.1 * gradJ_p


           ! *** Check convergence ***

           IF (iter>1 .AND. J_old - J_tot < eps_min) THEN
              converged = 1
              EXIT minloop
           END IF

        END DO minloop

        IF (mype==0 .AND. screen >0) THEN
           IF (converged==1) THEN
              WRITE (*,'(a, 8x, a, i5, a)') &
                   'PDAF', '--- optimization converged after ', iter,' iterations'
           ELSE
              WRITE (*,'(a, 8x, a, i5, a)') &
                   'PDAF', '--- optimization exited after maxiter ', maxiter,' iterations'
           END IF
        END IF

        DEALLOCATE(HVv_p, RiHVv_p, gradJ_p)

     END IF opt


! ***************************************************
! ***   Solving completed: Update state estimate  ***
! ***************************************************

     CALL PDAF_timeit(13, 'new')

     ! Get ensemble perturbation matrix X'=X-xmean subtract forecast state
     DO col = 1, dim_ens
        DO row = 1, dim_p
           ens_p(row, col) = ens_p(row, col) - state_p(row)
        END DO
     END DO

     fact = 1.0/SQRT(REAL(dim_ens-1))

     ! Transform control variable to state increment
     CALL gemvTYPE('n', dim_p, dim_ens, fact, ens_p, &
          dim_p, v_p, 1, 1.0, state_p, 1)

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
     DEALLOCATE(obs_p, deltay_p, HV_p)
     DEALLOCATE(v_p)
  END IF

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_analysis_cvt
