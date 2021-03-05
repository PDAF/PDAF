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
! !ROUTINE: PDAF_3dvar_analysis_transf --- incr. 3DVAR with variable transformation
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_analysis_transf(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, state_inc_p, forget, forget_ana, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, flag)

! !DESCRIPTION:
! Analysis step of incremental 3DVAR with variable transformation.
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
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE
  USE PDAF_mod_filter, &
       ONLY: type_trans, filterstr, obs_member, observe_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(out)   :: Ainv(dim_ens, dim_ens)  ! on entry: uninitialized
                                      ! on exit: weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)      ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  REAL, INTENT(out)   :: forget_ana   ! Forgetting factor actually used in analysis
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
! Called by: PDAF_3dvar_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_ETKF_Tright
! Calls: PDAF_generate_rndmat
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, iter, member, col, row ! Counters
  INTEGER :: type_min                  ! Type of minimizer
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
  REAL :: J_B, J_obs, J_tot, J_old     ! Cost function terms
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  INTEGER :: incremental_dummy         ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)         ! Dummy variable to avoid compiler warning
  REAL :: fact                         ! Scaling factor for transforming from v_p to x_p

  ! Variables for LBFG
  INTEGER, PARAMETER :: m = 5
  INTEGER, PARAMETER :: iprint = 0
  CHARACTER(len=60)  :: task, csave
  LOGICAL            :: lsave(4)
  INTEGER            :: isave(44)
  REAL, PARAMETER    :: factr  = 1.0e+7, pgtol  = 1.0e-5
  REAL               :: dsave(29)
  INTEGER,  ALLOCATABLE  :: nbd(:), iwa(:)
  REAL, ALLOCATABLE  :: lvec(:), uvec(:), wa(:)



! **********************
! *** INITIALIZATION ***
! **********************

  ! Type of minimizer
  type_min = 1   ! 0: Steepest descent, 1: LBGFS

  ! Settings for steepest descent
  maxiter = 100
  eps_min = 0.00001

  ! Settings for LBGFS
  ALLOCATE(nbd(dim_ens), lvec(dim_ens), uvec(dim_ens))
  ALLOCATE (iwa(3*dim_ens))
  ALLOCATE (wa(2*m*dim_ens + 5*dim_ens + 11*m*m + 8*m))
  nbd = 0  ! Values are unbounded
  task = 'START'


  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)


  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - 3DVAR transformed - intremental'
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


! ***************************
! ***   Iterative solving ***
! ***************************

     ! Prepare arrays for iterations
     ALLOCATE(v_p(dim_ens))
     ALLOCATE(HVv_p(dim_obs_p))
     ALLOCATE(RiHVv_p(dim_obs_p, 1))
     ALLOCATE(gradJ_p(dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_ens + 2*dim_obs_p)

     ! Initialize control vector
     v_p = 0.0


     WRITE (*,*) 'ITERATE ----------'

     minloop: DO iter = 1, maxiter

        IF (type_min==1) THEN
           IF (.NOT.(task(1:2).EQ.'FG'.OR.task.EQ.'NEW_X'.OR. &
                task.EQ.'START') ) THEN
              WRITE (*,*) 'Exit iterations, status ', task
              EXIT minloop
           END IF
        END IF
        IF (type_min == 0) THEN
           ! Steepest descent
           v_p = v_p - 0.1 * gradJ_p
        ELSE
           ! LBFGS
           CALL setulb (dim_ens, m, v_p, lvec, uvec, nbd, &
                J_tot, gradJ_p, factr, pgtol, &
                wa, iwa, task, iprint,&
                csave, lsave, isave, dsave )
        END IF


! ********************************
! ***   Evaluate cost function ***
! ********************************

        ! *******************************************
        ! ***   Observation part of cost function ***
        ! *******************************************

        CALL PDAF_timeit(10, 'new')

        CALL PDAF_timeit(31, 'new')

        ! Multiply HV deltav
        CALL gemvTYPE('n', dim_obs_p, dim_ens, 1.0, HV_p, &
             dim_obs_p, v_p, 1, 0.0, HVv_p, 1)

        ! HVv - deltay 
        HVv_p = HVv_p - deltay_p

        ! ***                RiHVv = Rinv HVv                
        ! *** This is implemented as a subroutine thus that
        ! *** Rinv does not need to be allocated explicitly.

        CALL PDAF_timeit(48, 'new')
        CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
        CALL PDAF_timeit(48, 'old')

        ! ***  Compute  J_obs ***

        CALL PDAF_timeit(51, 'new')

        J_obs = 0.0
        DO i = 1, dim_obs_p
           J_obs = J_obs + HVv_p(i)*RiHVv_p(i,1)
        END DO

        J_obs = 0.5*J_obs

        CALL PDAF_timeit(51, 'old')

        CALL PDAF_timeit(31, 'old')


        ! ******************************************
        ! ***   Background part of cost function ***
        ! ******************************************

        J_B = 0.0
        DO i=1, dim_p
           J_B = J_B + v_p(i)*v_p(i)
        END DO
        J_B = 0.5*J_B


        ! *****************************
        ! ***   Total cost function ***
        ! *****************************

        IF (iter>1) J_old = J_tot

        J_tot = J_B + J_obs

        WRITE (*,*) 'J: ', iter, J_B, J_obs, J_tot

        ! Exit condition for steepest descent
        IF (type_min==0 .AND. iter>1 .AND. J_old-J_tot < eps_min) EXIT minloop

        CALL PDAF_timeit(10, 'old')


! **************************
! ***   Compute gradient ***
! **************************

        CALL PDAF_timeit(20, 'new')

        ! Multiply HV deltav
        CALL gemvTYPE('t', dim_obs_p, dim_ens, 1.0, HV_p, &
             dim_obs_p, RiHVv_p, 1, 0.0, gradJ_p, 1)
 
        ! Complete gradient adding v_p
        gradJ_p = v_p + gradJ_p


! *******************************
! ***   Update control vector ***
! *******************************

        IF (type_min == 0) THEN
           ! Steepest descent
           v_p = v_p - 0.1 * gradJ_p
        END IF

        CALL PDAF_timeit(20, 'old')

     END DO minloop


! **************************************************
! ***   Solving completed: Update state estimate ***
! **************************************************

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
          dim_obs_p, v_p, 1, 1.0, state_p, 1)

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
     DEALLOCATE(v_p, HVv_p, RiHVv_p, gradJ_p)
  END IF

  IF (type_min == 1) THEN
     DEALLOCATE(nbd, lvec, uvec, iwa, wa)
  END IF

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_analysis_transf
