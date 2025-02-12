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
!> Optimization loops for ensemble 3D-var
!!
!! There is one subroutine for each solver:
!! * LBFGS: PDAF_hyb3dvar_optim_lbfgs
!! * CGplus: PDAF_hyb3dvar_optim_cgplus
!! * plain CG: PDAF_hyb3dvar_optim_cg
!!
!! In addition routine evaluation the cost function are included:
!! * For LBFGS & CGplus: PDAF_hyb3dvar_costf_cg_cvt
!! * For plain CG: PDAF_hyb3dvar_costf_cg_cvt
!!
MODULE PDAF_hyb3dvar_optim

CONTAINS
!> Optimization with LBFGS for Hyb3dVar
!!
!! Optimization routine for ensemble 3D-Var using the LBFGS solver
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_hyb3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cv_par_p, dim_cv_ens_p, &
     dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, beta_3dvar, screen)

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_3dvar, &
       ONLY: m_lbfgs_var, factr_lbfgs_var, pgtol_lbfgs_var, debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                  !< Current time step
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens               !< ensemble size
  INTEGER, INTENT(in) :: dim_cv_par_p          !< Size of control vector (parameterized)
  INTEGER, INTENT(in) :: dim_cv_ens_p          !< Size of control vector (ensemble)
  INTEGER, INTENT(in) :: dim_obs_p             !< PE-local dimension of observation vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in) :: obs_p(dim_obs_p)         !< Vector of observations
  REAL, INTENT(in) :: dy_p(dim_obs_p)          !< Background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p) !< Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p) !< Control vector (ensemble part)
  INTEGER, INTENT(in) :: screen                !< Verbosity flag
  INTEGER, INTENT(in) :: opt_parallel          !< Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta_3dvar               !< Hybrid weight

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   !< Provide product R^-1 A
       U_cvt, &                !< Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            !< Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            !< Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        !< Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         !< Linearized observation operator
       U_obs_op_adj            !< Adjoint observation operator

! *** local variables ***
  INTEGER :: iter                      ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: dim_cv_p                  ! Full size of control vector
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local full control vector

  ! Variables for LFBGS
  INTEGER            :: m = 5          ! Number of corrections used in limited memory matrix; 3<=m<=20 recommended
  INTEGER            :: iprint
  CHARACTER(len=60)  :: task, csave
  LOGICAL            :: lsave(4)
  INTEGER            :: isave(44)
  REAL               :: factr = 1.0e+7  ! Tolerance in termination test
  REAL               :: pgtol = 1.0e-5  ! Limit for stopping iterations
  REAL               :: dsave(29)
  INTEGER, ALLOCATABLE :: nbd(:), iwa(:)
  REAL, ALLOCATABLE  :: lvec(:), uvec(:), wa(:)


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_LBFGS -- START'

  ! Initialize overall dimension of control vector
  dim_cv_p = dim_cv_par_p + dim_cv_ens_p

  ! Set verbosity of solver
  IF (screen>0 .AND. screen<2) THEN
     iprint = -1
  ELSEIF (screen<3) THEN
     iprint = 0
     IF (mype>0) iprint = -1
  ELSE
     iprint = 99
  END IF
  IF (debug>0) iprint = 99

  ! Allocate arrays
  ALLOCATE(v_p(dim_cv_p))
  ALLOCATE(nbd(dim_cv_p), lvec(dim_cv_p), uvec(dim_cv_p))
  ALLOCATE (iwa(3*dim_cv_p))
  ALLOCATE (wa(2*m*dim_cv_p + 5*dim_cv_p + 11*m*m + 8*m))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 12*dim_cv_p + 2*m*dim_cv_p + 11*m*m + 8*m)

  ! Settings for LBGFS
  m = m_lbfgs_var
  factr = factr_lbfgs_var
  pgtol = pgtol_lbfgs_var
  nbd = 0  ! Values are unbounded
  task = 'START'
  iter = 1

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_LBFGS', debug, &
          'Solver config: m    ', m
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_LBFGS', debug, &
          'Solver config: factr', factr
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_LBFGS', debug, &
          'Solver config: pgtol', pgtol
  END IF

  ! Initialize control vector
  v_p = 0.0


! ***************************
! ***   Iterative solving ***
! ***************************

  ! Prepare arrays for iterations
  ALLOCATE(gradJ_p(dim_cv_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_cv_p)

  IF (mype==0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

  minloop: DO

     IF (.NOT.(task(1:2).EQ.'FG'.OR.task.EQ.'NEW_X'.OR. &
          task.EQ.'START') ) THEN
        IF (mype==0 .AND. screen > 0) &
             WRITE (*,'(a, 5x, a, a)') 'PDAF', '--- Exit optimization, status ', task
        EXIT minloop
     END IF

     ! LBFGS
     CALL PDAF_timeit(54, 'new')
     CALL setulb(dim_cv_p, m, v_p, lvec, uvec, nbd, &
          J_tot, gradJ_p, factr, pgtol, &
          wa, iwa, task, iprint,&
          csave, lsave, isave, dsave )
     CALL PDAF_timeit(54, 'old')


! ********************************
! ***   Evaluate cost function ***
! ********************************

     CALL PDAF_timeit(53, 'new')
     CALL PDAF_hyb3dvar_costf_cvt(step, iter, dim_p, dim_ens, &
          dim_cv_p, dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
          dy_p, v_par_p, v_ens_p, v_p, J_tot, gradJ_p, &
          U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
          U_obs_op_lin, U_obs_op_adj, opt_parallel, beta_3dvar)
     CALL PDAF_timeit(53, 'old')

     IF (mype==0 .AND. screen >2) &
          WRITE (*,'(a, 8x, a, i5, es12.4)') 'PDAF', '--- iter, J: ', iter, J_tot

     iter = iter + 1

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(gradJ_p)
  DEALLOCATE(nbd, lvec, uvec, iwa, wa)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_LBFGS -- END'

END SUBROUTINE PDAF_hyb3dvar_optim_lbfgs

!-------------------------------------------------------------------------------
!> Optimization with CG+ for Hyb3dVar
!!
!! Optimization routine for ensemble 3D-Var using the CG+ solver
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2021-03 - Lars Nerger - Initial code
!! Later revisions - see repository log
!!
SUBROUTINE PDAF_hyb3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cv_par_p, dim_cv_ens_p, &
     dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, beta_3dvar, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, comm_filter, npes_filter, MPIerr
  USE PDAF_3dvar, &
       ONLY: method_cgplus_var, irest_cgplus_var, eps_cgplus_var, debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                  !< Current time step
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens               !< ensemble size
  INTEGER, INTENT(in) :: dim_cv_par_p          !< Size of control vector (parameterized)
  INTEGER, INTENT(in) :: dim_cv_ens_p          !< Size of control vector (ensemble)
  INTEGER, INTENT(in) :: dim_obs_p             !< PE-local dimension of observation vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)        !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)         !< Background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p) !< Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p) !< Control vector (ensemble part)
  INTEGER, INTENT(in) :: screen                !< Verbosity flag
  INTEGER, INTENT(in) :: opt_parallel          !< Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta_3dvar               !< Hybrid weight

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   !< Provide product R^-1 A
       U_cvt, &                !< Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            !< Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            !< Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        !< Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         !< Linearized observation operator
       U_obs_op_adj            !< Adjoint observation operator

! *** local variables ***
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: optiter                   ! Additional iteration counter
  INTEGER :: dim_cv_p                  ! Process-local full size of control vector
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local full control vector

  ! Variables for CG+
  INTEGER :: method=2
  INTEGER :: irest=5
  REAL :: eps=1.0e-5
  INTEGER :: iprint(2), iflag, icall, mp, lp, i
  REAL, ALLOCATABLE :: d(:), gradJ_old_p(:), w(:)
  REAL :: tlev
  LOGICAL :: finish, update_J
  INTEGER :: ifinish_p, iupdate_J_p    ! Flags used for MPI_allreduce to determine exit status
  INTEGER :: ifinish, iupdate_J        ! Flags used for MPI_allreduce to determine exit status
  INTEGER :: iter, nfun
  COMMON /cgdd/    mp,lp
  COMMON /runinf/  iter,nfun


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CGPLUS -- START'

  ! Initialize overall dimension of control vector
  dim_cv_p = dim_cv_par_p + dim_cv_ens_p

  ! Settings for CG+
  method = method_cgplus_var  ! (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere (default=2)
  irest = irest_cgplus_var    ! (0) no restarts; (1) restart every n steps (default=1)
  EPS = eps_cgplus_var        ! Convergence constant (default=1.0e-5)
  icall = 0
  iflag = 0
  FINISH = .FALSE.
  update_J = .TRUE.
  ifinish_p = 0
  iupdate_J_p = 0
  optiter = 1

  ! Set verbosity of solver
  IF (screen>0 .AND. screen<2) THEN
     iprint(1) = -1
  ELSEIF (screen==2) THEN
     iprint(1) = 0
     IF (mype>0) iprint(1) = -1
  ELSE
     iprint(1) = 0
  END IF
  iprint(2) = 0  
  if (debug>0) iprint(1) = 0

  ! Allocate arrays
  ALLOCATE(v_p(dim_cv_p))
  ALLOCATE(d(dim_cv_p), w(dim_cv_p))
  ALLOCATE(gradJ_p(dim_cv_p), gradJ_old_p(dim_cv_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4*dim_cv_p)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: method  ', method
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: restarts', irest
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: EPS     ', EPS
  END IF

  ! Initialize numbers
  v_p = 0.0


! ***************************
! ***   Iterative solving ***
! ***************************

  IF (mype==0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

  minloop: DO

! ********************************
! ***   Evaluate cost function ***
! ********************************

     IF (update_J) THEN
        CALL PDAF_timeit(53, 'new')
        CALL PDAF_hyb3dvar_costf_cvt(step, optiter, dim_p, dim_ens, &
             dim_cv_p, dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
             dy_p, v_par_p, v_ens_p, v_p, J_tot, gradJ_p, &
             U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
             U_obs_op_lin, U_obs_op_adj, opt_parallel, beta_3dvar)
        CALL PDAF_timeit(53, 'old')
     END IF


! ***************************
! ***   Optimize with CG+ ***
! ***************************

     CALL PDAF_timeit(54, 'new')
     IF (opt_parallel==0) THEN
        CALL CGFAM(dim_cv_p, v_p, J_tot, gradJ_p, d, gradJ_old_p, IPRINT, EPS, W,  &
             iflag, IREST, METHOD, FINISH)
     ELSE
        CALL CGFAM_mpi(dim_cv_p, v_p, J_tot, gradJ_p, d, gradJ_old_p, IPRINT, EPS, W,  &
             iflag, IREST, METHOD, FINISH, comm_filter, npes_filter)
     END IF
     CALL PDAF_timeit(54, 'old')


     ! *** Check exit status ***

     ! iflag=
     !    0 : successful termination
     !    1 : return to evaluate F and G
     !    2 : return with a new iterate, try termination test
     !   -i : error

     update_J = .TRUE.
     IF (iflag <= 0 .OR. icall > 10000) EXIT minloop
     IF (iflag == 1 ) icall = icall + 1
     IF (iflag == 2) THEN
        ! Termination Test.
        tlev = eps*(1.0 + ABS(J_tot))
        i=0

        ! Process-local check
        checktest: DO
           i = i + 1
           IF(i > dim_cv_p) THEN
              FINISH = .TRUE.
              update_J = .FALSE.
              ifinish_p = 1
              iupdate_J_p = 1
              EXIT checktest
           ENDIF
           IF(ABS(gradJ_p(i)) > tlev) THEN
              update_J = .FALSE.
              iupdate_J_p = 1
              EXIT checktest
           ENDIF
        END DO checktest

        ! Global check
        IF (npes_filter > 1) THEN
           CALL MPI_ALLREDUCE(ifinish_p, ifinish, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
           CALL MPI_ALLREDUCE(iupdate_J_p, iupdate_J, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
           IF (ifinish>0) THEN
              FINISH = .TRUE.
              update_J = .FALSE.
           END IF
           IF (iupdate_J>0) THEN
              update_J = .FALSE.
           END IF
        END IF
     ENDIF

     ! Increment loop counter
     optiter = optiter+1

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  ! Initialize two partial control vectors
  v_par_p = v_p(1 : dim_cv_par_p)
  v_ens_p = v_p(dim_cv_par_p+1 : dim_cv_p)

  DEALLOCATE(gradJ_p, v_p)
  DEALLOCATE(d, gradJ_old_p, w)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CGPLUS -- END'

END SUBROUTINE PDAF_hyb3dvar_optim_cgplus

!-------------------------------------------------------------------------------
!> Optimization with parallelized CG for Hyb3dVar
!!
!! Optimization routine for ensemble 3D-Var using direct 
!! parallelized implementation of CG.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_hyb3dvar_optim_cg(step, dim_p, dim_ens, dim_cv_par_p, dim_cv_ens_p, &
     dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, beta_3dvar, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, Comm_filter, MPI_REALTYPE, MPI_SUM, MPIerr
  USE PDAF_3dvar, &
       ONLY: maxiter_cg_var, eps_cg_var, debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                  !< Current time step
  INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens               !< ensemble size
  INTEGER, INTENT(in) :: dim_cv_par_p          !< Size of control vector (parameterized)
  INTEGER, INTENT(in) :: dim_cv_ens_p          !< Size of control vector (ensemble)
  INTEGER, INTENT(in) :: dim_obs_p             !< PE-local dimension of observation vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)        !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)         !< Background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p) !< Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p) !< Control vector (ensemble part)
  INTEGER, INTENT(in) :: screen                !< Verbosity flag
  INTEGER, INTENT(in) :: opt_parallel          !< Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta_3dvar               !< Hybrid weight

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   !< Provide product R^-1 A
       U_cvt, &                !< Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            !< Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            !< Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        !< Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         !< Linearized observation operator
       U_obs_op_adj            !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i, iter                       ! Iteration counter
  INTEGER :: maxiter=200                   ! maximum number of iterations
  INTEGER, SAVE :: allocflag = 0           ! Flag whether first time allocation is done
  REAL :: J_tot, J_old                     ! Cost function
  REAL, ALLOCATABLE :: gradJ_par_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: gradJ_ens_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: hessJd_par_p(:)     ! Hessian times v
  REAL, ALLOCATABLE :: hessJd_ens_p(:)     ! Hessian times v
  REAL :: gprod_p, dprod_p, gprod_new_p    ! temporary variables for step size computation
  REAL :: gprod, dprod, gprod_new          ! temporary variables for step size computation
  REAL :: alpha, beta                      ! step sizes
  REAL, ALLOCATABLE :: d_par_p(:)          ! descent direction
  REAL, ALLOCATABLE :: d_ens_p(:)          ! descent direction
  REAL, ALLOCATABLE :: v2_par_p(:)         ! iterated control vector
  REAL, ALLOCATABLE :: v2_ens_p(:)         ! iterated control vector
  REAL, ALLOCATABLE :: gradJ2_par_p(:)     ! iterated gradient
  REAL, ALLOCATABLE :: gradJ2_ens_p(:)     ! iterated gradient
  REAL, ALLOCATABLE :: d2_par_p(:)         ! iterated descent direction
  REAL, ALLOCATABLE :: d2_ens_p(:)         ! iterated descent direction
  REAL :: eps=1.0e-6                       ! Convergence condition value


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CG -- START'

  maxiter = maxiter_cg_var    ! Maximum number of iterations (default=200)
  eps = eps_cg_var            ! Convergence limit (default=1.0e-6)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CG', debug, &
          'Solver config: maxiter ', maxiter
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CG', debug, &
          'Solver config: EPS     ', EPS
  END IF

  ! Prepare arrays for iterations
  ALLOCATE(gradJ_par_p(dim_cv_par_p))
  ALLOCATE(gradJ_ens_p(dim_cv_ens_p))
  ALLOCATE(hessJd_par_p(dim_cv_par_p))
  ALLOCATE(hessJd_ens_p(dim_cv_ens_p))
  ALLOCATE(d_par_p(dim_cv_par_p))
  ALLOCATE(d_ens_p(dim_cv_ens_p))
  ALLOCATE(v2_par_p(dim_cv_ens_p), gradJ2_par_p(dim_cv_ens_p), d2_par_p(dim_cv_ens_p))
  ALLOCATE(v2_ens_p(dim_cv_ens_p), gradJ2_ens_p(dim_cv_ens_p), d2_ens_p(dim_cv_ens_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 6*dim_cv_ens_p)

  ! Initialize numbers
  J_tot = 0.0
  d_par_p = 0.0
  d_ens_p = 0.0
  

! ***************************
! ***   Iterative solving ***
! ***************************

  IF (mype==0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

  minloop: DO iter = 1, maxiter

! ********************************
! ***   Evaluate cost function ***
! ********************************

     CALL PDAF_timeit(53, 'new')
     J_old = J_tot
     CALL PDAF_hyb3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens, &
          dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
          dy_p, v_par_p, v_ens_p, d_par_p, d_ens_p, &
          J_tot, gradJ_par_p, gradJ_ens_p, hessJd_par_p, hessJd_ens_p, &
          U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
          U_obs_op_lin, U_obs_op_adj, opt_parallel, beta_3dvar)
     CALL PDAF_timeit(53, 'old')

     IF (mype==0 .AND. screen > 2) &
          WRITE (*,'(a, 8x, a, i5, 1x, es14.6)') 'PDAF', '--- iter, J: ', iter, J_tot
     IF (iter>1 .AND. J_old - J_tot < eps) THEN
        IF (mype==0 .AND. screen >= 2) THEN
           WRITE (*,'(a, 8x, a)') 'PDAF', '--- CG solver converged'
           WRITE (*,'(a, 12x, a6, 4x, a4, 7x, a4)') 'PDAF', 'iter', 'eps','F'
           WRITE (*,'(a, 13x, i4, 2x, es10.3, 2x, es10.3/)') 'PDAF', iter, J_old-J_tot, J_tot
        END IF
        EXIT minloop
     END IF


! ***************************
! ***   Optimize with CG  ***
! ***************************

     CALL PDAF_timeit(54, 'new')

     ! Compute step size alpha
     IF (iter==1) THEN
        gprod_p = 0.0
        DO i=1, dim_cv_par_p
           gprod_p = gprod_p + gradJ_par_p(i)*gradJ_par_p(i)
        END DO
        DO i=1, dim_cv_ens_p
           gprod_p = gprod_p + gradJ_ens_p(i)*gradJ_ens_p(i)
        END DO
     
        IF (opt_parallel==1) THEN
           ! Get global value
           CALL MPI_Allreduce(gprod_p, gprod, 1, MPI_REALTYPE, MPI_SUM, &
                COMM_filter, MPIerr)
        ELSE
           gprod = gprod_p
        END IF
     END IF

     dprod_p = 0.0
     DO i=1, dim_cv_ens_p
        dprod_p = dprod_p + d_par_p(i)*hessJd_par_p(i)
     END DO
     DO i=1, dim_cv_ens_p
        dprod_p = dprod_p + d_ens_p(i)*hessJd_ens_p(i)
     END DO
     
     IF (opt_parallel==1) THEN
        ! Get global value
        CALL MPI_Allreduce(dprod_p, dprod, 1, MPI_REALTYPE, MPI_SUM, &
             COMM_filter, MPIerr)
     ELSE
        dprod = dprod_p
     END IF

     alpha = gprod / dprod

     ! Update control vector
     v2_par_p = v_par_p + alpha * d_par_p
     v2_ens_p = v_ens_p + alpha * d_ens_p

     ! Update gradient
     gradJ2_par_p = gradJ_par_p + alpha * hessJd_par_p
     gradJ2_ens_p = gradJ_ens_p + alpha * hessJd_ens_p

     ! Compute step size beta for update of descent direction
     gprod_new_p = 0.0
     DO i=1, dim_cv_par_p
        gprod_new_p = gprod_new_p + gradJ2_par_p(i)*gradJ2_par_p(i)
     END DO
     DO i=1, dim_cv_ens_p
        gprod_new_p = gprod_new_p + gradJ2_ens_p(i)*gradJ2_ens_p(i)
     END DO
     
     IF (opt_parallel==1) THEN
        ! Get global value
        CALL MPI_Allreduce(gprod_new_p, gprod_new, 1, MPI_REALTYPE, MPI_SUM, &
             COMM_filter, MPIerr)
     ELSE
        gprod_new = gprod_new_p
     END IF

     beta = gprod_new / gprod

     ! Update descent direction
     d2_par_p = - gradJ2_par_p + beta * d_par_p
     d2_ens_p = - gradJ2_ens_p + beta * d_ens_p
     
     ! prepare next iteration
     gradJ_par_p = gradJ2_par_p
     gradJ_ens_p = gradJ2_ens_p
     d_par_p = d2_par_p
     d_ens_p = d2_ens_p
     v_par_p = v2_par_p
     v_ens_p = v2_ens_p
     gprod = gprod_new

     CALL PDAF_timeit(54, 'old')

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(gradJ_ens_p)
  DEALLOCATE(gradJ_par_p)
  DEALLOCATE(hessJd_ens_p, d_ens_p, v2_ens_p, gradJ2_ens_p, d2_ens_p)
  DEALLOCATE(hessJd_par_p, d_par_p, v2_par_p, gradJ2_par_p, d2_par_p)
  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CG -- END'

END SUBROUTINE PDAF_hyb3dvar_optim_cg

!-------------------------------------------------------------------------------
!> Evaluate cost function and its gradient
!!
!! Routine to evaluate the cost function and its gradient
!! for the incremental hybrid 3D-Var with variable transformation.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_hyb3dvar_costf_cvt(step, iter, dim_p, dim_ens, &
     dim_cv_p, dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
     dy_p, v_par_p, v_ens_p, v_p, J_tot, gradJ, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
     U_obs_op_lin, U_obs_op_adj, opt_parallel, beta)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE
  USE PDAF_mod_filter, &
       ONLY: debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                   !< Current time step
  INTEGER, INTENT(in) :: iter                   !< Optimization iteration
  INTEGER, INTENT(in) :: dim_p                  !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                !< ensemble size
  INTEGER, INTENT(in) :: dim_cv_p               !< Size of control vector (full)
  INTEGER, INTENT(in) :: dim_cv_par_p           !< Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cv_ens_p           !< Size of control vector (ensemble part)
  INTEGER, INTENT(in) :: dim_obs_p              !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)         !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)          !< background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p)  !< Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p)  !< Control vector (ensemble part)
  REAL, INTENT(in)  :: v_p(dim_cv_p)            !< Control vector (full)
  REAL, INTENT(out) :: J_tot                    !< on exit: Value of cost function
  REAL, INTENT(out) :: gradJ(dim_cv_p)          !< on exit: PE-local gradient of J (full)
  INTEGER, INTENT(in) :: opt_parallel           !< Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta                      !< Hybrid weight

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   !< Provide product R^-1 A
       U_cvt, &                !< Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            !< Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            !< Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        !< Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         !< Linearized observation operator
       U_obs_op_adj            !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: Vv_p(:)         ! PE-local product V deltav (parameterized and total)
  REAL, ALLOCATABLE :: Vv_ens_p(:)     ! PE-local product V deltav (ensemble)
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local product HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: gradJ_par(:)    ! PE-local part of gradJ (parameterized)
  REAL, ALLOCATABLE :: gradJ_ens(:)    ! PE-local part of gradJ (ensemble)
  REAL :: J_B_p, J_B, J_obs_p, J_obs   ! Cost function terms
  REAL :: sbeta, sombeta               ! square-root of beta and one minus beta


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- START: iteration', iter

  ! Initialize parts of control vector
  v_par_p = v_p(1 : dim_cv_par_p)
  v_ens_p = v_p(dim_cv_par_p+1 : dim_cv_p)

  ! Allocate arrays
  ALLOCATE(Vv_p(dim_p))
  ALLOCATE(Vv_ens_p(dim_p))
  ALLOCATE(HVv_p(dim_obs_p))
  ALLOCATE(RiHVv_p(dim_obs_p, 1))
  ALLOCATE(gradJ_par(dim_cv_par_p))
  ALLOCATE(gradJ_ens(dim_cv_ens_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p + dim_cv_ens_p + dim_cv_par_p + 2*dim_p)

  ! Initialize numbers
  sbeta = SQRT(beta)
  sombeta = SQRT(1.0 - beta)


! *******************************************
! ***   Observation part of cost function ***
! *******************************************

  CALL PDAF_timeit(55, 'new')

  CALL PDAF_timeit(56, 'new')

  ! *** Apply V to control vector v_p ***

  Vv_p = 0.0
  Vv_ens_p = 0.0

  ! parameterized
  IF (dim_cv_par_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt'

     CALL PDAF_timeit(60, 'new')
     CALL U_cvt(iter, dim_p, dim_cv_par_p, v_par_p, Vv_p)
     CALL PDAF_timeit(60, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'state increment after CVT dX_par(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dX_par', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF
  END IF

  ! ensemble
  IF (dim_cv_ens_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_ens'

     CALL PDAF_timeit(61, 'new')
     CALL U_cvt_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, v_ens_p, Vv_ens_p)
     CALL PDAF_timeit(61, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'state increment after CVT dX_ens(1:min(dim_p,6))', Vv_ens_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dX_ens', MINVAL(Vv_ens_p), MAXVAL(Vv_ens_p)
     END IF
  END IF

  Vv_p = sombeta*Vv_p + sbeta*Vv_ens_p

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'hybrid sum state increment after CVT dX(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX dX', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply linearized observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'observed state after CVT (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX HdX', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! HVv - dy 
  CALL PDAF_timeit(51, 'new')
  HVv_p = HVv_p - dy_p
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'process local residual d (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX d', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVv = Rinv HVv                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'R^-1 d (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX R^-1 d', MINVAL(RiHVv_p), MAXVAL(RiHVv_p)
  END IF

  ! ***  Compute  J_obs ***

  CALL PDAF_timeit(51, 'new')

  J_obs_p = 0.0
  DO i = 1, dim_obs_p
     J_obs_p = J_obs_p + HVv_p(i)*RiHVv_p(i,1)
  END DO

  J_obs_p = 0.5*J_obs_p

  IF (debug>0) &
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'process local observation cost J_obs', J_obs_p

  ! Get global value
  CALL MPI_Allreduce(J_obs_p, J_obs, 1, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

  IF (debug>0) &
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'global observation cost J_obs', J_obs

  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(56, 'old')


! ******************************************
! ***   Background part of cost function ***
! ******************************************

  CALL PDAF_timeit(57, 'new')
  CALL PDAF_timeit(51, 'new')

  J_B_p = 0.0
  DO i = 1, dim_cv_par_p
     J_B_p = J_B_p + v_par_p(i)*v_par_p(i)
  END DO
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'process local CV cost parameterized part J_B_par', 0.5 * J_B_p

  DO i = 1, dim_cv_ens_p
     J_B_p = J_B_p + v_ens_p(i)*v_ens_p(i)
  END DO
  J_B_p = 0.5 * J_B_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'process local CV cost par+ens J_B', J_B_p

  IF (opt_parallel==1) THEN
     ! Get global value
     CALL MPI_Allreduce(J_B_p, J_B, 1, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     J_B = J_B_p
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'global CV cost J_B', J_B

  CALL PDAF_timeit(57, 'old')


! *****************************
! ***   Total cost function ***
! *****************************

  J_tot = J_B + J_obs

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(55, 'old')


! **************************
! ***   Compute gradient ***
! **************************

  CALL PDAF_timeit(58, 'new')

  ! Apply adjoint of observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_adj'

  CALL PDAF_timeit(65, 'new')
  Vv_p = 0.0
  CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
  CALL PDAF_timeit(65, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'H^TR^-1 d (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX H^TR^-1 d', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply V^T to vector
  IF (dim_cv_par_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj'

     CALL PDAF_timeit(62, 'new')
     CALL U_cvt_adj(iter, dim_p, dim_cv_par_p, Vv_p, gradJ_par)
     CALL PDAF_timeit(62, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'CVT_par(H^TR^-1 d) (1:min(dim_cv_par_p,6))', gradJ_par(1:min(dim_cv_par_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX CVT_par(H^TR^-1 d)', MINVAL(gradJ_par), MAXVAL(gradJ_par)
     END IF
  END IF

  IF (dim_cv_ens_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj_ens'

     CALL PDAF_timeit(63, 'new')
     CALL U_cvt_adj_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, gradJ_ens)
     CALL PDAF_timeit(63, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'CVT_ens(H^TR^-1 d) (1:min(dim_cv_ens_p,6))', gradJ_ens(1:min(dim_cv_ens_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX CVT_ens(H^TR^-1 d)', MINVAL(gradJ_ens), MAXVAL(gradJ_ens)
     END IF
  END IF

  ! Complete gradient adding v_p
  CALL PDAF_timeit(51, 'new')
  DO i = 1, dim_cv_par_p
     gradJ(i) = v_par_p(i) + sombeta*gradJ_par(i)
  END DO
  DO i = 1, dim_cv_ens_p
     gradJ(i + dim_cv_par_p) = v_ens_p(i) + sbeta*gradJ_ens(i)
  END DO
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'process local gradient gradJ (1:min(dim_cv_p,6))', gradJ(1:min(dim_cv_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX gradJ', MINVAL(gradJ), MAXVAL(gradJ)
  END IF

  CALL PDAF_timeit(58, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Vv_p, HVv_p, RiHVv_p, gradJ_par, gradJ_ens)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- END'

END SUBROUTINE PDAF_hyb3dvar_costf_cvt

!-------------------------------------------------------------------------------
!> Evaluate cost function, its gradient and Hessian
!!
!! Routine to evaluate the cost function, its gradient, and
!! the product of its Hessian time descent direction
!! for the incremental hybrid 3D-Var with variable transformation.
!!
!! The subroutine distinguishes two cases:
!! iter==1
!!   In this case all quantities are computed, the 
!!   descent direction is initialized from the gradient vector
!! iter>1
!!   In this case only the cost function value and the product
!!   of the Hessian times descent direction are computed.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2021-03 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
SUBROUTINE PDAF_hyb3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens, &
     dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
     dy_p, v_par_p, v_ens_p, d_par_p, d_ens_p, &
     J_tot, gradJ_par, gradJ_ens, hessJd_par, hessJd_ens, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
     U_obs_op_lin, U_obs_op_adj, opt_parallel, beta)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE
  USE PDAF_mod_filter, &
       ONLY: debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                   !< Current time step
  INTEGER, INTENT(in) :: iter                   !< Optimization iteration
  INTEGER, INTENT(in) :: dim_ens                !< ensemble size
  INTEGER, INTENT(in) :: dim_p                  !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cv_par_p           !< Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cv_ens_p           !< Size of control vector (ensemble part)
  INTEGER, INTENT(in) :: dim_obs_p              !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)         !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)          !< Background innovation
  REAL, INTENT(in)  :: v_par_p(dim_cv_par_p)    !< Control vector (parameterized part)
  REAL, INTENT(in)  :: v_ens_p(dim_cv_ens_p)    !< Control vector (ensemble part)
  REAL, INTENT(inout) :: d_par_p(dim_cv_par_p)  !< CG descent direction (parameterized part)
  REAL, INTENT(inout) :: d_ens_p(dim_cv_ens_p)  !< CG descent direction (ensemble part)
  REAL, INTENT(out) :: J_tot                    !< on exit: Value of cost function
  REAL, INTENT(out) :: gradJ_par(dim_cv_par_p)  !< on exit: gradient of J (parameterized part)
  REAL, INTENT(out) :: gradJ_ens(dim_cv_ens_p)  !< on exit: gradient of J (ensemble part)
  REAL, INTENT(out) :: hessJd_par(dim_cv_par_p) !< on exit: Hessian of J times d_p (parameterized part)
  REAL, INTENT(out) :: hessJd_ens(dim_cv_ens_p) !< on exit: Hessian of J times d_p (ensemble part)
  INTEGER, INTENT(in) :: opt_parallel           !< Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta                      !< Hybrid weight

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   !< Provide product R^-1 A
       U_cvt, &                !< Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            !< Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            !< Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        !< Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         !< Linearized observation operator
       U_obs_op_adj            !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: Vv_p(:)         ! PE-local product V deltav (parameterized and total)
  REAL, ALLOCATABLE :: Vv_ens_p(:)     ! PE-local product V deltav (ensemble)
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local produce HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL :: J_B_p, J_B, J_obs_p, J_obs   ! Cost function terms
  REAL :: sbeta, sombeta               ! square-root of beta and one minus beta


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- START: iteration', iter

  ! Allocate arrays
  ALLOCATE(Vv_p(dim_p))
  ALLOCATE(Vv_ens_p(dim_p))
  ALLOCATE(HVv_p(dim_obs_p))
  ALLOCATE(RiHVv_p(dim_obs_p, 1))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p + 2*dim_p)

  ! Initialize numbers
  sbeta = SQRT(beta)
  sombeta = SQRT(1.0 - beta)


! *******************************************
! ***   Observation part of cost function ***
! *******************************************

  CALL PDAF_timeit(55, 'new')

  CALL PDAF_timeit(56, 'new')

  ! *** Apply V to control vector v_p ***

  Vv_p = 0.0
  Vv_ens_p = 0.0

  ! parameterized
  IF (dim_cv_par_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt'

     CALL PDAF_timeit(60, 'new')
     CALL U_cvt(iter, dim_p, dim_cv_par_p, v_par_p, Vv_p)
     CALL PDAF_timeit(60, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'state increment after CVT dX_par(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dX_par', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF
  END IF

  ! ensemble
  IF (dim_cv_ens_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_ens'

     CALL PDAF_timeit(61, 'new')
     CALL U_cvt_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, v_ens_p, Vv_ens_p)
     CALL PDAF_timeit(61, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'state increment after CVT dX_ens(1:min(dim_p,6))', Vv_ens_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dX_ens', MINVAL(Vv_ens_p), MAXVAL(Vv_ens_p)
     END IF
  END IF

  Vv_p = sombeta*Vv_p + sbeta*Vv_ens_p

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'hybrid sum state increment after CVT dX(1:min(dim_p,6))', Vv_ens_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX dX', MINVAL(Vv_ens_p), MAXVAL(Vv_ens_p)
  END IF

  ! Apply linearized observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'observed state after CVT (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX HdX', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! HVv - dy
  CALL PDAF_timeit(51, 'new')
  HVv_p = HVv_p - dy_p
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'process local residual d (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX d', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVv = Rinv HVv                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'R^-1 d (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX R^-1 d', MINVAL(RiHVv_p), MAXVAL(RiHVv_p)
  END IF

  ! ***  Compute  J_obs ***

  CALL PDAF_timeit(51, 'new')

  J_obs_p = 0.0
  DO i = 1, dim_obs_p
     J_obs_p = J_obs_p + HVv_p(i)*RiHVv_p(i,1)
  END DO

  J_obs_p = 0.5*J_obs_p

  IF (debug>0) &
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'process local observation cost J_obs', J_obs_p

  ! Get global value
  CALL MPI_Allreduce(J_obs_p, J_obs, 1, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

  IF (debug>0) &
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'global observation cost J_obs', J_obs

  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(56, 'old')


! ******************************************
! ***   Background part of cost function ***
! ******************************************

  CALL PDAF_timeit(57, 'new')
  CALL PDAF_timeit(51, 'new')

  J_B_p = 0.0
  DO i = 1, dim_cv_par_p
     J_B_p = J_B_p + v_par_p(i)*v_par_p(i)
  END DO
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'process local CV cost parameterized part J_B_par', 0.5 * J_B_p

  DO i = 1, dim_cv_ens_p
     J_B_p = J_B_p + v_ens_p(i)*v_ens_p(i)
  END DO
  J_B_p = 0.5 * J_B_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'process local CV cost par+ens J_B', J_B_p

  IF (opt_parallel==1) THEN
     ! Get global value
     CALL MPI_Allreduce(J_B_p, J_B, 1, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     J_B = J_B_p
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
       'global CV cost J_B', J_B

  CALL PDAF_timeit(57, 'old')


! *****************************
! ***   Total cost function ***
! *****************************

  J_tot = J_B + J_obs

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(55, 'old')


! **************************
! ***   Compute gradient ***
! **************************

  ! Only at first iteration
  IF (iter==1) THEN

     CALL PDAF_timeit(58, 'new')

     ! Apply adjoint of observation operator
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_adj'

     CALL PDAF_timeit(65, 'new')
     Vv_p = 0.0
     CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
     CALL PDAF_timeit(65, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'H^TR^-1 d (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX H^TR^-1 d', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF

     ! Apply V^T to vector
     IF (dim_cv_par_p>0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj'

        CALL PDAF_timeit(62, 'new')
        CALL U_cvt_adj(iter, dim_p, dim_cv_par_p, Vv_p, gradJ_par)
        CALL PDAF_timeit(62, 'old')

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'CVT_par(H^TR^-1 d) (1:min(dim_cv_par_p,6))', gradJ_par(1:min(dim_cv_par_p,6))
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'MIN/MAX CVT_par(H^TR^-1 d)', MINVAL(gradJ_par), MAXVAL(gradJ_par)
        END IF
     END IF
     IF (dim_cv_ens_p>0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj_ens'

        CALL PDAF_timeit(63, 'new')
        CALL U_cvt_adj_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, gradJ_ens)
        CALL PDAF_timeit(63, 'old')

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'CVT_ens(H^TR^-1 d) (1:min(dim_cv_ens_p,6))', gradJ_ens(1:min(dim_cv_ens_p,6))
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'MIN/MAX CVT_ens(H^TR^-1 d)', MINVAL(gradJ_ens), MAXVAL(gradJ_ens)
        END IF
     END IF

     ! Complete gradient adding v_p
     CALL PDAF_timeit(51, 'new')
     IF (dim_cv_par_p>0) THEN
        gradJ_par = v_par_p + sombeta*gradJ_par

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'process local gradient gradJ_par (1:min(dim_cv_p,6))', gradJ_par(1:min(dim_cv_par_p,6))
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'MIN/MAX gradJ_par', MINVAL(gradJ_par), MAXVAL(gradJ_par)
        END IF
     END IF
     IF (dim_cv_ens_p>0) THEN
        gradJ_ens = v_ens_p + sbeta*gradJ_ens

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'process local gradient gradJ_ens (1:min(dim_cv_p,6))', gradJ_ens(1:min(dim_cv_ens_p,6))
           WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
                'MIN/MAX gradJ_ens', MINVAL(gradJ_ens), MAXVAL(gradJ_ens)
        END IF
     END IF
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(58, 'old')

  END IF


! *****************************************************
! ***   Compute Hessian times direction vector d_p  ***
! *****************************************************

  CALL PDAF_timeit(59, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- compute Hessian times direction cvp'

  ! Initialize descent direction d_p at first iteration
  IF (iter==1) THEN
     CALL PDAF_timeit(51, 'new')
       IF (dim_cv_par_p>0) d_par_p = - gradJ_par
       IF (dim_cv_ens_p>0) d_ens_p = - gradJ_ens
     CALL PDAF_timeit(51, 'old')
  END IF

  ! Apply V to control vector v_p
  Vv_p = 0.0
  Vv_ens_p = 0.0
  IF (dim_cv_par_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt'

     CALL PDAF_timeit(60, 'new')
     CALL U_cvt(-iter, dim_p, dim_cv_par_p, d_par_p, Vv_p)
     CALL PDAF_timeit(60, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'direction in state space after CVT dp_par(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dp_par', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF
  END IF
  IF (dim_cv_ens_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_ens'

     CALL PDAF_timeit(61, 'new')
     CALL U_cvt_ens(-iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, d_ens_p, Vv_ens_p)
     CALL PDAF_timeit(61, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'direction in state space after CVT dp_ens(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX dp_ens', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF
  END IF

  Vv_p = sombeta*Vv_p + sbeta*Vv_ens_p

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'hybrid sum state increment after CVT dp(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX dp', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'observed direction Hdp (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX Hdp', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVd = Rinv HVd                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  ! *** RiHVd is stored in RiHVv
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'R^-1 Hdp (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX R^-1 Hdp', MINVAL(RiHVv_p), MAXVAL(RiHVv_p)
  END IF

  ! Apply adjoint of observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call obs_op_adj'

  CALL PDAF_timeit(65, 'new')
  Vv_p = 0.0
  CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
  CALL PDAF_timeit(65, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'H^TR^-1 dp (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
          'MIN/MAX H^TR^-1 dp', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply V^T to vector
  IF (dim_cv_par_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj'

     CALL PDAF_timeit(62, 'new')
     CALL U_cvt_adj(-iter, dim_p, dim_cv_par_p, Vv_p, hessJd_par)
     CALL PDAF_timeit(62, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'CVT_par(H^TR^-1 dp) (1:min(dim_cv_p,6))', hessJd_par(1:min(dim_cv_par_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX CVT_par(H^TR^-1 dp)', MINVAL(hessJd_par), MAXVAL(hessJd_par)
     END IF
  END IF
  IF (dim_cv_ens_p>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- call cvt_adj_ens'

     CALL PDAF_timeit(63, 'new')
     CALL U_cvt_adj_ens(-iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, hessJd_ens)
     CALL PDAF_timeit(63, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'CVT_ens(H^TR^-1 dp) (1:min(dim_cv_p,6))', hessJd_ens(1:min(dim_cv_ens_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX CVT_ens(H^TR^-1 dp)', MINVAL(hessJd_ens), MAXVAL(hessJd_ens)
     END IF
  END IF

  ! Add d_p to complete Hessian times d_p
  CALL PDAF_timeit(51, 'new')
  IF (dim_cv_par_p>0) THEN
     hessJd_par = sombeta*hessJd_par + d_par_p

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'hybrid sum (Hessian dp)_par (1:min(dim_cv_p,6))', hessJd_par(1:min(dim_cv_par_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX (Hessian dp)_par', MINVAL(hessJd_par), MAXVAL(hessJd_par)
     END IF
  END IF
  IF (dim_cv_ens_p>0) THEN
     hessJd_ens = sbeta*hessJd_ens + d_ens_p

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'hybrid sum (Hessian dp)_ens (1:min(dim_cv_p,6))', hessJd_ens(1:min(dim_cv_par_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_costf_cvt:', debug, &
             'MIN/MAX (Hessian dp)_ens', MINVAL(hessJd_ens), MAXVAL(hessJd_ens)
     END IF
  END IF
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(59, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Vv_p, Vv_ens_p, HVv_p, RiHVv_p)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_costf_cvt -- END'

END SUBROUTINE PDAF_hyb3dvar_costf_cg_cvt

END MODULE PDAF_hyb3dvar_optim
