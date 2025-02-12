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
!
!> Optimization loops for 3D-var
!!
!! There is one subroutine for each solver:
!! * LBFGS: PDAF_3dvar_optim_lbfgs
!! * CGplus: PDAF_3dvar_optim_cgplus
!! * plain CG: PDAF_3dvar_optim_cg
!!
!! In addition routine evaluation the cost function are included:
!! * For LBFGS & CGplus: PDAF_3dvar_costf_cg_cvt
!! * For plain CG: PDAF_3dvar_costf_cg_cvt
!!
MODULE PDAF_3dvar_optim

CONTAINS
!> Optimization loop for LBFGS
!!
!! Optimization routine for 3D-Var using the LBFGS solver
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
SUBROUTINE PDAF_3dvar_optim_lbfgs(step, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, screen)

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
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_p             !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p        !< Size of control vector
  INTEGER, INTENT(in) :: dim_obs_p         !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)    !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)     !< Background innovation
  REAL, INTENT(inout) :: v_p(dim_cvec_p)   !< Control vector
  INTEGER, INTENT(in) :: opt_parallel      !< Whether to use a decomposed control vector
  INTEGER, INTENT(in) :: screen            !< Verbosity flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &               !< Provide product R^-1 A
       U_cvt, &                            !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                        !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &                     !< Linearized observation operator
       U_obs_op_adj                        !< Adjoint observation operator

! *** local variables ***
  INTEGER :: iter                      ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J

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
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_LBFGS -- START'

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
  ALLOCATE(nbd(dim_cvec_p), lvec(dim_cvec_p), uvec(dim_cvec_p))
  ALLOCATE (iwa(3*dim_cvec_p))
  ALLOCATE (wa(2*m*dim_cvec_p + 5*dim_cvec_p + 11*m*m + 8*m))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 11*dim_cvec_p + 2*m*dim_cvec_p + 11*m*m + 8*m)

  ! Settings for LBGFS
  m = m_lbfgs_var
  factr = factr_lbfgs_var
  pgtol = pgtol_lbfgs_var
  nbd = 0  ! Values are unbounded
  task = 'START'
  iter = 1

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_LBFGS', debug, &
          'Solver config: m    ', m
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_LBFGS', debug, &
          'Solver config: factr', factr
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_LBFGS', debug, &
          'Solver config: pgtol', pgtol
  END IF
  

! ***************************
! ***   Iterative solving ***
! ***************************

  ! Prepare arrays for iterations
  ALLOCATE(gradJ_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_cvec_p)

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
     CALL setulb(dim_cvec_p, m, v_p, lvec, uvec, nbd, &
          J_tot, gradJ_p, factr, pgtol, &
          wa, iwa, task, iprint,&
          csave, lsave, isave, dsave )
     CALL PDAF_timeit(54, 'old')


! ********************************
! ***   Evaluate cost function ***
! ********************************

     CALL PDAF_timeit(53, 'new')
     CALL PDAF_3dvar_costf_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p, &
          obs_p, dy_p, v_p, J_tot, gradJ_p, &
          U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          opt_parallel)
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
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_LBFGS -- END'

END SUBROUTINE PDAF_3dvar_optim_lbfgs

!-------------------------------------------------------------------------------
!> Optimization loop for CG+
!!
!! Optimization routine for 3D-Var using the CG+ solver
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
SUBROUTINE PDAF_3dvar_optim_cgplus(step, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, screen)

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
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_p             !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p        !< Size of control vector
  INTEGER, INTENT(in) :: dim_obs_p         !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)    !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)     !< Background innovation
  REAL, INTENT(inout) :: v_p(dim_cvec_p)   !< Control vector
  INTEGER, INTENT(in) :: opt_parallel      !< Whether to use a decomposed control vector
  INTEGER, INTENT(in) :: screen            !< Verbosity flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &               !< Provide product R^-1 A
       U_cvt, &                            !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                        !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &                     !< Linearized observation operator
       U_obs_op_adj                        !< Adjoint observation operator

! *** local variables ***
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  INTEGER :: optiter                   ! Additional iteration counter

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
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_CGPLUS -- START'

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
  ALLOCATE(d(dim_cvec_p), w(dim_cvec_p))
  ALLOCATE(gradJ_p(dim_cvec_p), gradJ_old_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4*dim_cvec_p)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_CGPLUS', debug, &
          'Solver config: method  ', method
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_CGPLUS', debug, &
          'Solver config: restarts', irest
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_CGPLUS', debug, &
          'Solver config: EPS     ', EPS
  END IF
  

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
        CALL PDAF_3dvar_costf_cvt(step, optiter, dim_p, dim_cvec_p, dim_obs_p, &
             obs_p, dy_p, v_p, J_tot, gradJ_p, &
             U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
             opt_parallel)
        CALL PDAF_timeit(53, 'old')
     END IF


! ***************************
! ***   Optimize with CG+ ***
! ***************************

     CALL PDAF_timeit(54, 'new')
     IF (opt_parallel==0) THEN
        CALL CGFAM(dim_cvec_p, v_p, J_tot, gradJ_p, D, gradJ_old_p, IPRINT, EPS, W,  &
             iflag, IREST, METHOD, FINISH)
     ELSE
        CALL CGFAM_mpi(dim_cvec_p, v_p, J_tot, gradJ_p, D, gradJ_old_p, IPRINT, EPS, W,  &
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
           IF(i > dim_cvec_p) THEN
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

  DEALLOCATE(gradJ_p)
  DEALLOCATE(d, gradJ_old_p, w)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_CGPLUS -- END'

END SUBROUTINE PDAF_3dvar_optim_cgplus

!-------------------------------------------------------------------------------
!> Optimization loop for parallelized CG
!!
!! Optimization routine for 3D-Var using direct
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
SUBROUTINE PDAF_3dvar_optim_cg(step, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, screen)

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
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: dim_p            !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p       !< Size of control vector
  INTEGER, INTENT(in) :: dim_obs_p        !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)   !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)    !< Background innovation
  REAL, INTENT(inout) :: v_p(dim_cvec_p)  !< Control vector
  INTEGER, INTENT(in) :: opt_parallel     !< Whether to use a decomposed control vector
  INTEGER, INTENT(in) :: screen           !< Verbosity flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &              !< Provide product R^-1 A
       U_cvt, &                           !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                       !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &                    !< Linearized observation operator
       U_obs_op_adj                       !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i, iter                   ! Iteration counter
  INTEGER :: maxiter=200               ! maximum number of iterations
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: J_tot, J_old                 ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: hessJd_p(:)     ! Hessian times v
  REAL :: gprod_p, dprod_p, gprod_new_p  ! temporary variables for step size computation
  REAL :: gprod, dprod, gprod_new      ! temporary variables for step size computation
  REAL :: alpha, beta                  ! step sizes
  REAL, ALLOCATABLE :: d_p(:)          ! descent direction
  REAL, ALLOCATABLE :: v_new_p(:)      ! iterated control vector
  REAL, ALLOCATABLE :: gradJ_new_p(:)  ! iterated gradient
  REAL, ALLOCATABLE :: d_new_p(:)      ! iterated descent direction
  REAL :: eps=1.0e-6                   ! Convergence condition value


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_CG -- START'

  maxiter = maxiter_cg_var    ! Maximum number of iterations (default=200)
  eps = eps_cg_var            ! Convergence limit (default=1.0e-6)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_CG', debug, &
          'Solver config: maxiter ', maxiter
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_optim_CG', debug, &
          'Solver config: EPS     ', EPS
  END IF

  ! Prepare arrays for iterations
  ALLOCATE(gradJ_p(dim_cvec_p))
  ALLOCATE(hessJd_p(dim_cvec_p))
  ALLOCATE(d_p(dim_cvec_p))
  ALLOCATE(v_new_p(dim_cvec_p), gradJ_new_p(dim_cvec_p), d_new_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 6*dim_cvec_p)

  ! Initialize numbers
  J_tot = 0.0
  d_p = 0.0
  

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
     CALL PDAF_3dvar_costf_cg_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p, &
          obs_p, dy_p, v_p, d_p, J_tot, gradJ_p, hessJd_p, &
          U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          opt_parallel)
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
        DO i=1, dim_cvec_p
           gprod_p = gprod_p + gradJ_p(i)*gradJ_p(i)
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
     DO i=1, dim_cvec_p
        dprod_p = dprod_p + d_p(i)*hessJd_p(i)
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
     v_new_p = v_p + alpha * d_p

     ! Update gradient
     gradJ_new_p = gradJ_p + alpha * hessJd_p

     ! Compute step size beta for update of descent direction
     gprod_new_p = 0.0
     DO i=1, dim_cvec_p
        gprod_new_p = gprod_new_p + gradJ_new_p(i)*gradJ_new_p(i)
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
     d_new_p = - gradJ_new_p + beta * d_p
     
     ! prepare next iteration
     gradJ_p = gradJ_new_p
     d_p = d_new_p
     v_p = v_new_p
     gprod = gprod_new

     CALL PDAF_timeit(54, 'old')

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(gradJ_p)
  DEALLOCATE(hessJd_p, d_p, v_new_p, gradJ_new_p, d_new_p)
  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_optim_CG -- END'

END SUBROUTINE PDAF_3dvar_optim_cg

!-------------------------------------------------------------------------------
!> Evaluate cost function and its gradient
!!
!! Routine to evaluate the cost function and its gradient
!! for the incremental 3D-Var with variable transformation
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
SUBROUTINE PDAF_3dvar_costf_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, J_tot, gradJ, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel)

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
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: iter             !< Optimization iteration
  INTEGER, INTENT(in) :: dim_p            !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p       !< PE-local size of control vector
  INTEGER, INTENT(in) :: dim_obs_p        !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)   !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)    !< background innovation
  REAL, INTENT(in)  :: v_p(dim_cvec_p)    !< control vector
  REAL, INTENT(out) :: J_tot              !< on exit: Value of cost function
  REAL, INTENT(out) :: gradJ(dim_cvec_p)  !< on exit: PE-local gradient of J
  INTEGER, INTENT(in) :: opt_parallel     !< Whether to use a decomposed control vector

! *** External subroutines ***
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &              !< Provide product R^-1 A
       U_cvt, &                           !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                       !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &                    !< Linearized observation operator
       U_obs_op_adj                       !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: Vv_p(:)         ! PE-local product V deltav
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local product HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradJ (partial sums)
  REAL :: J_B_p, J_B, J_obs_p, J_obs   ! Cost function terms


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- START: iteration', iter

  ! Allocate arrays
  ALLOCATE(Vv_p(dim_p))
  ALLOCATE(HVv_p(dim_obs_p))
  ALLOCATE(RiHVv_p(dim_obs_p, 1))
  ALLOCATE(gradJ_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p + dim_cvec_p + dim_p)


! *******************************************
! ***   Observation part of cost function ***
! *******************************************

  CALL PDAF_timeit(55, 'new')

  CALL PDAF_timeit(56, 'new')

  ! Apply V to control vector v_p
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt'

  CALL PDAF_timeit(60, 'new')
  CALL U_cvt(iter, dim_p, dim_cvec_p, v_p, Vv_p)
  CALL PDAF_timeit(60, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'state increment after CVT dX(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX dX', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply linearized observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'observed state after CVT (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX HdX', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! HVv - dy 
  CALL PDAF_timeit(51, 'new')
  HVv_p = HVv_p - dy_p
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'process local residual d (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX d', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVv = Rinv HVv                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'R^-1 d (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
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
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'process local observation cost J_obs', J_obs_p

  ! Get global value
  CALL MPI_Allreduce(J_obs_p, J_obs, 1, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

  IF (debug>0) &
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'global observation cost J_obs', J_obs

  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(56, 'old')


! ******************************************
! ***   Background part of cost function ***
! ******************************************

  CALL PDAF_timeit(57, 'new')
  CALL PDAF_timeit(51, 'new')

  J_B_p = 0.0
  DO i = 1, dim_cvec_p
     J_B_p = J_B_p + v_p(i)*v_p(i)
  END DO
  J_B_p = 0.5 * J_B_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
       'process local CV cost J_B', J_B_p

  IF (opt_parallel==1) THEN
     ! Get global value
     CALL MPI_Allreduce(J_B_p, J_B, 1, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     J_B = J_B_p
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
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
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_adj'

  CALL PDAF_timeit(65, 'new')
  Vv_p = 0.0
  CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
  CALL PDAF_timeit(65, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'H^TR^-1 d (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX H^TR^-1 d', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply V^T to vector
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt_adj'

  CALL PDAF_timeit(62, 'new')
  CALL U_cvt_adj(iter, dim_p, dim_cvec_p, Vv_p, gradJ)
  CALL PDAF_timeit(62, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'CVT(H^TR^-1 d) (1:min(dim_cvec_p,6))', gradJ(1:min(dim_cvec_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX CVT(H^TR^-1 d)', MINVAL(gradJ), MAXVAL(gradJ)
  END IF

  ! Complete gradient adding v_p
  CALL PDAF_timeit(51, 'new')
  gradJ = v_p + gradJ
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'process local gradient gradJ (1:min(dim_cvec_p,6))', gradJ(1:min(dim_cvec_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX gradJ', MINVAL(gradJ), MAXVAL(gradJ)
  END IF

  CALL PDAF_timeit(58, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Vv_p, HVv_p, RiHVv_p, gradJ_p)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- END'

END SUBROUTINE PDAF_3dvar_costf_cvt

!-------------------------------------------------------------------------------
!> Evaluate cost function, its gradient and Hessian
!!
!! Routine to evaluate the cost function, its gradient, and
!! the product of its Hessian time descent direction
!! for the incremental 3D-Var with variable transformation.
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
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_3dvar_costf_cg_cvt(step, iter, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, d_p, J_tot, gradJ, hessJd, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel)

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
  INTEGER, INTENT(in) :: step               !< Current time step
  INTEGER, INTENT(in) :: iter               !< CG iteration
  INTEGER, INTENT(in) :: dim_p              !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p         !< PE-local size of control vector
  INTEGER, INTENT(in) :: dim_obs_p          !< PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)     !< Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)      !< Background innovation
  REAL, INTENT(in)  :: v_p(dim_cvec_p)      !< Control vector
  REAL, INTENT(inout) :: d_p(dim_cvec_p)    !< CG descent direction
  REAL, INTENT(out) :: J_tot                !< on exit: Value of cost function
  REAL, INTENT(out) :: gradJ(dim_cvec_p)    !< on exit: gradient of J
  REAL, INTENT(out) :: hessJd(dim_cvec_p)   !< on exit: Hessian of J times d_p
  INTEGER, INTENT(in) :: opt_parallel       !< Whether to use a decomposed control vector

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &                !< Provide product R^-1 A
       U_cvt, &                             !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                         !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &                      !< Linearized observation operator
       U_obs_op_adj                         !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: Vv_p(:)         ! PE-local product V deltav
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local produce HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradJ (partial sums)
  REAL :: J_B_p, J_B, J_obs_p, J_obs   ! Cost function terms


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- START: iteration', iter

  ! Allocate arrays
  ALLOCATE(Vv_p(dim_p))
  ALLOCATE(HVv_p(dim_obs_p))
  ALLOCATE(RiHVv_p(dim_obs_p, 1))
  ALLOCATE(gradJ_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p + dim_cvec_p + dim_p)


! *******************************************
! ***   Observation part of cost function ***
! *******************************************

  CALL PDAF_timeit(55, 'new')

  CALL PDAF_timeit(56, 'new')

  ! Apply V to control vector v_p
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt'

  CALL PDAF_timeit(60, 'new')
  CALL U_cvt(iter, dim_p, dim_cvec_p, v_p, Vv_p)
  CALL PDAF_timeit(60, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'state increment after CVT dX(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX dX', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply linearized observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'observed increment HdX(1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX HdX', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! HVv - dy 
  CALL PDAF_timeit(51, 'new')
  HVv_p = HVv_p - dy_p
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'process local residual d (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX d', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVv = Rinv HVv                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'R^-1 d (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
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
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
       'process local observation cost J_obs', J_obs_p

  ! Get global value
  CALL MPI_Allreduce(J_obs_p, J_obs, 1, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
       'global observation cost J_obs', J_obs

  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(56, 'old')


! ******************************************
! ***   Background part of cost function ***
! ******************************************

  CALL PDAF_timeit(57, 'new')
  CALL PDAF_timeit(51, 'new')

  J_B_p = 0.0
  DO i = 1, dim_cvec_p
     J_B_p = J_B_p + v_p(i)*v_p(i)
  END DO
  J_B_p = 0.5 * J_B_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
       'process local CV cost J_B', J_B_p

  IF (opt_parallel==1) THEN
     ! Get global value
     CALL MPI_Allreduce(J_B_p, J_B, 1, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     J_B = J_B_p
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
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
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_adj'

     CALL PDAF_timeit(65, 'new')
     Vv_p = 0.0
     CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
     CALL PDAF_timeit(65, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'H^TR^-1 d (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'MIN/MAX H^TR^-1 d', MINVAL(Vv_p), MAXVAL(Vv_p)
     END IF

     ! Apply V^T to vector
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt_adj'

     CALL PDAF_timeit(62, 'new')
     CALL U_cvt_adj(iter, dim_p, dim_cvec_p, Vv_p, gradJ)
     CALL PDAF_timeit(62, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'CVT(H^TR^-1 d) (1:min(dim_cvec_p,6))', gradJ(1:min(dim_cvec_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'MIN/MAX CVT(H^TR^-1 d)', MINVAL(gradJ), MAXVAL(gradJ)
     END IF

     ! Complete gradient adding v_p
     CALL PDAF_timeit(51, 'new')
     gradJ = v_p + gradJ
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(58, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'process local gradient gradJ (1:min(dim_cvec_p,6))', gradJ(1:min(dim_cvec_p,6))
        WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
             'MIN/MAX gradJ', MINVAL(gradJ), MAXVAL(gradJ)
     END IF

  END IF


! *****************************************************
! ***   Compute Hessian times direction vector d_p  ***
! *****************************************************

  CALL PDAF_timeit(59, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- compute Hessian times direction cvp'

  ! Initialize descent direction d_p at first iteration
  IF (iter==1) THEN
     CALL PDAF_timeit(51, 'new')
     d_p = - gradJ
     CALL PDAF_timeit(51, 'old')
  END IF

  ! Apply V to control vector v_p
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt'

  CALL PDAF_timeit(60, 'new')
  CALL U_cvt(-iter, dim_p, dim_cvec_p, d_p, Vv_p)
  CALL PDAF_timeit(60, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'direction in state space after CVT dp(1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX dp', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_lin'

  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'observed direction Hdp (1:min(dim_obs_p,6))', HVv_p(1:min(dim_obs_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX Hdp', MINVAL(HVv_p), MAXVAL(HVv_p)
  END IF

  ! ***                RiHVd = Rinv HVd                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  ! *** RiHVd is stored in RiHVv
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call prodRinvA'

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'R^-1 Hdp (1:min(dim_obs_p,6))', RiHVv_p(1:min(dim_obs_p,6),1)
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX R^-1 Hdp', MINVAL(RiHVv_p), MAXVAL(RiHVv_p)
  END IF

  ! Apply adjoint of observation operator
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call obs_op_adj'

  CALL PDAF_timeit(65, 'new')
  Vv_p = 0.0
  CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
  CALL PDAF_timeit(65, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'H^TR^-1 dp (1:min(dim_p,6))', Vv_p(1:min(dim_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX H^TR^-1 dp', MINVAL(Vv_p), MAXVAL(Vv_p)
  END IF

  ! Apply V^T to vector
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- call cvt_adj'

  CALL PDAF_timeit(62, 'new')
  CALL U_cvt_adj(-iter, dim_p, dim_cvec_p, Vv_p, hessJd)
  CALL PDAF_timeit(62, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'CVT(H^TR^-1 dp) (1:min(dim_cvec_p,6))', hessJd(1:min(dim_cvec_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX CVT(H^TR^-1 dp)', MINVAL(hessJd), MAXVAL(hessJd)
  END IF

  ! Add d_p to complete Hessian times d_p
  CALL PDAF_timeit(51, 'new')
  hessJd = hessJd + d_p
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'Hessian times dp (1:min(dim_cvec_p,6))', hessJd(1:min(dim_cvec_p,6))
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_costf_cvt:', debug, &
          'MIN/MAX Hessian times dp', MINVAL(hessJd), MAXVAL(hessJd)
  END IF

  CALL PDAF_timeit(59, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Vv_p, HVv_p, RiHVv_p, gradJ_p)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_costf_cvt -- END'

END SUBROUTINE PDAF_3dvar_costf_cg_cvt

END MODULE PDAF_3dvar_optim
