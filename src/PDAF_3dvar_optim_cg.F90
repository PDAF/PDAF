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
! !ROUTINE: PDAF_3dvar_optim_cg_ --- Optimization loop for parallelized CG
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_optim_cg(step, dim_p, dim_cvec_p, dim_obs_p, &
     obs_p, dy_p, v_p, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, screen)

! !DESCRIPTION:
! Optimization routine for 3D-Var using direct
! parallelized implementation of CG.
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
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, Comm_filter, MPI_REALTYPE, MPI_SUM, MPIerr

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                  ! Current time step
  INTEGER, INTENT(in) :: dim_p                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_cvec_p            ! Size of control vector
  INTEGER, INTENT(in) :: dim_obs_p             ! PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)        ! Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)         ! Background innovation
  REAL, INTENT(inout) :: v_p(dim_cvec_p)       ! Control vector
  INTEGER, INTENT(in) :: opt_parallel          ! Whether to use a decomposed control vector
  INTEGER, INTENT(in) :: screen                ! Verbosity flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &              ! Provide product R^-1 A
       U_cvt, &                           ! Apply control vector transform matrix to control vector
       U_cvt_adj, &                       ! Apply adjoint control vector transform matrix
       U_obs_op_lin, &                    ! Linearized observation operator
       U_obs_op_adj                       ! Adjoint observation operator

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_analysis_cg_cvt
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER :: i, iter                   ! Iteration counter
  INTEGER :: maxiter                   ! maximum number of iterations
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
  REAL :: eps                          ! Convergence condition value


! **********************
! *** INITIALIZATION ***
! **********************

  maxiter = 200    ! Maximum number of iterations
  eps = 1.0e-6     ! Convergence limit

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

!   IF (mype==0 .AND. screen > 0) &
!        WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

       WRITE (*, *) mype, 'PDAF', '--- OPTIMIZE' 

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

END SUBROUTINE PDAF_3dvar_optim_cg
