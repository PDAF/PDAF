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
! !ROUTINE: PDAF_hyb3dvar_costf_cg_cvt --- Evaluate cost function, its gradient and Hessian
!
! !INTERFACE:
SUBROUTINE PDAF_hyb3dvar_costf_cg_cvt(step, iter, dim_p, dim_ens, &
     dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
     dy_p, v_par_p, v_ens_p, d_par_p, d_ens_p, &
     J_tot, gradJ_par, gradJ_ens, hessJd_par, hessJd_ens, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
     U_obs_op_lin, U_obs_op_adj, opt_parallel, beta)

! !DESCRIPTION:
! Routine to evaluate the cost function, its gradient, and
! the product of its Hessian time descent direction
! for the incremental hybrid 3D-Var with variable transformation.
!
! The subroutine distinguishes two cases:
! iter==1
!   In this case all quantities are computed, the 
!   descent direction is initialized from the gradient vector
! iter>1
!   In this case only the cost function value and the product
!   of the Hessian times descent direction are computed.
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
       ONLY: MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                   ! Current time step
  INTEGER, INTENT(in) :: iter                   ! Optimization iteration
  INTEGER, INTENT(in) :: dim_ens                ! ensemble size
  INTEGER, INTENT(in) :: dim_p                  ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_cv_par_p           ! Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cv_ens_p           ! Size of control vector (ensemble part)
  INTEGER, INTENT(in) :: dim_obs_p              ! PE-local dimension of observation vector
  REAL, INTENT(in)  :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)         ! Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)          ! Background innovation
  REAL, INTENT(in)  :: v_par_p(dim_cv_par_p)    ! Control vector (parameterized part)
  REAL, INTENT(in)  :: v_ens_p(dim_cv_ens_p)    ! Control vector (ensemble part)
  REAL, INTENT(inout) :: d_par_p(dim_cv_par_p)  ! CG descent direction (parameterized part)
  REAL, INTENT(inout) :: d_ens_p(dim_cv_ens_p)  ! CG descent direction (ensemble part)
  REAL, INTENT(out) :: J_tot                    ! on exit: Value of cost function
  REAL, INTENT(out) :: gradJ_par(dim_cv_par_p)  ! on exit: gradient of J (parameterized part)
  REAL, INTENT(out) :: gradJ_ens(dim_cv_ens_p)  ! on exit: gradient of J (ensemble part)
  REAL, INTENT(out) :: hessJd_par(dim_cv_par_p) ! on exit: Hessian of J times d_p (parameterized part)
  REAL, INTENT(out) :: hessJd_ens(dim_cv_ens_p) ! on exit: Hessian of J times d_p (ensemble part)
  INTEGER, INTENT(in) :: opt_parallel           ! Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta                      ! Hybrid weight

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   ! Provide product R^-1 A
       U_cvt, &                ! Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            ! Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            ! Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        ! Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         ! Linearized observation operator
       U_obs_op_adj            ! Adjoint observation operator

! !CALLING SEQUENCE:
! Called by: PDAF_hyb3dvar_analysis_cvt
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
!EOP

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
     CALL PDAF_timeit(60, 'new')
     CALL U_cvt(iter, dim_p, dim_cv_par_p, v_par_p, Vv_p)
     CALL PDAF_timeit(60, 'old')
  END IF

  ! ensemble
  IF (dim_cv_ens_p>0) THEN
     CALL PDAF_timeit(61, 'new')
     CALL U_cvt_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, v_ens_p, Vv_ens_p)
     CALL PDAF_timeit(61, 'old')
  END IF

  Vv_p = sombeta*Vv_p + sbeta*Vv_ens_p

  ! Apply linearized observation operator
  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  ! HVv - dy
  CALL PDAF_timeit(51, 'new')
  HVv_p = HVv_p - dy_p
  CALL PDAF_timeit(51, 'old')

  ! ***                RiHVv = Rinv HVv                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  ! ***  Compute  J_obs ***

  CALL PDAF_timeit(51, 'new')

  J_obs_p = 0.0
  DO i = 1, dim_obs_p
     J_obs_p = J_obs_p + HVv_p(i)*RiHVv_p(i,1)
  END DO

  J_obs_p = 0.5*J_obs_p

  ! Get global value
  CALL MPI_Allreduce(J_obs_p, J_obs, 1, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

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
  DO i = 1, dim_cv_ens_p
     J_B_p = J_B_p + v_ens_p(i)*v_ens_p(i)
  END DO

  IF (opt_parallel==1) THEN
     ! Get global value
     CALL MPI_Allreduce(J_B_p, J_B, 1, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     J_B = J_B_p
  END IF

  J_B = 0.5*J_B

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
     CALL PDAF_timeit(65, 'new')
     Vv_p = 0.0
     CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
     CALL PDAF_timeit(65, 'old')

     ! Apply V^T to vector
     IF (dim_cv_par_p>0) THEN
        CALL PDAF_timeit(62, 'new')
        CALL U_cvt_adj(iter, dim_p, dim_cv_par_p, Vv_p, gradJ_par)
        CALL PDAF_timeit(62, 'old')
     END IF
     IF (dim_cv_ens_p>0) THEN
        CALL PDAF_timeit(63, 'new')
        CALL U_cvt_adj_ens(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, gradJ_ens)
        CALL PDAF_timeit(63, 'old')
     END IF

     ! Complete gradient adding v_p
     CALL PDAF_timeit(51, 'new')
     IF (dim_cv_par_p>0) THEN
        gradJ_par = v_par_p + sombeta*gradJ_par
     END IF
     IF (dim_cv_ens_p>0) THEN
        gradJ_ens = v_ens_p + sbeta*gradJ_ens
     END IF
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(58, 'old')

  END IF


! *****************************************************
! ***   Compute Hessian times direction vector d_p  ***
! *****************************************************

  CALL PDAF_timeit(59, 'new')

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
     CALL PDAF_timeit(60, 'new')
     CALL U_cvt(-iter, dim_p, dim_cv_par_p, d_par_p, Vv_p)
     CALL PDAF_timeit(60, 'old')
  END IF
  IF (dim_cv_ens_p>0) THEN
     CALL PDAF_timeit(61, 'new')
     CALL U_cvt_ens(-iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, d_ens_p, Vv_ens_p)
     CALL PDAF_timeit(61, 'old')
  END IF

  Vv_p = sombeta*Vv_p + sbeta*Vv_ens_p

  ! Apply observation operator
  CALL PDAF_timeit(64, 'new')
  CALL U_obs_op_lin(step, dim_p, dim_obs_p, Vv_p, HVv_p)
  CALL PDAF_timeit(64, 'old')

  ! ***                RiHVd = Rinv HVd                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  ! *** RiHVd is stored in RiHVv

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  ! Apply adjoint of observation operator
  CALL PDAF_timeit(65, 'new')
  Vv_p = 0.0
  CALL U_obs_op_adj(step, dim_p, dim_obs_p, RiHVv_p, Vv_p)
  CALL PDAF_timeit(65, 'old')

  ! Apply V^T to vector
  IF (dim_cv_par_p>0) THEN
     CALL PDAF_timeit(62, 'new')
     CALL U_cvt_adj(-iter, dim_p, dim_cv_par_p, Vv_p, hessJd_par)
     CALL PDAF_timeit(62, 'old')
  END IF
  IF (dim_cv_ens_p>0) THEN
     CALL PDAF_timeit(63, 'new')
     CALL U_cvt_adj_ens(-iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, hessJd_ens)
     CALL PDAF_timeit(63, 'old')
  END IF

  ! Add d_p to complete Hessian times d_p
  CALL PDAF_timeit(51, 'new')
    IF (dim_cv_par_p>0) hessJd_par = sombeta*hessJd_par + d_par_p
    IF (dim_cv_ens_p>0) hessJd_ens = sbeta*hessJd_ens + d_ens_p
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(59, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Vv_p, Vv_ens_p, HVv_p, RiHVv_p)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_hyb3dvar_costf_cg_cvt
