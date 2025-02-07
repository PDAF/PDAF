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
!$Id$
!BOP
!
! !ROUTINE: PDAF_hyb3dvar_costf_cvt --- Evaluate cost function and its gradient
!
! !INTERFACE:
SUBROUTINE PDAF_hyb3dvar_costf_cvt(step, iter, dim_p, dim_ens, &
     dim_cv_p, dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
     dy_p, v_par_p, v_ens_p, v_p, J_tot, gradJ, &
     U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
     U_obs_op_lin, U_obs_op_adj, opt_parallel, beta)

! !DESCRIPTION:
! Routine to evaluate the cost function and its gradient
! for the incremental hybrid 3D-Var with variable transformation.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
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
  USE PDAF_mod_filter, &
       ONLY: debug

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                   ! Current time step
  INTEGER, INTENT(in) :: iter                   ! Optimization iteration
  INTEGER, INTENT(in) :: dim_p                  ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                ! ensemble size
  INTEGER, INTENT(in) :: dim_cv_p               ! Size of control vector (full)
  INTEGER, INTENT(in) :: dim_cv_par_p           ! Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cv_ens_p           ! Size of control vector (ensemble part)
  INTEGER, INTENT(in) :: dim_obs_p              ! PE-local dimension of observation vector
  REAL, INTENT(in)  :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)         ! Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)          ! background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p)  ! Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p)  ! Control vector (ensemble part)
  REAL, INTENT(in)  :: v_p(dim_cv_p)            ! Control vector (full)
  REAL, INTENT(out) :: J_tot                    ! on exit: Value of cost function
  REAL, INTENT(out) :: gradJ(dim_cv_p)          ! on exit: PE-local gradient of J (full)
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
