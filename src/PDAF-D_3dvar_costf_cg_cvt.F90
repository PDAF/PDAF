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
! !ROUTINE: PDAF_3dvar_costf_cg_cvt --- Evaluate cost function, its gradient and Hessian
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_costf_cg_cvt(step, iter, dim_ens, dim_obs_p, &
     obs_p, dy_p, HV_p, v_p, d_p, &
     J_tot, gradJ, hessJd, &
     U_prodRinvA, screen)

! !DESCRIPTION:
! Routine to evaluate the cost function, its gradient, and
! the product of its Hessian time descent direction
! for the incremental 3D-Var with variable transformation.
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
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: iter         ! CG iteration
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local dimension of observation vector
  REAL, INTENT(in)  :: obs_p(dim_obs_p)         ! Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)          ! background innovation
  REAL, INTENT(in)  :: HV_p(dim_obs_p,dim_ens)  ! on exit: PE-local forecast state
  REAL, INTENT(in)  :: v_p(dim_ens)             ! control vector
  REAL, INTENT(inout)  :: d_p(dim_ens)          ! CG descent direction
  REAL, INTENT(out) :: J_tot                    ! on exit: Value of cost function
  REAL, INTENT(out) :: gradJ(dim_ens)         ! on exit: gradient of J
  REAL, INTENT(out) :: hessJd(dim_ens)        ! on exit: Hessian of J times d_p
  INTEGER, INTENT(in) :: screen       ! Verbosity flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA             ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_analysis_cvt
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
!EOP

! *** local variables ***
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HVv_p(:)        ! PE-local produce HV deltav
  REAL, ALLOCATABLE :: RiHVv_p(:,:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradJ (partial sums)
  REAL, ALLOCATABLE :: hessJd_p(:)     ! PE-local part of hessJd (partial sums)
  REAL :: J_B, J_obs_p, J_obs          ! Cost function terms


! **********************
! *** INITIALIZATION ***
! **********************

  ! Allocate arrays
  ALLOCATE(HVv_p(dim_obs_p))
  ALLOCATE(RiHVv_p(dim_obs_p, 1))
  ALLOCATE(gradJ_p(dim_ens))
  ALLOCATE(hessJd_p(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p + 2*dim_ens)


! *******************************************
! ***   Observation part of cost function ***
! *******************************************

  CALL PDAF_timeit(10, 'new')

  CALL PDAF_timeit(31, 'new')

  ! Multiply HV deltav
  CALL gemvTYPE('n', dim_obs_p, dim_ens, 1.0, HV_p, &
       dim_obs_p, v_p, 1, 0.0, HVv_p, 1)

  ! HVv - dy 
  HVv_p = HVv_p - dy_p

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

  CALL PDAF_timeit(31, 'old')


! ******************************************
! ***   Background part of cost function ***
! ******************************************

  J_B = 0.0
  DO i = 1, dim_ens
     J_B = J_B + v_p(i)*v_p(i)
  END DO
  J_B = 0.5*J_B


! *****************************
! ***   Total cost function ***
! *****************************

  J_tot = J_B + J_obs

  CALL PDAF_timeit(10, 'old')


! **************************
! ***   Compute gradient ***
! **************************

  ! Only at first iteration
  IF (iter==1) THEN

     CALL PDAF_timeit(20, 'new')

     ! Multiplication HV * deltav
     CALL gemvTYPE('t', dim_obs_p, dim_ens, 1.0, HV_p, &
          dim_obs_p, RiHVv_p, 1, 0.0, gradJ_p, 1)

     ! Get vector with global values
     CALL MPI_Allreduce(gradJ_p, gradJ, dim_ens, MPI_REALTYPE, MPI_SUM, &
          COMM_filter, MPIerr)
 
     ! Complete gradient adding v_p
     gradJ = v_p + gradJ

     CALL PDAF_timeit(20, 'old')

  END IF

! *****************************************************
! ***   Compute Hessian times direction vector d_p  ***
! *****************************************************

  ! Initialize descent direction d_p at first iteration
  IF (iter==1) THEN
     d_p = - gradJ
  END IF

  ! Multiply HV d_p (we store this in HVv_p)
  CALL gemvTYPE('n', dim_obs_p, dim_ens, 1.0, HV_p, &
       dim_obs_p, d_p, 1, 0.0, HVv_p, 1)

  ! ***                RiHVd = Rinv HVd                
  ! *** This is implemented as a subroutine thus that
  ! *** Rinv does not need to be allocated explicitly.
  ! *** RiHVd is stored in RiHVv

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, HVv_p, RiHVv_p)
  CALL PDAF_timeit(48, 'old')

  ! Multiply HV RiHVv_p
  CALL gemvTYPE('t', dim_obs_p, dim_ens, 1.0, HV_p, &
       dim_obs_p, RiHVv_p, 1, 0.0, hessJd_p, 1)

  ! Get vector with global values
  CALL MPI_Allreduce(hessJd_p, hessJd, dim_ens, MPI_REALTYPE, MPI_SUM, &
       COMM_filter, MPIerr)

  ! Add d_p to complete Hessian times d_p
  hessJd = hessJd + d_p


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(HVv_p, RiHVv_p, gradJ_p, hessJd_p)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_costf_cg_cvt
