!$Id$
!> Apply adjoint covariance operator to a state vector
!!
!! The routine is called during the analysis step.
!! It has to apply the adjoint covariance operator 
!! (transpose of square root of P) to a vector in
!! state space.
!!
!! For domain decomposition, the action is for
!! the PE-local sub-domain of the state. Thus the
!! covariance operator is applied to a sub-state.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cvec, Vv_p, v_p)

  USE mod_assimilation, &
       ONLY: Vmat_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter          !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
  REAL, INTENT(in)    :: Vv_p(dim_p)   !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec) !< PE-local result vector


! ***********************
! *** Compute V^T v_p ***
! ***********************

  ! Transform control variable to state increment
  CALL dgemv('t', dim_p, dim_cvec, 1.0, Vmat_p, &
       dim_p, Vv_p, 1, 0.0, v_p, 1)

END SUBROUTINE cvt_adj_pdaf
