!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!> Apply adjoint ensemble covariance operator to a state vector
!!
!! The routine is called during the analysis step.
!! It has to apply the adjoint covariance operator 
!! (transpose of square root of P) to a vector in
!! state space.
!!
!! For ensemble 3D-Var the ensemble representation
!! of the covariance operator is used.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, Vv_p, v_p)

  USE mod_assimilation, &
       ONLY: Vmat_ens_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter               !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p              !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
  INTEGER, INTENT(in) :: dim_cvec_ens       !< Number of columns in HV_p
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local ensemble
  REAL, INTENT(in)    :: Vv_p(dim_p)        !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec_ens)  !< PE-local result vector


! ***********************
! *** Compute V^T v_p ***
! ***********************

  ! Transform control variable to state increment
  CALL dgemv('t', dim_p, dim_cvec_ens, 1.0, Vmat_ens_p, &
       dim_p, Vv_p, 1, 0.0, v_p, 1)

END SUBROUTINE cvt_adj_ens_pdaf
