!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cov_op_cvec_pdaf --- Generate matrix of localized ensemble perturbations
!
! !INTERFACE:
SUBROUTINE cov_op_cvec_pdaf(dim_p, dim_cvec, v_p, Vv_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: ensemble 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to apply the covariance operator (square
! root of P) to the control vector or the descent
! direction vector of CG.
!
! For domain decomposition, the action is on the
! PE-local part of the control vector and has to 
! provide the sub-state vector for the PE-local 
! domain.
!
! Implementation for the 1D dummy model.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: Vmat_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p         ! PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec      ! Dimension of control vector
  REAL, INTENT(in)    :: v_p(dim_cvec) ! PE-local model state
  REAL, INTENT(inout) :: Vv_p(dim_p)   ! PE-local result vector
!EOP


! *********************
! *** Compute V v_p ***
! *********************

  ! Transform control variable to state increment
  CALL dgemv('n', dim_p, dim_cvec, 1.0, Vmat_p, &
       dim_p, v_p, 1, 0.0, Vv_p, 1)

END SUBROUTINE cov_op_cvec_pdaf
