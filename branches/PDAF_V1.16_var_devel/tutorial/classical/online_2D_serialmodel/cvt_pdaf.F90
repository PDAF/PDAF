!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvt_pdaf --- Generate matrix of localized ensemble perturbations
!
! !INTERFACE:
SUBROUTINE cvt_pdaf(iter, dim_p, dim_cvec, v_p, Vv_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to apply the covariance operator (square
! root of P) to a vector in control space.
!
! For domain decomposition, the action is on
! the control vector for the PE-local part of
! the sub-state vector for the PE-local domain.
!
! This code variant uses an explicit array holding
! the covariance operator as a matrix.
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
  INTEGER, INTENT(in) :: iter          ! Iteration of optimization
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

END SUBROUTINE cvt_pdaf
