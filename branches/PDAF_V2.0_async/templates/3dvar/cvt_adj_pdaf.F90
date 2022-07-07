!$Id: cvt_adj_pdaf.F90 906 2021-12-01 17:26:32Z lnerger $
!> Apply adjoint covariance operator to a state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in 3D-Var and Hybrid 3D-Var.
!!
!! The routine is called during the analysis step of
!! 3D-Var or hybrid 3D-Var. It has to apply the 
!! adjoint covariance operator (transpose of square
!! root of B) to a vector in state space.
!!
!! For domain decomposition, the action is for
!! the PE-local sub-domain of the state. Thus the
!! covariance operator is applied to a sub-state.
!! In addition the control vector can also be 
!! distributed (in case of type_opt=12 or 13).
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cvec, Vv_p, v_p)

!  USE mod_assimilation, &     ! Assimilation variables
!       ONLY: Vmat_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter          !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
  REAL, INTENT(in)    :: Vv_p(dim_p)   !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec) !< PE-local result vector

! *** local variables ***



! ***************************************************
! *** Apply covariance operator to a state vector ***
! *** by computing Vmat^T Vv_p                    ***
! ***************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE cvt_adj_pdaf.F90: Implement adjoint control vector transform!'

!  v_p = ??



END SUBROUTINE cvt_adj_pdaf
