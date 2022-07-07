!$Id: cvt_pdaf.F90 906 2021-12-01 17:26:32Z lnerger $
!> Apply covariance operator to a control vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in 3D-Var and Hybrid 3D-Var.
!! 
!! The routine is called during the analysis step of
!! 3D-Var or hybrid 3D-Var. It has to apply the 
!! covariance operator (square root of P) to a vector 
!! in control space.
!!
!! For domain decomposition, the action is on
!! the control vector for the PE-local part of
!! the sub-state vector for the PE-local domain.
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
SUBROUTINE cvt_pdaf(iter, dim_p, dim_cvec, v_p, Vv_p)

!  USE mod_assimilation, &     ! Assimilation variables
!       ONLY: Vmat_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter          !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
  REAL, INTENT(in)    :: v_p(dim_cvec) !< PE-local control vector
  REAL, INTENT(inout) :: Vv_p(dim_p)   !< PE-local result vector

! *** local variables ***



! ***************************************************
! *** Apply covariance operator to control vector ***
! *** by computing Vmat v_p                       ***
! ***************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE cvt_pdaf.F90: Implement control vector transform!'

!  Vv_p = ??

END SUBROUTINE cvt_pdaf
