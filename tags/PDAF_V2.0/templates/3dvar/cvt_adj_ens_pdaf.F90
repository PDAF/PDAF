!$Id: cvt_adj_ens_pdaf.F90 906 2021-12-01 17:26:32Z lnerger $
!> Apply adjoint ensemble covariance operator to a state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in 3D Ensemble Var and Hybrid 3D-Var.
!!
!! The routine is called during the analysis step of
!! ensemble 3D-Var or hybrid 3D-Var. It has to apply
!! the adjoint ensemble covariance operator (square
!! root of B) to a vector in control space.
!!
!! For domain decomposition, the action is for
!! the PE-local sub-domain of the state. Thus the
!! covariance operator is applied to a sub-state.
!! In addition the control vector can also be 
!! distributed (in case of type_opt=12 or 13).
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, &
     Vv_p, v_p)

!  USE mod_assimilation, &     ! Assimilation variables
!       ONLY: Vmat_ens_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter               !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p              !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
  INTEGER, INTENT(in) :: dim_cvec_ens       !< Number of columns in HV_p
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local ensemble
  REAL, INTENT(in)    :: Vv_p(dim_p)        !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec_ens)  !< PE-local result vector

! *** local variables ***



! ***************************************************
! *** Apply covariance operator to a state vector ***
! *** by computing Vmat^T Vv_p                    ***
! *** Here, Vmat is represented by the ensemble   ***
! ***************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE cvt_adj_ens_pdaf.F90: Implement adjoint ensemble control vector transform!'

  ! Here one could use Vmat_ens_p as ensemble covariance 
  ! Operator if this was readily computed in cvt_ens_pdaf.

!  v_p = ??

END SUBROUTINE cvt_adj_ens_pdaf
