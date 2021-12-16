!$Id: cvt_ens_pdaf.F90 906 2021-12-01 17:26:32Z lnerger $
!> Apply ensemble covariance operator to a control vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in 3D Ensemble Var and Hybrid 3D-Var.
!!
!! The routine is called during the analysis step of
!! ensemble 3D-Var or hybrid 3D-Var. It has to apply
!! the ensemble covariance operator (square root of P)
!! to a vector in control space.
!!
!! For domain decomposition, the action is on
!! the control vector for the PE-local part of
!! the sub-state vector for the PE-local domain.
!! In addition the control vector can also be 
!! distributed (in case of type_opt=12 or 13).
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, &
     v_p, Vv_p)

!   USE mod_assimilation, &     ! Assimilation variables
!        ONLY: mcols_cvec_ens, Vmat_ens_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter               !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p              !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
  INTEGER, INTENT(in) :: dim_cvec_ens       !< Dimension of control vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local ensemble
  REAL, INTENT(in) :: v_p(dim_cvec_ens)     !< PE-local control vector
  REAL, INTENT(inout) :: Vv_p(dim_p)        !< PE-local state increment

! *** local variables ***



! *************************************************
! *** Convert control vector to state increment ***
! *** by computing   Vmat v_p                   ***
! *** Here, Vmat is represented by the ensemble ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE cvt_ens_pdaf.F90: Implement ensemble control vector transform!'

  firstiter: IF (iter==1) THEN
     ! At the beginning of the iteration one could initialize
     ! the scaled ensemble perturbation matrix and apply
     ! localization to initialize Vmat_ens_p. Then, this 
     ! matrix can be used here and in cvt_adj_ens_pdaf
     ! throughout the iterative optimization.
  END IF firstiter

!  Vv_p = ??

END SUBROUTINE cvt_ens_pdaf
