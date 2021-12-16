!$Id$
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
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, v_p, Vv_p)

  USE mod_assimilation, &     ! Assimilation variables
       ONLY: mcols_cvec_ens, Vmat_ens_p

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
  INTEGER :: i, member, row          ! Counters
  REAL :: fact                       ! Scaling factor
  REAL :: invdimens                  ! Inverse ensemble size


! *************************************************
! *** Convert control vector to state increment ***
! *** by computing   Vmat v_p                   ***
! *** Here, Vmat is represented by the ensemble ***
! *************************************************

  ! At beginning of iterations
  firstiter: IF (iter==1) THEN

     ! *** Generate control vector transform matrix ***

     fact = 1.0/SQRT(REAL(dim_cvec_ens-1))

     IF (ALLOCATED(Vmat_ens_p)) DEALLOCATE(Vmat_ens_p)
     ALLOCATE(Vmat_ens_p(dim_p, dim_cvec_ens))

     Vv_p = 0.0
     invdimens = 1.0 / REAL(dim_ens)
     DO member = 1, dim_ens
        DO row = 1, dim_p
           Vv_p(row) = Vv_p(row) + invdimens * ens_p(row, member)
        END DO
     END DO

     DO member = 1, dim_ens
        Vmat_ens_p(:,member) = fact*(ens_p(:,member) - Vv_p(:))
     END DO

     ! Fill additional columns (if Vmat_ens_p holds multiple sets of localized ensembles)
     ! This simulates what would be done with localization (without actually localizing here)
     DO i = 2, mcols_cvec_ens
        DO member = (i-1)*dim_ens+1, i*dim_ens
           Vmat_ens_p(:,member) = Vmat_ens_p(:,member-(i-1)*dim_ens)
        END DO
     END DO

  END IF firstiter

  ! Transform control variable to state increment
  CALL dgemv('n', dim_p, dim_cvec_ens, 1.0, Vmat_ens_p, &
       dim_p, v_p, 1, 0.0, Vv_p, 1)

END SUBROUTINE cvt_ens_pdaf
