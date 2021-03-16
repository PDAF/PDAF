!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvec2state_ens_pdaf --- Convert control vector to state increment
!
! !INTERFACE:
SUBROUTINE cvec2state_ens_pdaf(step, dim_p, dim_ens, dim_cvec_ens, enspert_p, v_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: 3D-Var
!
! The routine is called during the analysis step
! of the ensemble 3D-Var. It has to transform
! the control vector after the optimization into
! an increment vector in state space.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! Implementation for the 2D online example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: mcols_cvec_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            ! Ensemble size
  INTEGER, INTENT(in) :: dim_cvec_ens       ! Number of columns in HV_p
  REAL, INTENT(in) :: enspert_p(dim_p, dim_ens) ! PE-local ensemble
  REAL, INTENT(in) :: v_p(dim_cvec_ens)     ! PE-local control vector
  REAL, INTENT(inout) :: state_p(dim_p)     ! PE-local state increment
!EOP

! *** local variables ***
  INTEGER :: i, member               ! Counters
  REAL :: fact                       ! Scaling factor
  REAL, ALLOCATABLE :: Vmat_p(:,:)   ! Extended ensemble perturbation matrix


! *************************************************
! *** Convert control vector to state increment ***
! *************************************************

  ALLOCATE(Vmat_p(dim_p, dim_cvec_ens))
  
  DO member = 1, dim_ens
     Vmat_p(:,member) = enspert_p(:,member)
  END DO

  DO i = 2, mcols_cvec_ens
     DO member = (i-1)*dim_ens+1, i*dim_ens
        Vmat_p(:,member) = enspert_p(:,member-(i-1)*dim_ens)
     END DO
  END DO
  
  fact = 1.0/SQRT(REAL(dim_cvec_ens-1))

  ! Transform control variable to state increment
  CALL dgemv('n', dim_p, dim_cvec_ens, fact, Vmat_p, &
       dim_p, v_p, 1, 1.0, state_p, 1)

  DEALLOCATE(Vmat_p)

END SUBROUTINE cvec2state_ens_pdaf
