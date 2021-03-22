!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvt_adj_ens_pdaf --- Apply adjoint covariance operator
!
! !INTERFACE:
SUBROUTINE cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vv_p, v_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: ensemble 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to apply the adjoint covariance operator 
! (transpose of square root of P) to a vector in
! state space.
!
! For ensemble-var the ensemble representation
! of the covariance operator is used.
!
! For domain decomposition, the action is for
! the PE-local sub-domain of the state. Thus the
! covariance operator is applied to a sub-state.
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
       ONLY: mcols_cvec_ens, dim_cvec_ens, off_cv_p, type_opt
  USE mod_parallel, &
       ONLY: MPI_REAL8, COMM_filter, MPI_SUM, MPIerr, mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: iter               ! Iteration of optimization
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            ! Ensemble size
  INTEGER, INTENT(in) :: dim_cv_ens_p       ! PE-local dimension of control vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens) ! PE-local ensemble
  REAL, INTENT(in)    :: Vv_p(dim_p)        ! PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cv_ens_p)  ! PE-local result vector
!EOP

! *** local variables ***
  INTEGER :: i, member, row          ! Counters
  REAL :: fact                       ! Scaling factor
  REAL, ALLOCATABLE :: Vmat_p(:,:)   ! Extended ensemble perturbation matrix
  REAL, ALLOCATABLE :: state_p(:)    ! Ensemble mean state
  REAL, ALLOCATABLE :: v_g(:)        ! Global control vector
  REAL, ALLOCATABLE :: v_g_part(:)   ! Global control vector (partial sums)
  REAL :: invdimens                  ! Inverse ensemble size



! *****************************************************
! *** Compute V^T x_p with x_p is some state vector ***
! *****************************************************

  ! *** Generate control vector transform matrix ***

  ALLOCATE(Vmat_p(dim_p, dim_ens))
  ALLOCATE(state_p(dim_p))

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  DO member = 1, dim_ens
     Vmat_p(:,member) = ens_p(:,member) - state_p(:)
  END DO

  ! Fill additional columns (if Vmat_p holds multiple sets of localized ensenbles)
  DO i = 2, mcols_cvec_ens
     DO member = (i-1)*dim_ens+1, i*dim_ens
        Vmat_p(:,member) = ens_p(:,member-(i-1)*dim_ens)
     END DO
  END DO

  ! Initialize scaling factor
  fact = 1.0/SQRT(REAL(dim_cvec_ens-1))

  ALLOCATE(v_g_part(dim_cvec_ens))

  IF (type_opt/=3) THEN

     ! Transform control variable to state increment
     CALL dgemv('t', dim_p, dim_cv_ens_p, fact, Vmat_p, &
          dim_p, Vv_p, 1, 0.0, v_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(v_g_part, v_p, dim_cvec_ens, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)

  ELSE

     ! Initialize distributed vector on control space
     ALLOCATE(v_g(dim_cvec_ens))

     ! Transform control variable to state increment 
     ! - global vector of partial sums
     CALL dgemv('t', dim_p, dim_cvec_ens, fact, Vmat_p, &
          dim_p, Vv_p, 1, 0.0, v_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(v_g_part, v_g, dim_cvec_ens, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)
     
     ! Select PE-local part of control vector
     DO i = 1, dim_cv_ens_p
        v_p(i) = v_g(i + off_cv_p(mype_filter+1))
     END DO

     DEALLOCATE(v_g)

  END IF


! *** Clean up ***

  DEALLOCATE(Vmat_p, state_p, v_g_part)

END SUBROUTINE cvt_adj_ens_pdaf
