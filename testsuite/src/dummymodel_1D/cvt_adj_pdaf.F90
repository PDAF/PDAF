!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvt_adj_pdaf --- Apply adjoint covariance operator
!
! !INTERFACE:
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cv_p, Vcv_p, cv_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to apply the adjoint covariance operator 
! (transpose of square root of P) to a vector in
! state space.
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
       ONLY: Vmat_p, dim_cvec, off_cv_p, type_opt
  USE mod_parallel, &
       ONLY: MPI_REAL8, COMM_filter, MPI_SUM, MPIerr, mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: iter           ! Iteration of optimization
  INTEGER, INTENT(in) :: dim_p          ! PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cv_p       ! PE-local dimension of control vector
  REAL, INTENT(in)    :: Vcv_p(dim_p)   ! PE-local input vector
  REAL, INTENT(inout) :: cv_p(dim_cvec) ! PE-local result vector
!EOP

! *** local variables ***
  INTEGER :: i                        ! Counters
  REAL, ALLOCATABLE :: cv_g(:)        ! Global control vector
  REAL, ALLOCATABLE :: cv_g_part(:)   ! Global control vector (partial sums)


! *****************************************************
! *** Compute Vmat^T x_p with x_p some state vector ***
! *****************************************************

  ALLOCATE(cv_g_part(dim_cvec))

  IF (type_opt==12 .OR. type_opt==13) THEN

     ! Initialize distributed vector on control space
     ALLOCATE(cv_g(dim_cvec))

     ! Transform control variable to state increment 
     ! - global vector of partial sums
     CALL dgemv('t', dim_p, dim_cvec, 1.0, Vmat_p, &
          dim_p, Vcv_p, 1, 0.0, cv_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(cv_g_part, cv_g, dim_cvec, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)
     
     ! Select PE-local part of control vector
     DO i = 1, dim_cv_p
        cv_p(i) = cv_g(i + off_cv_p(mype_filter+1))
     END DO

     DEALLOCATE(cv_g)

  ELSE

     ! Transform control variable to state increment
     CALL dgemv('t', dim_p, dim_cv_p, 1.0, Vmat_p, &
          dim_p, Vcv_p, 1, 0.0, cv_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(cv_g_part, cv_p, dim_cvec, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)

  END IF


! *** Clean up ***

  DEALLOCATE(cv_g_part)

END SUBROUTINE cvt_adj_pdaf
