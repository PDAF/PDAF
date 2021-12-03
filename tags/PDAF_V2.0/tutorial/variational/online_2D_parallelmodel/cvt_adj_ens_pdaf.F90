!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!> Apply adjoint ensemble covariance operator to a state vector
!!
!! The routine is called during the analysis step.
!! It has to apply the adjoint covariance operator 
!! (transpose of square root of P) to a vector in
!! state space.
!!
!! For ensemble 3D-Var the ensemble representation
!! of the covariance operator is used.
!!
!! For domain decomposition, the action is for
!! the PE-local sub-domain of the state. Thus the
!! covariance operator is applied to a sub-state.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, Vcv_p, cv_p)

  USE mpi                     ! MPI
  USE mod_assimilation, &     ! Assimilation variables
       ONLY: Vmat_ens_p, dim_cvec_ens, off_cv_ens_p, type_opt
  USE mod_parallel_pdaf, &    ! PDAF parallelization variables
       ONLY: COMM_filter, MPIerr, mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter                 !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p                !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens              !< Ensemble size
  INTEGER, INTENT(in) :: dim_cv_ens_p         !< PE-local dimension of control vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)   !< PE-local ensemble
  REAL, INTENT(in)    :: Vcv_p(dim_p)         !< PE-local input vector
  REAL, INTENT(inout) :: cv_p(dim_cv_ens_p)   !< PE-local result vector

! *** local variables ***
  INTEGER :: i                       ! Counters
  REAL, ALLOCATABLE :: cv_g(:)       ! Global control vector
  REAL, ALLOCATABLE :: cv_g_part(:)  ! Global control vector (partial sums)


! *****************************************************
! *** Compute Vmat^T x_p with x_p some state vector ***
! *****************************************************

  ALLOCATE(cv_g_part(dim_cvec_ens))

  IF (type_opt==12 .OR. type_opt==13) THEN

     ! For the domain-decomposed solvers

     ! Initialize distributed vector on control space
     ALLOCATE(cv_g(dim_cvec_ens))

     ! Transform control variable to state increment 
     ! - global vector of partial sums
     CALL dgemv('t', dim_p, dim_cvec_ens, 1.0, Vmat_ens_p, &
          dim_p, Vcv_p, 1, 0.0, cv_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(cv_g_part, cv_g, dim_cvec_ens, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)
     
     ! Select PE-local part of control vector
     DO i = 1, dim_cv_ens_p
        cv_p(i) = cv_g(i + off_cv_ens_p(mype_filter+1))
     END DO

     DEALLOCATE(cv_g)

  ELSE

     ! Without domain-decomposition

     ! Transform control variable to state increment
     CALL dgemv('t', dim_p, dim_cv_ens_p, 1.0, Vmat_ens_p, &
          dim_p, Vcv_p, 1, 0.0, cv_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(cv_g_part, cv_p, dim_cvec_ens, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)

  END IF


! *** Clean up ***

  DEALLOCATE(cv_g_part)

END SUBROUTINE cvt_adj_ens_pdaf
