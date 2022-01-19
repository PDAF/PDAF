!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
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
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cvec_p, Vv_p, v_p)

  USE MPI                     ! MPI
  USE mod_assimilation, &     ! Assimilation variables
       ONLY: Vmat_p, dim_cvec, off_v_p, type_opt
  USE mod_parallel_pdaf, &    ! PDAF parallelization variables
       ONLY: COMM_filter, MPIerr, mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter            !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p           !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec_p      !< Dimension of control vector
  REAL, INTENT(in)    :: Vv_p(dim_p)     !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec_p) !< PE-local result vector

! *** local variables ***
  INTEGER :: i                        ! Counter
  REAL, ALLOCATABLE :: v_g(:)         ! Global control vector
  REAL, ALLOCATABLE :: v_g_part(:)    ! Global control vector (partial sums)


! ***************************************************
! *** Apply covariance operator to a state vector ***
! *** by computing Vmat^T Vv_p                    ***
! ***************************************************

  ALLOCATE(v_g_part(dim_cvec))

  IF (type_opt==12 .OR. type_opt==13) THEN

     ! For the domain-decomposed solvers

     ! Initialize distributed vector on control space
     ALLOCATE(v_g(dim_cvec))

     ! Transform control variable to state increment 
     ! - global vector of partial sums
     CALL dgemv('t', dim_p, dim_cvec, 1.0, Vmat_p, &
          dim_p, Vv_p, 1, 0.0, v_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(v_g_part, v_g, dim_cvec, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)
     
     ! Select PE-local part of control vector
     DO i = 1, dim_cvec_p
        v_p(i) = v_g(i + off_v_p(mype_filter+1))
     END DO

     DEALLOCATE(v_g)

  ELSE

     ! Without domain-decomposition

     ! Transform control variable to state increment
     CALL dgemv('t', dim_p, dim_cvec_p, 1.0, Vmat_p, &
          dim_p, Vv_p, 1, 0.0, v_g_part, 1)

     ! Get global vector with global sums
     CALL MPI_Allreduce(v_g_part, v_p, dim_cvec, MPI_REAL8, MPI_SUM, &
          COMM_filter, MPIerr)

  END IF


! *** Clean up ***

  DEALLOCATE(v_g_part)

END SUBROUTINE cvt_adj_pdaf
