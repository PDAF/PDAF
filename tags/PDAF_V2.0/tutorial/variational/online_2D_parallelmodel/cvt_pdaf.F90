!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!> Apply covariance operator to a control vector
!!
!! The routine is called during the analysis step of
!! 3D-Var or hybrid 3D-Var. It has to apply the 
!! covariance operator (square root of P) to a vector 
!! in control space.
!!
!! For domain decomposition, the action is on
!! the control vector for the PE-local part of
!! the sub-state vector for the PE-local domain.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_pdaf(iter, dim_p, dim_cv_p, cv_p, Vcv_p)

  USE mpi                     ! MPI
  USE mod_assimilation, &     ! Assimilation variables
       ONLY: Vmat_p, dim_cvec, dims_cv_p, off_cv_p, type_opt
  USE mod_parallel_pdaf, &    ! PDAF parallelization variables
       ONLY: COMM_filter, MPIerr

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter           !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p          !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cv_p       !< Dimension of control vector
  REAL, INTENT(in)    :: cv_p(dim_cv_p) !< PE-local model state
  REAL, INTENT(inout) :: Vcv_p(dim_p)   !< PE-local result vector

! *** local variables ***
  REAL, ALLOCATABLE :: v_g(:)        ! Global control vector


! *************************
! *** Compute Vmat cv_p ***
! *************************

  IF (type_opt==12 .OR. type_opt==13) THEN

     ! Gather global control vector
     ALLOCATE(v_g(dim_cvec))

     CALL MPI_AllGatherV(cv_p, dim_cv_p, MPI_REAL8, &
          v_g, dims_cv_p, off_cv_p, MPI_REAL8, &
          COMM_filter, MPIerr)

     ! Transform control variable to state increment
     CALL dgemv('n', dim_p, dim_cvec, 1.0, Vmat_p, &
          dim_p, v_g, 1, 0.0, Vcv_p, 1)

     DEALLOCATE(v_g)

  ELSE
     ! Transform control variable to state increment
     CALL dgemv('n', dim_p, dim_cv_p, 1.0, Vmat_p, &
          dim_p, cv_p, 1, 0.0, Vcv_p, 1)

  END IF

END SUBROUTINE cvt_pdaf
