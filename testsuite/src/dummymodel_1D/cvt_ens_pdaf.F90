!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvt_ens_pdaf --- Convert control vector to state increment
!
! !INTERFACE:
SUBROUTINE cvt_ens_pdaf(iter, dim_p, dim_ens, dim_cv_ens_p, ens_p, v_p, Vv_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: ensemble 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to apply the covariance operator (square
! root of P) to a vector in control space.
!
! For ensemble-var the ensemble representation
! of the covariance operator is used.
!
! For domain decomposition, the action is on
! the control vector for the PE-local part of
! the sub-state vector for the PE-local domain.
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
       ONLY: mcols_cvec_ens, dim_cvec_ens, dims_cv_p, off_cv_p, type_opt
  USE mod_parallel, &
       ONLY: MPI_REAL8, COMM_filter, MPIerr, mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: iter               ! Iteration of optimization
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            ! Ensemble size
  INTEGER, INTENT(in) :: dim_cv_ens_p       ! Number of columns in HV_p
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens) ! PE-local ensemble
  REAL, INTENT(inout) :: v_p(dim_cv_ens_p)  ! PE-local control vector
  REAL, INTENT(inout) :: Vv_p(dim_p)        ! PE-local result vector
!EOP

! *** local variables ***
  INTEGER :: i, member, row          ! Counters
  REAL :: fact                       ! Scaling factor
  REAL, ALLOCATABLE :: Vmat_p(:,:)   ! Extended ensemble perturbation matrix
  REAL, ALLOCATABLE :: v_g(:)        ! Global control vector
  REAL :: invdimens                  ! Inverse ensemble size


! *********************
! *** Compute V v_p ***
! *********************
  
  ! *** Generate control vector transform matrix ***

  ALLOCATE(Vmat_p(dim_p, dim_ens))
  
  Vv_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        Vv_p(row) = Vv_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO
  
  DO member = 1, dim_ens
     Vmat_p(:,member) = ens_p(:,member) - Vv_p(:)
  END DO

  ! Fill additional columns (if Vmat_p holds multiple sets of localized ensenbles)
  DO i = 2, mcols_cvec_ens
     DO member = (i-1)*dim_ens+1, i*dim_ens
        Vmat_p(:,member) = ens_p(:,member-(i-1)*dim_ens)
     END DO
  END DO


  fact = 1.0/SQRT(REAL(dim_cvec_ens-1))

  IF (type_opt/=3) THEN

     ! Transform control variable to state increment
     CALL dgemv('n', dim_p, dim_cv_ens_p, fact, Vmat_p, &
          dim_p, v_p, 1, 0.0, Vv_p, 1)

  ELSE

     ! Gather global control vector
     ALLOCATE(v_g(dim_cvec_ens))
  
     CALL MPI_AllGatherV(v_p, dim_cv_ens_p, MPI_REAL8, &
          v_g, dims_cv_p, off_cv_p, MPI_REAL8, &
          COMM_filter, MPIerr)

     ! Transform control variable to state increment
     CALL dgemv('n', dim_p, dim_cvec_ens, fact, Vmat_p, &
          dim_p, v_g, 1, 0.0, Vv_p, 1)

     DEALLOCATE(v_g)

  END IF

  DEALLOCATE(Vmat_p)

END SUBROUTINE cvt_ens_pdaf
