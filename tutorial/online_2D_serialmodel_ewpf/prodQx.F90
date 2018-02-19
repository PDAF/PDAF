!$Id: obs_op_pdaf.F90 3 2013-09-05 10:28:51Z lnerger $
!BOP
!
! !ROUTINE: prodQx --- product of Q with state vector
!
! !INTERFACE:
SUBROUTINE prodQx(dim, x, model_state, Qx)

! !DESCRIPTION:
! User-supplied routine for PDAF.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner
! Later revisions - see svn log


  IMPLICIT NONE

  INTEGER, INTENT(in) :: dim            ! State dimension
  REAL, INTENT(in) :: x(dim)            ! State to multiply with
  REAL, INTENT(in) :: model_state(dim)  ! current model state (for linearization)
  REAL, INTENT(OUT) :: Qx(dim)          ! output model state

! Local variables
  INTEGER :: i
  REAL, ALLOCATABLE :: Qmat(:,:)


! Initialize covariance matrix Q
  ALLOCATE(Qmat(dim, dim))  ! Usualy we would avoid allocating the full matrix Q
  
  Qmat = 0.0
  DO i=1, dim
     Qmat(i,i) = 1.0
  END DO

! Apply Q to x

! Qx = ....
  Qx = x   ! TEMPORARY


  DEALLOCATE(Qmat)

END SUBROUTINE prodQx






