!$Id: mod_model.F90 61 2019-02-01 08:49:36Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides shared variables for the Lorenz05b model.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:
!    ! Control model run - available as command line options
  INTEGER :: dim_state           ! Model state dimension
  INTEGER :: step_null           ! Initial time step of assimilation

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: step_final          ! Final time step
  REAL    :: dt                  ! Time step size
  REAL, ALLOCATABLE :: x(:)      ! Array holding model field
  REAL    :: forcing             ! Model parameter
  INTEGER :: k_avg               ! Model parameter
!EOP

END MODULE mod_model
