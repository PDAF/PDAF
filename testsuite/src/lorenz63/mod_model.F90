!$Id: mod_model.F90 784 2009-12-07 10:30:03Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides shared variables for the dummy model.
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
  INTEGER :: step_null           ! Initial time step of assimilation

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: step_final          ! Final time step
  REAL    :: dt                  ! Time step size
  REAL, ALLOCATABLE :: x(:)      ! Array holding model field
  REAL    :: gamma               ! Model parameter
  REAL    :: rho                 ! Model parameter
  REAL    :: beta                ! Model parameter
!EOP

END MODULE mod_model
