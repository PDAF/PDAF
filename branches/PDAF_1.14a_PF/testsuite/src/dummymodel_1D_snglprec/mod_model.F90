!$Id: mod_model.F90 1711 2016-12-20 10:17:58Z lnerger $
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
  INTEGER :: dim_state           ! Model state dimension
  INTEGER :: step_null           ! Initial time step of assimilation

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
  INTEGER :: step_final          ! Final time step
  REAL    :: dt                  ! Time step size
  REAL, ALLOCATABLE :: field(:)          ! Array holding model field
  INTEGER, ALLOCATABLE :: local_dims(:)  ! Array for local state dimensions

  ! Declare the observation array here, for flexibility
  REAL, ALLOCATABLE :: observation(:)    ! Array of observations
!EOP

END MODULE mod_model
