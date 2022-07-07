!$Id: mod_model.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides shared variables for the 3D dummy model.
!
! !REVISION HISTORY:
! 2008-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:  
!    ! Control model run - available as command line options
  INTEGER :: dims(3)             ! Model state dimensions
  INTEGER :: step_null           ! Initial time step of assimilation

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: dim_l(3)            ! Array for dimensions of PE-local domain
  INTEGER :: step_final          ! Final time step
  REAL    :: dt                  ! Time step size
  REAL, ALLOCATABLE :: field(:,:,:)        ! Array holding model field
  INTEGER, ALLOCATABLE :: dims_l_all(:,:)  ! Array for local state dimensions
!EOP

END MODULE mod_model
