!$Id: mod_modeltime.F90 746 2009-08-04 12:16:28Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_modeltime

! !DESCRIPTION:
! This module provides variables for the time information
! of a model.
!
! !REVISION HISTORY:
! 2009-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:  
!    ! Control model run - available as command line options
  INTEGER :: total_steps         ! Total number of time steps to perform

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time                ! model time
!EOP

END MODULE mod_modeltime
