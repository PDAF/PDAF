!$Id: finalize_pdaf.F90 1415 2013-09-25 14:33:26Z lnerger $
!>  Finalize PDAF
!!
!! This routine calls routines for output on
!! timing and memory use, to deallocate
!! PDAF_internal arrays, and to finalize MPI.
!!
!! __Revision history__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE finalize_pdaf()

  USE PDAF_interfaces_module, &   ! PDAF interface definitions
       ONLY: PDAF_print_info, PDAF_deallocate
  USE mod_parallel, &             ! Parallelization
       ONLY: mype_world

  IMPLICIT NONE    


! *** Show allocated memory for PDAF ***
  IF (mype_world==0) CALL PDAF_print_info(2)

! *** Print PDAF timings onto screen ***
  IF (mype_world==0) CALL PDAF_print_info(3)

! *** Deallocate PDAF arrays ***
  CALL PDAF_deallocate()

END SUBROUTINE finalize_pdaf
