!$Id: finalize_pdaf.F90 1857 2017-12-14 18:26:01Z lnerger $
!BOP
!
! !ROUTINE: finalize_pdaf --- Finalize PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()

! !DESCRIPTION:
! This routine call MPI_finalize
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_model, &
       ONLY: mype_world

  IMPLICIT NONE    
  
! !CALLING SEQUENCE:
! Called by: main program
!EOP

! *** Show allocated memory for PDAF ***
  IF (mype_world==0) CALL PDAF_print_info(2)

! *** Print PDAF timings onto screen ***
  IF (mype_world==0) CALL PDAF_print_info(3)

! *** Deallocate PDAF arrays
  CALL PDAF_deallocate()

END SUBROUTINE finalize_pdaf
