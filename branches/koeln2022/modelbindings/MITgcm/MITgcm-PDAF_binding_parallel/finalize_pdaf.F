C     $Id$
C
C     !ROUTINE: finalize_pdaf - Initialize ensemble assimilation with PDAF
C
C     !INTERFACE:
      SUBROUTINE finalize_pdaf()

C     This routine collects the initialization of variables for PDAF.
C     In addition, the initialization routine PDAF_init is called
C     such that the internal initialization of PDAF is performed.
C     
C     !REVISION HISTORY:
C     2017-12 - L. Nerger - Initial version
C     Later revisions - see svn log
C
C     !USES:
      USE mod_parallel_model, ONLY : mype_world

      IMPLICIT NONE
  
C     !CALLING SEQUENCE:
C     Called by: the_main_loop


      IF (mype_world == 0) THEN
         ! Show timings for PDAF
         CALL PDAF_print_info(3)

         ! Show allocated memory for PDAF
         CALL PDAF_print_info(2)
      END IF

      WRITE (*,'(a,5x,a)')
     &     'PDAFuser:', 'Assimilation with PDAF completed!'


      call PDAF_deallocate()
  
      END SUBROUTINE finalize_pdaf
