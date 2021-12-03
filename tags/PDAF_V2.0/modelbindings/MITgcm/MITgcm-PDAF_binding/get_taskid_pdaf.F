C     $Id: init_parallel_pdaf.F 7 2018-01-31 18:11:38Z lnerger $
C     
C     !ROUTINE: get_taskid_pdaf --- Initialize communicators for PDAF

C     !INTERFACE:
      SUBROUTINE get_taskid_pdaf(mpi_task_id)

C     !DESCRIPTION:
C     Routine to return the value of the task ID.
C     This is particular for MITgcm sicne otherwise
C     it tries to write into the same scratch file
C     with sevaral tasks.
C
C     !REVISION HISTORY:
C     2021-03 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_parallel_pdaf,
     &      ONLY: task_id

      IMPLICIT NONE    
  
C     !ARGUMENTS:
      INTEGER, INTENT(out)   :: mpi_task_id

C     Set task ID of MITgcm to task_id 
C     initialized by init_parallel_pdaf
      mpi_task_id = task_id
      
      END SUBROUTINE get_taskid_pdaf
