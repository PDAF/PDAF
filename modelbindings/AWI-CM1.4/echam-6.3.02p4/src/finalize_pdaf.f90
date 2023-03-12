!BOP
!
! !ROUTINE: finalize_pdaf - Finalize ensemble assimilation with PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()

! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
!
! !REVISION HISTORY:
! 2017-07 - L. Nerger - initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_world, MPIerr, MPI_COMM_WORLD, mype_submodel
  
  IMPLICIT NONE
  
!EOP

  IF (mype_submodel == 0) THEN
     ! Show timings for PDAF
     CALL PDAF_print_info(3)

     ! Show allocated memory for PDAF
     CALL PDAF_print_info(2)

     WRITE (*,'(5x,a)') 'ECHAM-PDAF: Assimilation with PDAF completed!'
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIerr)


END SUBROUTINE finalize_pdaf
