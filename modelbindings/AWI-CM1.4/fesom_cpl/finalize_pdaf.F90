!$Id: $
!>  Routine to finalize ensemble assimilation with PDAF
!!
!! This routine Calls routines that perform screen output
!! of PDAF's internal timing and memory information.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE finalize_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_print_info, PDAF_deallocate
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_submodel, MPIerr, MPI_COMM_WORLD
  
  IMPLICIT NONE

  IF (mype_submodel == 0) THEN

     ! Show timings for PDAF
     CALL PDAF_print_info(3)

     ! Show allocated memory for PDAF
     CALL PDAF_print_info(2)

     WRITE (*,'(5x,a)') 'FESOM-PDAF: Assimilation with PDAF completed!'
  END IF

  ! *** Deallocate PDAF arrays ***
  CALL PDAF_deallocate()

  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIerr)

END SUBROUTINE finalize_pdaf
