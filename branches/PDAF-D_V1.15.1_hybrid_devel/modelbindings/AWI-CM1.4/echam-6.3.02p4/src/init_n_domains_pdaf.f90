!$Id: init_n_domains_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
!
! !INTERFACE:
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to set the number of local analysis 
! domains for the PE-local domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_n_domains)
! Called by: PDAF_lestkf_update  (as U_init_n_domains)
! Called by: PDAF_letkf_update   (as U_init_n_domains)
!EOP


! ************************************
! *** Initialize number of domains ***
! ************************************

  if (mype_filter==0) write (*,*) 'ECHAM-PDAF TEMPLATE init_n_domains_pdaf'

  n_domains_p = 0

END SUBROUTINE init_n_domains_pdaf
