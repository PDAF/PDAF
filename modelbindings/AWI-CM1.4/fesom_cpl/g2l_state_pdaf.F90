!$Id: g2l_state_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, &
       ONLY: index_local_domain
!!DEBUG
USE mod_parallel_pdaf, &
       ONLY: mype_filter


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
!EOP
  
! *** Local variables *** 
  INTEGER :: i, i_full


! *************************************
! *** Initialize local state vector ***
! *************************************

  DO i = 1, dim_l
     i_full = index_local_domain(i)
     state_l(i) = state_p(i_full)
     IF(state_l(i) > 50.0) WRITE(*,*) 'warning_extreme_state_l_g2l',state_l(i),mype_filter,domain_p
  ENDDO

END SUBROUTINE g2l_state_pdaf
