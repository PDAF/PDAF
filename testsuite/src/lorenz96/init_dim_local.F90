!$Id: init_dim_local.F90 836 2010-01-29 17:19:57Z lnerger $
!BOP
!
! !ROUTINE: init_dim_local --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_local(step, domain, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during analysis step
! in the llop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! This variant is for the Lorenz96 model without
! parallelization. We simply consider each single 
! grid point  as a local analysis domain.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step    ! Current time step
  INTEGER, INTENT(in)  :: domain  ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l   ! Local state dimension

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
!EOP


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  ! Simply one here
  dim_l = 1

END SUBROUTINE init_dim_local
