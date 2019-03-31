!$Id: init_dim_l_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the llop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! The routine is called by each filter process.
!
! Version for the dummy model with domain 
! decomposition. We simply consider each single 
! grid point  as a local analysis domain. In 3D 
! this could be a water column.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
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

END SUBROUTINE init_dim_l_pdaf
