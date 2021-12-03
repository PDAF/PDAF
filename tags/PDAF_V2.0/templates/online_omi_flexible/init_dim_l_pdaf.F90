!$Id$
!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during analysis step
!! in PDAF_X_update in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model  state on the current analysis
!! domain.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l, id_lstate_in_pstate

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_l_pdaf.F90: Set local state dimension here!'

!   dim_l = ??


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Global coordinates of local analysis domain

!  coords_l = ??


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

!  id_lstate_in_pstate = ??

END SUBROUTINE init_dim_l_pdaf
