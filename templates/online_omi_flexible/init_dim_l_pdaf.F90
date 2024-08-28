!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF/LKNETF
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

  USE PDAFlocal, &             ! Routine to provide local indices to PDAF
       ONLY: PDAFlocal_set_indices
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_l_pdaf.F90: Set local state dimension here!'

!  dim_l = ??


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
  ALLOCATE(id_lstate_in_pstate(dim_l))

!  id_lstate_in_pstate = ??

  ! Provide the index vector to PDAF
  CALL PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  DEALLOCATE(id_lstate_in_pstate)

END SUBROUTINE init_dim_l_pdaf
