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
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE PDAFlocal, &             ! Routine to provide local indices to PDAF
       ONLY: PDAFlocal_set_indices
  USE mod_model, &             ! Model variables
       ONLY: ny, nx_p
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l
  USE mod_parallel_pdaf, &     ! assimilation parallelization variables
       ONLY: mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER :: i                       ! Counters
  INTEGER :: off_p                   ! Process-local offset in global state vector
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) !< Indices of local state vector in PE-local global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  dim_l = 1


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Global coordinates of local analysis domain
  ! We use grid point indices as coordinates, but could e.g. use meters
  off_p = 0
  DO i = 1, mype_filter
     off_p = off_p + nx_p*ny
  END DO
  coords_l(1) = REAL(CEILING(REAL(domain_p+off_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p+off_p) - (coords_l(1)-1)*REAL(ny)


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! Here the local domain is a single grid point and variable given by DOMAIN_P
  id_lstate_in_pstate(1) = domain_p

  ! Provide the index vector to PDAF
  CALL PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  DEALLOCATE(id_lstate_in_pstate)

END SUBROUTINE init_dim_l_pdaf
