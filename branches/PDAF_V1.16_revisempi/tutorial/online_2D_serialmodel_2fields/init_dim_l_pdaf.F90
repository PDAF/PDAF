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
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! 2013-02 - Lars Nerger - Initial code
!! Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_model, &             ! Model variables
       ONLY: ny
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: coords_l, id_lstate_in_pstate, n_fields, off_fields

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER :: i                     ! Counter


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  dim_l = n_fields


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! We use grid point indices as coordinates, but could e.g. use meters
  coords_l(1) = REAL(CEILING(REAL(domain_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p) - (coords_l(1)-1)*REAL(ny)


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! Here the local domain is a single grid point
  ! The variables given by DOMAIN_P + offsets
  DO i=1, n_fields
     id_lstate_in_pstate(i) = domain_p + off_fields(i)
  END DO

END SUBROUTINE init_dim_l_pdaf
