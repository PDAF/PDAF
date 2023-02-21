!$Id: init_dim_l_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
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
!! with parallelization.
!!
!! 2013-02 - Lars Nerger - Initial code
!! Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: ny, local_dims, coords_l, id_lstate_in_pstate
  USE mod_parallel, &          ! Parallelization
       ONLY: mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER :: i                     ! Counters
  INTEGER :: off_p                 ! Process-local offset in global state vector


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
     off_p = off_p + local_dims(i)
  END DO
  coords_l(1) = REAL(CEILING(REAL(domain_p+off_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p+off_p) - (coords_l(1)-1)*REAL(ny)


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! Here the local domain is a single grid point and variable given by DOMAIN_P
  id_lstate_in_pstate(1) = domain_p

END SUBROUTINE init_dim_l_pdaf
