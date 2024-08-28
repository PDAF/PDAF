!$Id: init_dim_l_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the loop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! Implementation for the 2D offline example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAFlocal, &             ! Routine to provide local indices to PDAF
       ONLY: PDAFlocal_set_indices
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: ny, local_dims, coords_l
  USE mod_parallel, &          ! assimilation parallelization variables
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! !LOCAL VARIABLES:
  INTEGER :: i                       ! Counters
  INTEGER :: off_p                   ! Process-local offset in global state vector
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
!EOP


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
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! Here the local domain is a single grid point and variable given by DOMAIN_P
  id_lstate_in_pstate(1) = domain_p

  ! Provide the index vector to PDAF
  CALL PDAFlocal_set_indices(dim_l, id_lstate_in_pstate)

  ! Deallocate index array
  DEALLOCATE(id_lstate_in_pstate)

END SUBROUTINE init_dim_l_pdaf
