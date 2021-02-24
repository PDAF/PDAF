!$Id$
!>  Routine to set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!!
!! The routine is called during analysis step
!! in the loop over all local analysis domains.
!! It has to set the dimension of local model 
!! state on the current analysis domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate, off_fields_p, coords_l, pi
  USE mo_geoloc, &
       ONLY: philat_2d, philon_2d
  USE mo_decomposition, &
       ONLY: dc=>local_decomposition

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** Local variables ***
  INTEGER :: i, j, k                           ! Counters
  INTEGER, ALLOCATABLE :: nod_local_domain(:)  ! Grid node indices for current atm column
  INTEGER :: nlat, nlon, nlev                  ! ATM levels


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  nlev = dc%nlev  
  nlat = dc%nglat
  nlon = dc%nglon  

  ! *** Local state dimension
  
  dim_l = 6*nlev + 1


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  k = 1
  blkloop: DO i = 1,dc%ngpblks
     DO j = 1,dc%nproma
        coords_l(1) = philon_2d(j,i)
        coords_l(2) = philat_2d(j,i)
        IF (k == domain_p) EXIT blkloop
        k = k + 1
     END DO
  END DO blkloop

  ! Convert to radians
  coords_l(1) = coords_l(1)  * pi / 180.0
  coords_l(2) = coords_l(2) * pi / 180.0


! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate arrays
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(nod_local_domain(nlev))
  ALLOCATE(id_lstate_in_pstate(dim_l))


  ! *** node indices for local domain

  DO i = 1, nlev
    nod_local_domain(i) = domain_p + (i-1) * nlat * nlon 
  END DO

  ! *** indices for full state vector ***
  ! Currently only one variable -3D temperature + uv- in the state vector

  ! T
  id_lstate_in_pstate(1:nlev) = nod_local_domain(1:nlev) + off_fields_p(1)

  ! lsp
  id_lstate_in_pstate(nlev+1) = nod_local_domain(1) + off_fields_p(2)  

  !vorticity
  id_lstate_in_pstate(nlev+2:2*nlev+1) = nod_local_domain(1:nlev) + off_fields_p(3)

  !divergence
  id_lstate_in_pstate(2*nlev+2:3*nlev+1) = nod_local_domain(1:nlev) + off_fields_p(4)

  !humidity
  id_lstate_in_pstate(3*nlev+2:4*nlev+1) = nod_local_domain(1:nlev) + off_fields_p(5)
  
  ! u
  id_lstate_in_pstate(4*nlev+2:5*nlev+1) = nod_local_domain(1:nlev) + off_fields_p(6)

  ! v
  id_lstate_in_pstate(5*nlev+2:6*nlev+1) = nod_local_domain(1:nlev) + off_fields_p(7)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(nod_local_domain)

END SUBROUTINE init_dim_l_pdaf
