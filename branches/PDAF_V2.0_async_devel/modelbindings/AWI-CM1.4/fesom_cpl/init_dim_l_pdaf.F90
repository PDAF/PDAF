!$Id: init_dim_l_pdaf.F90 2153 2020-03-05 16:18:50Z lnerger $
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
       ONLY: id_lstate_in_pstate, off_fields_p, coords_l
  USE o_mesh, &
       ONLY: num_layers_below_nod2d, nod3d_below_nod2d, coord_nod2D
  USE g_rotate_grid, &
       ONLY: r2g

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** Local variables ***
  INTEGER :: nlay                              ! Number of layers for current domain
  INTEGER, ALLOCATABLE :: nod_local_domain(:)  ! Grid node indices for current water column


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  nlay = num_layers_below_nod2d(domain_p) + 1

  ! *** Local state dimension
  dim_l = 5*nlay + 1


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Get location of current water column (basis point)
  CALL r2g(coords_l(1), coords_l(2), coord_nod2d(1, domain_p), coord_nod2d(2, domain_p))


! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate arrays
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(nod_local_domain(nlay))
  ALLOCATE(id_lstate_in_pstate(dim_l))


  ! *** node indices for local domain
  nod_local_domain(1:nlay) = nod3d_below_nod2d(1:nlay, domain_p)

  ! *** indices for full state vector ***

  ! SSH
  id_lstate_in_pstate(1) = nod_local_domain(1) + off_fields_p(1)

  ! U
  id_lstate_in_pstate(1+1:nlay+1) = &
       nod_local_domain(1:nlay) + off_fields_p(2)

  ! V
  id_lstate_in_pstate(nlay+1+1 : 2*nlay+1) = &
       nod_local_domain(1:nlay) + off_fields_p(3)

  ! W
  id_lstate_in_pstate(2*nlay+1+1 : 3*nlay+1) = &
       nod_local_domain(1:nlay) + off_fields_p(4)

  ! Temperature
  id_lstate_in_pstate(3*nlay+1+1 : 4*nlay+1) = &
       nod_local_domain(1:nlay) + off_fields_p(5)

  ! Salinity
  id_lstate_in_pstate(4*nlay+1+1 : 5*nlay+1) = &
       nod_local_domain(1:nlay) + off_fields_p(6)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(nod_local_domain)

END SUBROUTINE init_dim_l_pdaf
