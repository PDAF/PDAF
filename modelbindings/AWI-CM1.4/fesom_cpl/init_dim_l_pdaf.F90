!$Id: init_dim_l_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
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
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, &
       ONLY: index_local_domain, offset
  USE o_mesh, &
       ONLY: num_layers_below_nod2d, nod3d_below_nod2d

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
!EOP

! *** Local variables ***
  integer :: layer                      ! Counter
  INTEGER :: nlay                       ! Num of layers for current domain
  INTEGER, ALLOCATABLE :: nod_local_domain(:)


! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  ! Number of layers
  nlay = num_layers_below_nod2d(domain_p) + 1

  ! Local state dimension
  dim_l = 5*nlay + 1


! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate arrays
  IF (ALLOCATED(index_local_domain)) DEALLOCATE(index_local_domain)
  ALLOCATE(index_local_domain(dim_l))
  ALLOCATE(nod_local_domain(nlay))


  ! *** node indices for local domain
  nod_local_domain(1:nlay) = nod3d_below_nod2d(1:nlay, domain_p)

  ! *** indices for full state vector ***

  ! SSH
  index_local_domain(1) = nod_local_domain(1) + offset(1)

  ! U
  index_local_domain(1+1:nlay+1) = &
       nod_local_domain(1:nlay) + offset(2)

  ! V
  index_local_domain(nlay+1+1 : 2*nlay+1) = &
       nod_local_domain(1:nlay) + offset(3)

  ! W
  index_local_domain(2*nlay+1+1 : 3*nlay+1) = &
       nod_local_domain(1:nlay) + offset(4)

  ! Temperature
  index_local_domain(3*nlay+1+1 : 4*nlay+1) = &
       nod_local_domain(1:nlay) + offset(5)

  ! Salinity
  index_local_domain(4*nlay+1+1 : 5*nlay+1) = &
       nod_local_domain(1:nlay) + offset(6)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(nod_local_domain)

END SUBROUTINE init_dim_l_pdaf
