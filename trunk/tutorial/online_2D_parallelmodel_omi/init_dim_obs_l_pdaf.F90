!$Id$
!>  Set dimension of local observation vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over
!! all local analysis domains. It has to set 
!! the dimension of the local observation vector 
!! for the current local analysis domain.
!! 
!! For multiple observation types the order of the
!! calls over the observation types is critical.
!! It has to be consystent with those in obs_obs_f_pdaf.
!!
!! Implementation for the 2D online example
!! with parallelization.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

  USE mod_model, &             ! Model variables
       ONLY: ny, nx_p
  USE mod_parallel_pdaf, &     ! assimilation parallelization variables
       ONLY: mype_filter
  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: init_dim_obs_l_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: i                       ! Counters
  REAL :: coords_l(2)                ! Coordinates of local analysis domain
  INTEGER :: off_p                   ! Process-local offset in global state vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Global coordinates of local analysis domain
  ! We use grid point indices as coordinates, but could e.g. use meters
  off_p = 0
  DO i = 1, mype_filter
     off_p = off_p + nx_p*ny
  END DO
  coords_l(1) = REAL(CEILING(REAL(domain_p+off_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p+off_p) - (coords_l(1)-1)*REAL(ny)

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_dim_obs_l_pdafomi(domain_p, step, coords_l, dim_obs_f, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdaf

