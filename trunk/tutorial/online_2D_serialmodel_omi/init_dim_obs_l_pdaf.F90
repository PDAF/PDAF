!$Id$
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analysis domain.
! 
! For multiple observation types the order of the
! calls over the observation types is critical.
! It has to be consystent with those in obs_obs_f_pdaf.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code for PDAF_OMI
! Later revisions - see repository log
!
! !USES:
  USE mod_model, &
       ONLY: ny
  USE interface_pdafomi, &
       ONLY: init_dim_obs_l_pdafomi

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
!EOP


! *** local variables ***
  REAL :: coords_l(2)                    ! Coordinates of local analysis domain


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Coordinates of local analysis domain 
  ! We use grid point indices as coordinates, but could e.g. use meters
  coords_l(1) = REAL(CEILING(REAL(domain_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p) - (coords_l(1)-1)*REAL(ny)

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL init_dim_obs_l_pdafomi(domain_p, step, coords_l, dim_obs_f, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdaf

