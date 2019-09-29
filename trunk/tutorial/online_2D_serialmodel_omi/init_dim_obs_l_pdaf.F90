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
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range
   USE mod_model, &
        ONLY: nx, ny
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, init_dim_obs_l_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, init_dim_obs_l_B

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
  REAL :: coords_l(2)                 ! Coordinates of local analysis domain
  INTEGER :: dim_obs_l_A, dim_obs_l_B ! Dimension of each observatiob
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Coordinates of local analysis domain 
  ! We use grid point indices as coordinates, but could e.g. use meters
  coords_l(1) = REAL(CEILING(REAL(domain_p)/REAL(ny)))
  coords_l(2) = REAL(domain_p) - (coords_l(1)-1)*REAL(ny)

  ! Initialize offsets with zero
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  IF (assim_A) CALL init_dim_obs_l_A(coords_l, local_range, dim_obs_l_A, offset_obs_l, offset_obs_f)
  IF (assim_B) CALL init_dim_obs_l_B(coords_l, local_range, dim_obs_l_B, offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_A + dim_obs_l_B

END SUBROUTINE init_dim_obs_l_pdaf

