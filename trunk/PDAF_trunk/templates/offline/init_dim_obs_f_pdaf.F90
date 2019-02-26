!$Id$
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_lseik\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to determine the dimension of the 
! observation vector according to the current 
! time step for all observations required for 
! the analyses in the loop over all local 
! analysis domains on the PE-local state domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY : nx, ny, local_dims, dim_obs_p, &
!        obs_f, obs_index_p, coords_obs_f
!   USE mod_parallel, &
!        ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
! Called by: PDAF_letkf_update   (as U_init_dim_obs)
! Called by: PDAF_lnetf_update   (as U_init_dim_l)
!EOP

! *** Local variables


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_obs_f_pdaf.F90: Set full observation dimension here!'

! dim_obs_f = ?

END SUBROUTINE init_dim_obs_f_pdaf

