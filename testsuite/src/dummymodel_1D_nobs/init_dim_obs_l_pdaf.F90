!$Id: init_dim_obs_l_pdaf.F90 1249 2012-01-27 13:44:54Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analsis domain.
!
! The routine is called by each filter process.
!
! Implementation for the dummy model with domain 
! decomposition. In this variant a local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated. (See also the variant 
! using a global observation domain)
!
! !REVISION HISTORY:
! 2007-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: local_range, dim_obs, local_dims_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain     ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
!EOP


! *** local variables ***
  INTEGER :: i                       ! Counter
  INTEGER :: domain_g                ! Global domain index


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! local analysis grid point  
  dim_obs_l = 1 + 2*local_range

  IF (dim_obs_f<dim_obs_l) THEN
     dim_obs_l = dim_obs_f
  END IF


!   ! Get domain index in global grid
!   domain_g = domain
!   DO i = 1, mype_filter
!      domain_g = domain_g + local_dims_obs(i)
!   ENDDO
! 
!   ! size left sided
!   IF (domain_g > local_range) THEN
!      dim_obs_l = dim_obs_l + local_range
!   ELSE
!      dim_obs_l = dim_obs_l + domain_g - 1
!   ENDIF
! 
!   ! size right sided
!   IF (domain_g + local_range <= dim_obs) THEN
!      dim_obs_l = dim_obs_l + local_range
!   ELSE
!      dim_obs_l = dim_obs_l + dim_obs - domain_g
!   ENDIF

END SUBROUTINE init_dim_obs_l_pdaf

