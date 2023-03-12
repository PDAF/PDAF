!$Id: init_dim_obs_l_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
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
       ONLY: dim_state, local_dims
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: local_range

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
  dim_obs_l = 1

  ! Get domain index in global grid
  domain_g = domain
  DO i = 1, mype_filter
     domain_g = domain_g + local_dims(i)
  ENDDO

  ! size left sided
  IF (domain_g > INT(local_range)) THEN
     dim_obs_l = dim_obs_l + INT(local_range)
  ELSE
     dim_obs_l = dim_obs_l + domain_g - 1
  ENDIF

  ! size right sided
  IF (domain_g + INT(local_range) <= dim_state) THEN
     dim_obs_l = dim_obs_l + INT(local_range)
  ELSE
     dim_obs_l = dim_obs_l + dim_state - domain_g
  ENDIF

END SUBROUTINE init_dim_obs_l_pdaf

