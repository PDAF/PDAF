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
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
        ONLY: coords_l, distance_l, id_lobs_in_fobs
!        local_range, coords_obs_f
!   USE mod_parallel, &
!        ONLY: mype_filter

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
! Called by: PDAF_lnetf_update   (as U_init_dim_obs_l)
!EOP

! *** local variables ***


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_obs_l_pdaf.F90: Set local observation dimension here!'

! dim_obs_l = ?


! *************************************
! *** Initialize array of distances ***
! *************************************

  ! Allocate array
  IF (ALLOCATED(distance_l)) DEALLOCATE(distance_l)
  ALLOCATE(distance_l(dim_obs_l))

!  distance_l = ?


! ********************************************************
! *** Initialize array of indices of local observation ***
! *** vector elements in full observation vector       ***
! ********************************************************

  ! Allocate array
  IF (ALLOCATED(id_lobs_in_fobs)) DEALLOCATE(id_lobs_in_fobs)
  ALLOCATE(id_lobs_in_fobs(dim_obs_l))

!  id_lobs_in_fobs = ?

END SUBROUTINE init_dim_obs_l_pdaf

