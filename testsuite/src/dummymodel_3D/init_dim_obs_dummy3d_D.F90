!$Id: init_dim_obs_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK/SEIK/EnKF):
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! The routine is called by all filter processes.
!
! For the domain-decomposed dummy model, the full
! model state is observed. Thus, the number of
! observations equals the PE-local state dimension.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter
  USE mod_model, &
       ONLY: dims_l_all
  USE mod_assimilation, &
        ONLY: n_obs, local_dims_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! Local variables
  INTEGER :: i   ! Counter
  INTEGER, SAVE :: allocflag = 0


! *** Initialize observation dimension ***

  ! Array of local observation domains
  IF (allocflag == 0) THEN
     allocate(local_dims_obs(npes_filter))
     allocflag = 1
  END IF

  DO i = 1, npes_filter
     local_dims_obs(i) = dims_l_all(1, i) * dims_l_all(2, i)
  END DO

  ! dimension for local domain
  dim_obs_p = local_dims_obs(mype_filter + 1)
  
END SUBROUTINE init_dim_obs

