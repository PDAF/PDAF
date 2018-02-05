!$Id: init_obs_l_pdaf.F90 1857 2017-12-14 18:26:01Z lnerger $
!BOP
!
! !ROUTINE: init_obs_l_pdaf --- Initialize local observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain in 
! PDAF\_lseik\_analysis.  It has to initialize 
! the local vector of observations for the 
! current local analysis domain.
!
! Implementation for the 2D online example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_f, obs_index_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain index
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_init_obs_l)
! Called by: PDAF_lestkf_analysis  (as U_init_obs_l)
! Called by: PDAF_letkf_analysis   (as U_init_obs_l)
!EOP


! *** local variables ***
  INTEGER :: i          ! counter


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  DO i = 1, dim_obs_l
     observation_l(i) = obs_f(obs_index_l(i))
  END DO

END SUBROUTINE init_obs_l_pdaf

