!$Id$
!BOP
!
! !ROUTINE: get_obs_f_pdaf --- get vector of synthetic observations from PDAF
!
! !INTERFACE:
SUBROUTINE get_obs_f_pdaf(step, dim_obs, observation)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! The routine is called when synthetic observations
! are generated with PDAF. With the call, the user
! is provided with a generated observation vetor. 
! This can then e.g. be written to a file.
!
! The routine is called by all filter processes.
!
! Version for the dummy model with domain 
! decomposition. Here, the state is fully observed.
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE obs_gp_pdafomi, &
       ONLY: file_syntobs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL, INTENT(out)   :: observation(dim_obs) ! Observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs   (as U_get_obs_f)
!EOP


! *************************
! *** store observation ***
! *************************

  CALL write_syn_obs(step, file_syntobs, dim_obs, observation, 1)

END SUBROUTINE get_obs_f_pdaf

