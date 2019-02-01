!$Id$
!BOP
!
! !ROUTINE: get_obs_full --- get vector of dynthetic observations from PDAF
!
! !INTERFACE:
SUBROUTINE get_obs_full(step, dim_obs, observation)

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
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL, INTENT(out)   :: observation(dim_obs) ! Observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs   (as U_get_obs_f)
!EOP

! *** local variables ***


! *************************
! *** store observation ***
! *************************

  CALL write_syn_obs(step, 'twinobs.nc',dim_obs, observation, 0)

END SUBROUTINE get_obs_full

