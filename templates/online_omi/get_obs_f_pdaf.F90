!$Id$
!> Get vector of synthetic observations from PDAF
!!
!! User-supplied routine for PDAF.
!! The routine is called by PDAF_get_obs when
!! synthetic observations are generated with PDAF. 
!! With the call, the user is provided with a
!! generated observation vector. This can then e.g. 
!! be written to a file.
!!
!! The routine is called by all filter processes.
!!
!! This is a full implementation in which the state is 
!! just written into a file using the provided routine
!! write_syn_obs. It can be used without changes for the 
!! local filters LESTKF/LESTKF/LSEIK/LNETF if the full 
!! observation vector includes all available observations.
!!
!! __Revision history:__
!! * 2019-01 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!
SUBROUTINE get_obs_f_pdaf(step, dim_obs_f, observation_f)

  USE mod_assimilation, &     ! Variables for assimilation
       ONLY: file_syntobs
  USE mod_parallel_pdaf, &    ! Parallelization
       ONLY: mype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Currrent time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector


! *********************************
! *** write observation to file ***
! *********************************

  IF (mype_filter==0) THEN
     CALL write_syn_obs(step, file_syntobs, dim_obs_f, observation_f, 1)
  END IF

END SUBROUTINE get_obs_f_pdaf

