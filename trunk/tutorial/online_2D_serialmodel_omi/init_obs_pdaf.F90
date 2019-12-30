!$Id$
!>  Initialize observation vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF/NETF
!!
!! The routine is called during the analysis step. 
!! It has to provide the PE-local observation vector 
!! for the current time step.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE init_obs_pdaf(step, dim_obs_p, observation_p)

  USE mod_assimilation, &     ! Assimilation variables
       ONLY: obs_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p        !< PE-local dimension of obs. vector
  REAL, INTENT(out)   :: observation_p(dim_obs_p) !< PE-local observation vector


! ******************************
! *** Initialize observation ***
! ******************************
  
  observation_p = obs_p

END SUBROUTINE init_obs_pdaf

