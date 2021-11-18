C     $Id$
C     !BOP
C
C     !ROUTINE: next_observation_pdaf --- Initialize information on next observation
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The subroutine is called before each forecast phase
C     by PDAF\_get\_state. It has to initialize the number
C     of time steps until the next available observation
C     (nsteps) and the current model time (time). In
C     addition the exit flag (exit) has to be initialized.
C     It indicates if the data assimilation process is
C     completed such that the ensemble loop in the model
C     routine can be exited.
C
C     The routine is called by all processes.
C
C     !REVISION HISTORY:
C     2013-09 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
#include "CPP_OPTIONS.h"
C     !INTERFACE:
      SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
C     !USES:
      USE mod_assimilation,
     &     ONLY: myTime, nsteps_between_obs

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: stepnow   ! Number of the current time step
      INTEGER, INTENT(out) :: nsteps    ! Number of time steps until next obs
      INTEGER, INTENT(out) :: doexit    ! Whether to exit forecasting (1 for exit)
      REAL(kind=8), INTENT(out) :: time ! Current model (physical) time

C     !CALLING SEQUENCE:
C     Called by: PDAF_get_state   (as U_next_obs)
!EOP


C     *****************************************************
C     *** Set number of time steps between observations ***
C     *****************************************************

      nsteps = nsteps_between_obs


C     *** Set exit flag and store time ***
C     (both not used in fully-parallel implementation variant)

      doexit = 0
      time = myTime

      END
