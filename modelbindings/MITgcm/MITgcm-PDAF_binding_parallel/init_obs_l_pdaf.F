C     $Id$
C
C     !ROUTINE: init_obs_l_pdaf --- Initialize local observation vector
C
C     !INTERFACE:
      SUBROUTINE init_obs_l_pdaf(domain_p, step,
     &      dim_obs_l, observation_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during the analysis step
C     on each local analysis domain in 
C     PDAF\_lseik\_analysis.  It has to initialize 
C     the local vector of observations for the 
C     current local analysis domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code 
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs_f, obs_index_l
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain index
      INTEGER, INTENT(in) :: step       ! Current time step
      INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
      REAL*8, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_analysis  (as U_init_obs_l)
C     Called by: PDAF_letkf_analysis   (as U_init_obs_l)

C     *** local variables ***
      INTEGER :: i              ! counter


C     *******************************************
C     *** Initialize local observation vector ***
C     *******************************************

      ! This is generic because we just use obs_index_l
      
      DO i = 1, dim_obs_l
         observation_l(i) = obs_f(obs_index_l(i))
      END DO

      END SUBROUTINE init_obs_l_pdaf

