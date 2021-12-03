!$Id: init_dim_obs_pdaf.F90 1366 2013-04-24 16:25:05Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! Implementation for the 2D offline example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY : nx, ny, local_dims, obs_p, obs_index_p
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
!EOP

! *** Local variables
  INTEGER :: i, j                     ! Counters
  INTEGER :: cnt, cnt0, cnt_p, cnt0_p ! Counters
  INTEGER :: off_p                    ! process-local offset in state vector
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file


! ****************************************
! *** Initialize observation dimension ***
! ****************************************

  ! Determine offset in state vector for this process
  off_p = 0
  DO i= 1, mype_filter
     off_p = off_p + local_dims(i)
  END DO

  ! Read observation field form file
  ALLOCATE(obs_field(ny, nx))

  OPEN (12, file='../../inputs_offline/obs.txt', status='old')
  DO i = 1, ny
     READ (12, *) obs_field(i, :)
  END DO
  CLOSE (12)

  ! Count observations
  cnt = 0
  cnt0 = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (cnt0 > off_p .AND. cnt0 <= off_p+local_dims(mype_filter+1)) THEN
           IF (obs_field(i,j) > -999.0) cnt = cnt + 1
        END IF
     END DO
  END DO

  ! Set number of observations
  dim_obs_p = cnt

  ! Initialize vector of observations and index array
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))

  cnt0 = 0
  cnt_p = 0
  cnt0_p = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (cnt0 > off_p .AND. cnt0 <= off_p+local_dims(mype_filter+1)) THEN
           cnt0_p = cnt0_p + 1
           IF (obs_field(i,j) > -999.0) THEN
              cnt_p = cnt_p + 1
              obs_index_p(cnt_p) = cnt0_p
              obs_p(cnt_p) = obs_field(i, j)
           END IF
        END IF
     END DO
  END DO


! *** Clean up ***

  DEALLOCATE(obs_field)

END SUBROUTINE init_dim_obs_pdaf

