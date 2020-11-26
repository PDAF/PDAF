!$Id: init_dim_obs_f_pdaf.F90 1861 2017-12-19 07:38:48Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_lseik\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to determine the dimension of the 
! observation vector according to the current 
! time step for all observations required for 
! the analyses in the loop over all local 
! analysis domains on the PE-local state domain.
!
! Implementation for the 2D online example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: dim_obs_p, obs_f, obs_index_p, coords_obs_f
  USE mod_model, &
       ONLY: nx, ny, nx_p
  USE mod_parallel_pdaf, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
! Called by: PDAF_letkf_update   (as U_init_dim_obs)
!EOP

! *** Local variables
  INTEGER :: i, j                     ! Counters
  INTEGER :: cnt0, cnt_p, cnt0_p      ! Counters
  INTEGER :: status=0                 ! Status flag
  INTEGER :: off_p                    ! Process-local offset in state vector
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file
  REAL, ALLOCATABLE :: obs_p(:)       ! Process-local observation vector
  REAL, ALLOCATABLE :: coords_obs_p(:,:) ! Coordinates of process-local observations
  CHARACTER(len=2) :: stepstr         ! String for time step


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Determine offset in state vector for this process
  off_p = 0
  DO i= 1, mype_filter
     off_p = off_p + nx_p*ny
  END DO

  ! Read observation field form file
  ALLOCATE(obs_field(ny, nx))

  IF (step<10) THEN
     WRITE (stepstr, '(i1)') step
  ELSE
     WRITE (stepstr, '(i2)') step
  END IF

  OPEN (12, file='../../inputs_online/obs_step'//TRIM(stepstr)//'.txt', status='old')
  DO i = 1, ny
     READ (12, *) obs_field(i, :)
  END DO
  CLOSE (12)

  ! Count process-local observations
  cnt0 = 0
  cnt_p = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (cnt0 > off_p .AND. cnt0 <= off_p+nx_p*ny) THEN
           IF (obs_field(i,j) > -999.0) cnt_p = cnt_p + 1
        END IF
     END DO
  END DO

  ! Set number of local observations
  dim_obs_p = cnt_p

  ! Initialize local vector and coordinates of observations and index array
  ! The coordinates are grid point indices here, but could e.g. be meters
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))
  ALLOCATE(coords_obs_p(2, dim_obs_p))

  cnt0 = 0
  cnt_p = 0
  cnt0_p = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (cnt0 > off_p .AND. cnt0 <= off_p + nx_p*ny) THEN
           cnt0_p = cnt0_p + 1
           IF (obs_field(i,j) > -999.0) THEN
              cnt_p = cnt_p + 1
              obs_index_p(cnt_p) = cnt0_p
              obs_p(cnt_p) = obs_field(i, j)
              coords_obs_p(1, cnt_p) = REAL(j)
              coords_obs_p(2, cnt_p) = REAL(i)
           END IF
        END IF
     END DO
  END DO


! *** Gather full observation information ***

  ! Full observation dimension
  CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

  IF (ALLOCATED(obs_f)) DEALLOCATE(obs_f)
  IF (ALLOCATED(coords_obs_f)) DEALLOCATE(coords_obs_f)
  ALLOCATE(obs_f(dim_obs_f))
  ALLOCATE(coords_obs_f(2, dim_obs_f))

  ! Get full observation vector
  CALL PDAF_gather_obs_f(obs_p, obs_f, status)

  ! Get full array of coordinates
  CALL PDAF_gather_obs_f2(coords_obs_p, coords_obs_f, 2, status)


! *** Clean up ***

  DEALLOCATE(obs_field, obs_p, coords_obs_p)

END SUBROUTINE init_dim_obs_f_pdaf

