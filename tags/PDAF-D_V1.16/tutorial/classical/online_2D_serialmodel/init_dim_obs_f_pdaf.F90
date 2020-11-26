!$Id: init_dim_obs_f_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
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
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY : obs_f, obs_index_p, coords_obs_f
  USE mod_model, &
       ONLY: nx, ny

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
  INTEGER :: cnt, cnt0                ! Counters
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file
  CHARACTER(len=2) :: stepstr         ! String for time step


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

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

  ! Count observations
  cnt = 0
  DO j = 1, nx
     DO i= 1, ny
        IF (obs_field(i,j) > -999.0) THEN
           cnt = cnt + 1
        END IF
     END DO
  END DO

  ! Set number of observations
  dim_obs_f = cnt

  ! Initialize vector of observations and index array
  ! The coordinates are grid point indices here, but could e.g. be meters
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_f)) DEALLOCATE(obs_f)
  IF (ALLOCATED(coords_obs_f)) DEALLOCATE(coords_obs_f)
  ALLOCATE(obs_index_p(dim_obs_f))
  ALLOCATE(obs_f(dim_obs_f))
  ALLOCATE(coords_obs_f(2, dim_obs_f))

  cnt = 0
  cnt0 = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (obs_field(i,j) > -999.0) THEN
           cnt = cnt + 1
           obs_index_p(cnt) = cnt0
           obs_f(cnt) = obs_field(i, j)
           coords_obs_f(1, cnt) = REAL(j)
           coords_obs_f(2, cnt) = REAL(i)
        END IF
     END DO
  END DO


! *** Clean up ***

  DEALLOCATE(obs_field)

END SUBROUTINE init_dim_obs_f_pdaf

