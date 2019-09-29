!$Id$
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
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, init_dim_obs_f_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, init_dim_obs_f_B

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
  INTEGER :: dim_obs_f_A, dim_obs_f_B ! Observation dimensions
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file
  CHARACTER(len=2) :: stepstr         ! String for time step


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_A = 0
  dim_obs_f_B = 0

  ! Call observation specific routines
  ! The routines are independent, so it is not important in
  ! which order they are called
  IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
  IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)

  dim_obs_f = dim_obs_f_A + dim_obs_f_B


END SUBROUTINE init_dim_obs_f_pdaf

