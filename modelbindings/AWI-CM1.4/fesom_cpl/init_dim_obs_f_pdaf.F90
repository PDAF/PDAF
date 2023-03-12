!$Id: init_dim_obs_f_pdaf.F90 2293 2020-05-11 14:52:41Z lnerger $
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
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, &
       ONLY: assim_sst, write_en4data

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
! Called by: PDAF_letkf_update   (as U_init_dim_obs)
!EOP


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize counting of observations
  dim_obs_f = 0

  IF (assim_sst) THEN
     ! Assimilate SST observations
     CALL init_dim_obs_f_sst_pdaf(step, dim_obs_f)
  END IF

  IF (write_en4data) THEN
     ! Generate file holdig EN4 profile data on FESOM mesh
     CALL write_profiles_pdaf(step, dim_obs_f)
  END IF

END SUBROUTINE init_dim_obs_f_pdaf

