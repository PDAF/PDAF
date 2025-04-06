!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!! 
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_state, do_omi_obsstats
  USE PDAF, &                  ! PDAF and PDAF-OMI diagnostic routines
       ONLY: PDAF_diag_stddev, PDAFomi_diag_obs_rmsd, PDAFomi_diag_stats

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** local variables ***
  INTEGER :: i, j, member, domain     ! Counters
  INTEGER :: pdaf_status              ! status flag
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: ens_stddev                  ! ensemble STDDEV = estimated RMS error
  INTEGER :: nobs                     ! Number of observations in diagnostics
  REAL, POINTER :: obsRMSD(:)         ! Array of observation RMS deviations
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  REAL, POINTER :: obsstats(:,:)      ! Array of observation statistics
  ! Variables for parallelization - global fields
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
  REAL,ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prepoststep_pdaf.F90: Implement functionality of prepoststep '

  dim_state = dim_p ! FOR TESTING - valid without domain decomposition

  IF (mype_filter == 0) THEN
     IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analyze forecasted state ensemble'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
     END IF
  END IF


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  CALL PDAF_diag_stddev(dim_p, dim_ens, state_p, ens_p, &
        ens_stddev, 1, COMM_filter, pdaf_status)

 
! **************************************
! *** Compute observation statistics ***
! **************************************

!TEMPLATE: We include two statistics here, which are optional and partly redundant
  IF (do_omi_obsstats) THEN
     ! Compute RMS deviation between observation and observed ensemble mean
     CALL PDAFomi_diag_obs_rmsd(nobs, obsrmsd, 1/(mype_filter+1))

     ! Compute statistics on deviation between observation and observed ensemble
     CALL PDAFomi_diag_stats(nobs, obsstats, 1/(mype_filter+1))
  END IF


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
     WRITE (*, '(12x, a, es12.4)') &
          'RMS error according to sampled standard deviation: ', ens_stddev
  END IF


! *******************
! *** File output ***
! *******************

  ! Here, one could implement file output


! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_pdaf
