C     $Id$
C
C     !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The routine is called for global filters (e.g. SEIK)
C     before the analysis and after the ensemble transformation.
C     For local filters (e.g. LSEIK) the routine is called
C     before and after the loop over all local analysis
C     domains.
C     The routine provides full access to the state
C     estimate and the state ensemble to the user.
C     Thus, user-controlled pre- and poststep
C     operations can be performed here. For example
C     the forecast and the analysis states and ensemble
C     covariance matrix can be analyzed, e.g. by
C     computing the estimated variances.
C     For the offline mode, this routine is the place
C     in which the writing of the analysis ensemble
C     can be performed.
C
C     If a user considers to perform adjustments to the
C     estimates (e.g. for balances), this routine is
C     the right place for it.
C
C     !REVISION HISTORY:
C     2013-09 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
#include "CPP_OPTIONS.h"
C     !INTERFACE:
      SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p,
     &     dim_obs_p, state_p, Uinv, ens_p, flag)
C     !USES:
      USE mod_assimilation,
     &     ONLY: dim_state, myIter, myThid, iter_null,
     &     nfields, off_fields, dim_fields, write_ens
      USE mod_parallel_pdaf, ONLY: mype_filter, npes_filter,
     &     COMM_filter, MPI_DOUBLE_PRECISION, MPIerr, MPIstatus

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step         ! Current time step (negative for call after forecast)
      INTEGER, INTENT(in) :: dim_p        ! PE-local state dimension
      INTEGER, INTENT(in) :: dim_ens      ! Size of state ensemble
      INTEGER, INTENT(in) :: dim_ens_p    ! PE-local size of ensemble
      INTEGER, INTENT(in) :: dim_obs_p    ! PE-local dimension of observation vector
      REAL(kind=8),  INTENT(inout) :: state_p(dim_p)              ! PE-local forecast/analysis state
      ! The array 'state_p' is not generally not initialized in the case of SEIK.
      ! It can be used freely here.
      REAL(kind=8),  INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1)  ! Inverse of matrix U
      REAL(kind=8),  INTENT(inout) :: ens_p(dim_p, dim_ens)       ! PE-local state ensemble
      INTEGER, INTENT(in) :: flag         ! PDAF status flag

C     !CALLING SEQUENCE:
C     Called by: PDAF_get_state      (as U_prepoststep)
C     Called by: PDAF_X_update       (as U_prepoststep)
C     Calls: MPI_send
C     Calls: MPI_recv

C     *** local variables ***

# include "SIZE.h"
# include "EEPARAMS.h"
# include "DYNVARS.h"
      INTEGER :: i, j, member, domain    ! Counters
      INTEGER :: ifield                  ! Counter
      INTEGER, SAVE :: allocflag = 0     ! Flag for memory counting
      REAL(kind=8) :: invdim_ens                   ! Inverse ensemble size
      REAL(kind=8) :: invdim_ensm1                 ! Inverse of ensemble size minus 1
      REAL(kind=8), ALLOCATABLE :: rmserror_est(:) ! estimated RMS error
      REAL(kind=8), ALLOCATABLE :: variance_p(:)   ! model state variances
      CHARACTER(len=3) :: ensstr          ! String for ensemble member
      CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
C     !$  ! Variables for parallelization - global fields
C     !$  INTEGER :: offset                   ! Row-offset according to domain decomposition
C     !$  REAL(kind=8), ALLOCATABLE :: variance(:)    ! local variance
C     !$  REAL(kind=8), ALLOCATABLE :: ens(:,:)       ! global ensemble
C     !$  REAL(kind=8), ALLOCATABLE :: state(:)       ! global state vector
C     !$  REAL(kind=8),ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
C     !$  REAL(kind=8),ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain

      CHARACTER*(MAX_LEN_MBUF) suff
      CHARACTER*(MAX_LEN_MBUF) fcstanls
      INTEGER :: koffset, bj, bi


C     **********************
C     *** INITIALIZATION ***
C     **********************

      IF (mype_filter == 0) THEN
        IF (step-iter_null==0) THEN
           WRITE (*, '(a, 8x, a)') 'PDAFuser: ',
     &          'Analize initial state ensemble'
           anastr = 'ini'
        ELSE
           IF (step<0) THEN
              WRITE (*, '(a, 8x, a)') 'PDAFuser: ',
     &             'Analize and write forecasted state ensemble'
              anastr = 'for'
           ELSE
              WRITE (*, '(a, 8x, a)') 'PDAFuser: ',
     &             'Analize and write assimilated state ensemble'
              anastr = 'ana'
           END IF
        END IF
      END IF

C     Allocate fields
      ALLOCATE(variance_p(dim_p))
      ALLOCATE(rmserror_est(nfields))

C     Initialize numbers
      rmserror_est  = 0.0
      invdim_ens    = 1.0d0 / REAL(dim_ens, 8)
      invdim_ensm1  = 1.0d0 / REAL(dim_ens - 1, 8)


C     ********************************************
C     *** Compute estimated variances          ***
C     ********************************************

C     *** Compute mean state ***
      IF (mype_filter == 0) WRITE (*, '(a, 8x, a)')
     &  'PDAFuser: ','--- compute ensemble mean'

      state_p = 0.0
      DO member = 1, dim_ens
        DO i = 1, dim_p
          state_p(i) = state_p(i) + ens_p(i, member)
        END DO
      END DO
      state_p(:) = invdim_ens * state_p(:)


C     *** Compute sampled variances ***

      variance_p(:) = 0.0
      DO member = 1, dim_ens
        DO j = 1, dim_p
          variance_p(j) = variance_p(j)
     &           + (ens_p(j, member) - state_p(j))
     &           * (ens_p(j, member) - state_p(j))
        END DO
      END DO
      variance_p(:) = invdim_ensm1 * variance_p(:)
C     
C     
C     ******************************************************
C     *** Assemble global variance vector on filter PE 0 ***
C     ******************************************************
C     
C        WRITE (*,*) 'TEMPLATE prepoststep_ens_pdaf.F90: Initialize variance, either directly or with MPI'
C     
C        PE0_a: IF (mype_filter /= 0) THEN
C     
C           ! send sub-fields from PEs /=0
C           CALL MPI_send(variance_p(1 : dim_p), dim_p, &
C                MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)
C     
C        ELSE PE0_a
C           ! receive and assemble variance field
C     
C           ! On PE 0 init variance directly
C           variance(1 : dim_p) = variance_p(1 : dim_p)
C     
C           ! Receive part of variance field from PEs > 0 into
C           ! correct part of global variance
C     
C           offset = 0
C     
C           DO i = 2, npes_filter
C              ! Increment offset
C              offset = offset + nx_p*ny
C     
C              ! Receive variance part
C              CALL MPI_recv(variance(1 + offset), nx_p*ny, &
C                   MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
C           END DO
C     
C        END IF PE0_a
C     
C        DEALLOCATE(variance_p)


C     !************************************************************
C     !*** Compute RMS errors according to sampled covar matrix ***
C     !************************************************************

      DO ifield = 1, nfields
         DO i = off_fields(ifield)+1,
     &        dim_fields(ifield)+off_fields(ifield)
            rmserror_est(ifield) = rmserror_est(ifield) + variance_p(i)
         ENDDO
         rmserror_est(ifield)
     &        = SQRT(rmserror_est(ifield) / REAL(dim_fields(ifield), 8))
      ENDDO
      
      WRITE(suff,'(I10.10)') myIter
      CALL WRITE_FLD_XY_RL('ENS.',suff,variance_p,myIter,myThid)

      DEALLOCATE(variance_p)


C     *********************
C     *** Screen Output ***
C     *********************

C     !Output RMS errors given by sampled covar matrix
      IF (mype_filter == 0) THEN
        WRITE (*, '(a, 12x, a/, a,10x, 5(a,7x)/, a, 5x, 5es12.4, a, a)')
     &        'PDAFuser:','RMS errors according to sampled variance: ',
     &        'PDAFuser:', 'etaN ', 'theta', 'salt ', 'uvel ', 'vvel ',
     &        'PDAFuser:', rmserror_est, '  RMSe-',anastr
        write (*,*) 'RMSe_est: ', rmserror_est
      END IF


C     *******************
C     *** File output ***
C     *******************

C     map back to geometric/geographic fields
      CALL distribute_state_pdaf(dim_p, state_p)

C     write mean state
      WRITE(suff,'(I10.10)') myIter
      IF ( step-iter_null == 0 ) THEN
       fcstanls = "_initial."
      ELSE
         IF ( step > 0 ) THEN
C     save analysis
            fcstanls = "_analysis."
         ELSE
C     save forcast
            fcstanls = "_forecast."
         ENDIF
      ENDIF
      CALL WRITE_FLD_XY_RL('ETAN'//fcstanls,suff,etaN,myIter,myThid)
      CALL WRITE_FLD_XYZ_RL('THETA'//fcstanls,suff,theta,myIter,myThid)
      CALL WRITE_FLD_XYZ_RL('SALT'//fcstanls,suff,salt,myIter,myThid)
      CALL WRITE_FLD_XYZ_RL('UVEL'//fcstanls,suff,uvel,myIter,myThid)
      CALL WRITE_FLD_XYZ_RL('VVEL'//fcstanls,suff,vvel,myIter,myThid)

C     Write ensemble
      IF (write_ens) THEN

         WRITE(suff,'(I10.10)') myIter

         DO i=1,dim_ens
C           map back to geometric/geographic fields
            CALL distribute_state_pdaf(dim_p, ens_p(:,i))

            WRITE(ensstr,'(I3.3)') i
            CALL WRITE_FLD_XY_RL('ETAN_'//ensstr//fcstanls,
     &           suff,etaN,myIter,myThid)
            CALL WRITE_FLD_XYZ_RL('THETA_'//ensstr//fcstanls,
     &           suff,theta,myIter,myThid)
            CALL WRITE_FLD_XYZ_RL('SALT_'//ensstr//fcstanls,
     &           suff,salt,myIter,myThid)
            CALL WRITE_FLD_XYZ_RL('UVEL_'//ensstr//fcstanls,
     &           suff,uvel,myIter,myThid)
            CALL WRITE_FLD_XYZ_RL('VVEL_'//ensstr//fcstanls,
     &           suff,vvel,myIter,myThid)

         ENDDO
      ENDIF

      
C     ********************
C     *** finishing up ***
C     ********************

      DEALLOCATE(rmserror_est)
      
      END
