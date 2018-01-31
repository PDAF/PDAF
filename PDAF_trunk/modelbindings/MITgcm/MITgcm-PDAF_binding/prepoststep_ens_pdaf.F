C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
C
C     !INTERFACE:
      SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, 
     &     dim_obs_p, state_p, Uinv, ens_p, flag)

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
      _RL,  INTENT(inout) :: state_p(dim_p)              ! PE-local forecast/analysis state
      ! The array 'state_p' is not generally not initialized in the case of SEIK.
      ! It can be used freely here.
      _RL,  INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1)  ! Inverse of matrix U
      _RL,  INTENT(inout) :: ens_p(dim_p, dim_ens)       ! PE-local state ensemble
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
      _RL :: invdim_ens                  ! Inverse ensemble size
      _RL :: invdim_ensm1                ! Inverse of ensemble size minus 1
      _RL, ALLOCATABLE :: rmserror_est(:)                ! estimated RMS error
      _RL, ALLOCATABLE :: variance_p(:)  ! model state variances
      CHARACTER(len=3) :: ensstr          ! String for ensemble member
      CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
!!$  ! Variables for parallelization - global fields
!!$  INTEGER :: offset                   ! Row-offset according to domain decomposition
!!$  REAL, ALLOCATABLE :: variance(:)    ! local variance
!!$  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
!!$  REAL, ALLOCATABLE :: state(:)       ! global state vector
!!$  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
!!$  REAL,ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain

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
!!$
!!$
!!$! ******************************************************
!!$! *** Assemble global variance vector on filter PE 0 ***
!!$! ******************************************************
!!$
!!$  WRITE (*,*) 'TEMPLATE prepoststep_ens_pdaf.F90: Initialize variance, either directly or with MPI'
!!$
!!$!   PE0_a: IF (mype_filter /= 0) THEN
!!$! 
!!$!      ! send sub-fields from PEs /=0
!!$!      CALL MPI_send(variance_p(1 : dim_p), dim_p, &
!!$!           MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)
!!$! 
!!$!   ELSE PE0_a
!!$!      ! receive and assemble variance field
!!$! 
!!$!      ! On PE 0 init variance directly
!!$!      variance(1 : dim_p) = variance_p(1 : dim_p)
!!$! 
!!$!      ! Receive part of variance field from PEs > 0 into 
!!$!      ! correct part of global variance
!!$! 
!!$!      offset = 0
!!$! 
!!$!      DO i = 2, npes_filter
!!$!         ! Increment offset
!!$!         offset = offset + nx_p*ny
!!$! 
!!$!         ! Receive variance part
!!$!         CALL MPI_recv(variance(1 + offset), nx_p*ny, &
!!$!              MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
!!$!      END DO
!!$!       
!!$!   END IF PE0_a
!!$
!!$  DEALLOCATE(variance_p)

      
C     ************************************************************
C     *** Compute RMS errors according to sampled covar matrix ***
C     ************************************************************

      DO ifield = 1, nfields
         DO i = off_fields(ifield)+1,
     &        dim_fields(ifield)+off_fields(ifield)
            rmserror_est(ifield) = rmserror_est(ifield) + variance_p(i)
         ENDDO
         rmserror_est(ifield)
     &        = SQRT(rmserror_est(ifield) / REAL(dim_fields(ifield), 8))
      ENDDO

      DEALLOCATE(variance_p)

       
C     *********************
C     *** Screen Output ***
C     *********************

C     Output RMS errors given by sampled covar matrix
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
