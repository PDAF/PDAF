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
     &     write_ens, nfields, off_fields_p, dim_fields_p, dim_fields
      USE mod_parallel_pdaf, ONLY: mype_filter, npes_filter,
     &     COMM_filter, MPI_REAL8, MPIerr, MPIstatus, MPI_SUM

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
      _RL, ALLOCATABLE :: var_est_p(:)   ! estimated process-local total variance
      _RL, ALLOCATABLE :: rmse_est_g(:)  ! estimated RMS error
      _RL, ALLOCATABLE :: variance_p(:)  ! model state variances
      CHARACTER(len=3) :: ensstr          ! String for ensemble member
      CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
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
      ALLOCATE(var_est_p(nfields))
      ALLOCATE(rmse_est_g(nfields))
      
C     Initialize numbers
      var_est_p  = 0.0
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

      
C     ************************************************************
C     *** Compute RMS errors according to sampled covar matrix ***
C     ************************************************************

C     Process-local variances
      DO ifield = 1, nfields
         DO i = off_fields_p(ifield)+1,
     &        dim_fields_p(ifield)+off_fields_p(ifield)
            var_est_p(ifield) = var_est_p(ifield) + variance_p(i)
         ENDDO
      ENDDO

      DEALLOCATE(variance_p)

C     Compute global total variance using MPI_Allreduce
      CALL MPI_Allreduce(var_est_p, rmse_est_g, nfields, MPI_REAL8,
     &      MPI_SUM, COMM_filter, MPIerr)
      
      DO ifield = 1, nfields
         rmse_est_g(ifield)
     &        = SQRT(rmse_est_g(ifield) / REAL(dim_fields(ifield), 8))
      ENDDO
         
       
C     *********************
C     *** Screen Output ***
C     *********************

C     Output RMS errors given by sampled covar matrix
      IF (mype_filter == 0) THEN
        WRITE (*, '(a, 12x, a/, a,10x, 5(a,7x)/, a, 5x, 5es12.4, a, a)') 
     &        'PDAFuser:','RMS errors according to sampled variance: ',
     &        'PDAFuser:', 'etaN ', 'theta', 'salt ', 'uvel ', 'vvel ',
     &        'PDAFuser:', rmse_est_g, '  RMSe-',anastr
        write (*,*) 'RMSe_est: ', rmse_est_g
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

      DEALLOCATE(var_est_p, rmse_est_g)

      END
