C     $Id$
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
C!BOP
C
C     !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
C
C     !INTERFACE:
      SUBROUTINE init_pdaf( nIter0, myTimeLoc, myIterLoc, myThidLoc )

C     !DESCRIPTION:
C     This routine collects the initialization of variables for PDAF.
C     In addition, the initialization routine PDAF_init is called
C     such that the internal initialization of PDAF is performed.
C     This variant is for the online mode of PDAF.
C
C     This routine is generic. However, it assumes a constant observation
C     error (rms_obs_etan). Further, with parallelization the local state
C     dimension dim_state_p is used.
C
C     !REVISION HISTORY:
C     2008-10 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
C     Interface definitions for PDAF
      USE pdaf_interfaces_module,
     &     ONLY: PDAF_init
C     Parallelization variables for model
      USE mod_parallel_model,
     &     ONLY: mype_world, COMM_model, abort_parallel
C     Parallelization variables fro assimilation
      USE mod_parallel_pdaf, ONLY: n_modeltasks, task_id, COMM_filter,
     &     COMM_couple, filterpe
C     Variables for assimilation
      USE mod_assimilation,  ONLY: dim_state_p, dim_state, screen,
     &     filtertype, subtype, dim_ens, iter_null,
     &     nsteps_between_obs, type_forget, forget,
     &     locweight, local_range, srange, type_trans, type_sqrt,
     &     myIter, myThid, MyTime, write_ens, 
     &     n_fields, id, dim_fields, off_fields
      USE obs_etan_pdafomi,
     &      ONLY: assim_etan, rms_obs_etan

      IMPLICIT NONE

#include "SIZE.h"

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: nIter0
      INTEGER, INTENT(in) :: myIterLoc
      INTEGER, INTENT(in) :: myThidLoc
      REAL(kind=8), INTENT(in) :: myTimeLoc

C     !CALLING SEQUENCE:
C     Called by: the_main_loop
C     Calls: init_pdaf_parse
C     Calls: init_pdaf_info
C     Calls: PDAF_init
C     Calls: PDAF_get_state
C!EOP

C     Local variables
      INTEGER :: i                       ! Counter
      INTEGER :: filter_param_i(7)       ! Integer parameter array for filter
      REAL(kind=8) :: filter_param_r(2)  ! Real parameter array for filter
      INTEGER :: status_pdaf             ! PDAF status flag
      INTEGER :: doexit, steps           ! Not used in this implementation
      REAL(kind=8) :: timenow            ! Not used in this implementation

C     External subroutines
      EXTERNAL :: init_ens_pdaf    ! Ensemble initialization
      EXTERNAL :: next_observation_pdaf, ! Provide time step, model time,
                                   ! and dimension of next observation
     &     distribute_state_pdaf,  ! Distribute a state vector to model fields
     &     prepoststep_ens_pdaf    ! User supplied pre/poststep routine

C     Namelist for PDAF settings
      NAMELIST /pdaf_nml/ screen, filtertype, subtype,
     &     iter_null, nsteps_between_obs,
     &     type_forget, forget, type_trans, type_sqrt,
     &     locweight, local_range, srange, write_ens,
     &     assim_etan, rms_obs_etan


C     *** Store information from MITgcm to be used in call-back routines ***
      myIter = myIterLoc
      myThid = myThidLoc
      myTime = myTimeLoc


C     ***************************
C     ***   Initialize PDAF   ***
C     ***************************

      IF (mype_world == 0) THEN
         WRITE (*,'(/a,1x,a)') 'PDAFuser:',
     &        'INITIALIZE PDAF - ONLINE MODE'
      END IF

      ! Store initial iteration
      iter_null = myIter


C     *****************************
C     ***   Define state vector ***
C     *****************************

C     Define state dimension
      dim_state_p = sNx*sNy + 4*sNx*sNy*Nr  ! Process local dimension
      dim_state   = Nx*Ny + 4*Nx*Ny*Nr      ! Global dimension

C     Number of model fields in state vector
      n_fields = 5  ! Number of fields in state vector     

C     Specify field IDs in state vector ('id' allows to give the field a name)
      id%etan = 1
      id%theta= 2
      id%salt = 3
      id%uvel = 4
      id%vvel = 5

C     Define field dimensions
      allocate(dim_fields(n_fields))

      dim_fields(id%etan) = Nx*Ny    ! EtaN
      dim_fields(id%theta)= Nx*Ny*Nr ! Theta
      dim_fields(id%salt) = Nx*Ny*Nr ! Salt
      dim_fields(id%uvel) = Nx*Ny*Nr ! uVel
      dim_fields(id%vvel) = Nx*Ny*Nr ! vVel

C     Offsets of fields in state vector
      allocate(off_fields(n_fields))

      off_fields(1) = 0
      DO i = 2, n_fields
         off_fields(i) = off_fields(i-1) + dim_fields(i-1)
      END DO


C     **********************************************************
C     ***   CONTROL OF PDAF - used in call to PDAF_init      ***
C     **********************************************************

C     *** IO options ***
      screen = 2    ! Write screen output (1) for output, (2) add timings

C     *** Filter specific variables
      filtertype = 7            ! Type of filter
                    !   SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                    !   ESTKF (6), LESTKF (7), LEnKF (8), NETF (9), LNETF (10), PF (12)
      subtype = 0               ! subtype of filter:
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   EnKF:
                          !     (0) analysis for large observation dimension
                          !     (1) analysis for small observation dimension
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) LETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   ESTKF:
                          !     (0) standard ESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LESTKF:
                          !     (0) standard LESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   LEnKF:
                          !     (0) Standard form of EnKF with covariance localization
                          !   NETF:
                          !     (0) standard NETF 
                          !   LNETF:
                          !     (0) standard LNETF
                          !   PF:
                          !     (0) standard PF with resampling
      type_trans = 0            ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ESTKF/LESTKF:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! NETF/LNETF:
                           ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           ! (1) use identity transformation
      type_forget = 0           ! Type of forgetting factor in ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
      forget  = 1.0             ! Forgetting factor
      type_sqrt = 0             ! Type of transform matrix square-root
                                !   (0) symmetric square root, (1) Cholesky decomposition


C     *********************************************************************
C     ***   Settings for analysis steps  - used in call-back routines   ***
C     *********************************************************************

C     Forecast length
      nsteps_between_obs = 10   ! Number of time steps between analysis/assimilation steps

C     Whether to write full ensemble files
      write_ens = .false.

C     Observation error standard deviation
      rms_obs_etan = 0.10       ! For sea surface height etan

C     Localization settings for LESTKF/LETKF
      locweight = 4             ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
      local_range = 50000.0     ! Localization radius
      srange = local_range      ! Support radius for 5th-order polynomial
                    ! or range for 1/e for exponential weighting


C     ***********************************
C     *** Some optional functionality ***
C     ***********************************

C     *** Read namelist file for PDAF ***

      OPEN (500,file='pdaf.nml')
      READ (500,NML=pdaf_nml)
      CLOSE (500)

C     *** Parse command line options   ***
C     *** This is optional, but useful ***

      call init_pdaf_parse()


C     *** Initial Screen output ***
C     *** This is optional      ***

      IF (mype_world == 0) call init_pdaf_info()


C     *****************************************************
C     *** Call PDAF initialization routine on all PEs.  ***
C     ***                                               ***
C     *** For all filters, first the arrays of integer  ***
C     *** and real number parameters are initialized.   ***
C     *** Subsequently, PDAF_init is called.            ***
C     *****************************************************

      filter_param_i(1) = dim_state_p ! State dimension
      filter_param_i(2) = dim_ens     ! Size of ensemble
      filter_param_i(3) = 0           ! Smoother lag (not implemented here)
      filter_param_i(4) = 0           ! Incremental analysis (not implemented here)
      filter_param_i(5) = type_forget ! Type of forgetting factor
      filter_param_i(6) = type_trans  ! Type of ensemble transformation
      filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
      filter_param_r(1) = forget      ! Forgetting factor

      CALL PDAF_init(filtertype, subtype, nIter0,
     &     filter_param_i, 2,
     &     filter_param_r, 2,
     &     COMM_model, COMM_filter, COMM_couple,
     &     task_id, n_modeltasks, filterpe, init_ens_pdaf,
     &     screen, status_pdaf)


C     *** Check whether initialization of PDAF was successful ***
      IF (status_pdaf /= 0) THEN
         WRITE (*,'(/1x,a6,i3,a43,i4,a1/)')
     &        'ERROR ', status_pdaf,
     &        ' in initialization of PDAF ',
     &        '- stopping           ! (PE ', mype_world,')'
         CALL abort_parallel()
      END IF


C     **********************************
C     *** Prepare ensemble forecasts ***
C     **********************************

      CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf,
     &     distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

      END
