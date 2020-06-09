! Copyright (c) 2004-2019 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAFomi_obs_f.F90 333 2019-12-31 16:19:13Z lnerger $

!> PDAF-OMI routines for full observations
!!
!! This module contains subroutines to handle full observations. 
!! Further, it contains routines to restrict the global full vector of observations
!! to those observations that are relevant for a process-local model subdomain.
!! The routines are
!!
!! * PDAFomi_gather_obs_f \n
!!        Gather full observation information
!! * PDAFomi_gather_obsstate_f \n
!!        Gather a full observed state vector (used in observation operators)
!! * PDAFomi_init_obs_f \n
!!        Initialize full vector of observations for adaptive forgetting factor
!! * PDAFomi_init_obsvar_f \n
!!        Compute mean observation error variance for adaptive forgetting factor
!! * PDAFomi_prodRinvA \n
!!        Multiply an intermediate matrix of the global filter analysis
!!        with the inverse of the observation error covariance matrix
!! * PDAFomi_likelihood \n
!!        Compute likelihood for an ensemble member
!! * PDAFomi_add_obs_err \n
!!        Add observation error to some matrix
!! * PDAFomi_init_obscovar \n
!!        Initialize global observation error covariance matrix
!! * PDAFomi_deallocate_obs \n
!!        Deallocate arrays in observation type
!! * PDAFomi_set_domain_limits \n
!!        Set min/max coordinate locations of decomposed grid
!! * PDAFomi_get_domain_limits_unstr \n
!!        Find min/max coordinate locations in unstructured grid
!! * PDAFomi_get_local_ids_obs_f \n
!!        Find observations inside or close to process domain
!! * PDAFomi_limit_obs_f \n
!!        Reduce full observation vector to part relevant for local process domain
!!
!! The routine PDAFomi_get_domain_limits_unstr assumed geographic coordinates in radians
!! and within the range -pi to +pi for longitude (- is westward) and -pi/2 to +pi/2 for
!! latitude.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
MODULE PDAFomi_obs_f

  USE PDAF_mod_filtermpi, &
       ONLY: mype_filter, COMM_FILTER, MPI_INTEGER, MPIerr, MPI_MIN, MPI_MAX

  IMPLICIT NONE
  SAVE

! *** Module internal variables
  INTEGER :: debug=0                    !< Debugging flag

  REAL, ALLOCATABLE :: domain_limits(:) !< Limiting coordinates (NSWE) for process domain
  REAL, PARAMETER :: r_earth=6.3675e6   !< Earth radius in meters
  REAL, PARAMETER :: pi=3.141592653589793   !< Pi

! *** Data type to define the full observations by internally shared variables of the module
  type obs_f
     ! ---- Mandatory variables to be set in init_dim_obs_f ----
     INTEGER :: doassim=0                 !< Whether to assimilate this observation type
     INTEGER :: disttype                  !< Type of distance computation to use for localization
     INTEGER :: ncoord                    !< Number of coordinates use for distance computation
     INTEGER :: dim_obs_p                 !< number of PE-local observations
     INTEGER :: dim_obs_f                 !< number of full observations
     INTEGER, ALLOCATABLE :: id_obs_p(:,:) !< indices of process-local observed field in state vector
     REAL, ALLOCATABLE :: obs_f(:)        !< Full observed field
     REAL, ALLOCATABLE :: ocoord_f(:,:)   !< Coordinates of full observation vector
     REAL, ALLOCATABLE :: ivar_obs_f(:)   !< Inverse variance of full observations
     ! ---- Variables with predefined values - they can be changed in init_dim_obs_f  ----
     INTEGER :: obs_err_type=0            !< Type of observation error: (0) Gauss, (1) Laplace
     ! ---- Optional variables - they can be set in init_dim_obs_f ----
     REAL, ALLOCATABLE :: icoeff_p(:,:)   !< Interpolation coefficients for obs. operator (optional)
     REAL, ALLOCATABLE :: domainsize(:)   !< Size of domain for periodicity (<=0 for no periodicity) (optional)
     ! ---- Optional variables set in obs_op_f when not using global full observation ---
     LOGICAL :: use_global_obs=.true.     !< Whether to use (T) global full obs. 
                                          !< or (F) obs. restricted to those relevant for a process domain
     INTEGER :: dim_obs_g                 !< global number of observations
     INTEGER, ALLOCATABLE :: id_obs_f_lim(:) !< Indices of domain-relevant full obs. in global vector of obs.
     ! ---- Mandatory variable to be set in obs_op_f ---
     INTEGER :: off_obs_f                 !< Offset of this observation in overall full obs. vector
  end type obs_f

!$OMP THREADPRIVATE(debug)


!-------------------------------------------------------------------------------
  
CONTAINS
!> Gather full observational information
!!
!! This routine uses different gather routines from PDAF to obtain
!! combined full observational information.
!!
!! The observation-type specific variables that are initialized here are
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!!
!! If the full observations are restricted to those relevant for a process domain, 
!! there are further initialized 
!! * thisobs\%dim_obs_g  - Number of global observations
!! * thisobs\%id_obs_f_lim - Ids of full observations in global observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_gather_obs_f(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
       ncoord, lradius, dim_obs_f)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs   !< Data type with full observation
    INTEGER, INTENT(in) :: dim_obs_p        !< Number of process-local observation
    REAL, INTENT(in) :: obs_p(:)            !< Vector of process-local observations
    REAL, INTENT(in) :: ivar_obs_p(:)       !< Vector of process-local inverse observation error variance
    REAL, INTENT(in) :: ocoord_p(:,:)       !< Array of process-local observation coordinates
    INTEGER, INTENT(in) :: ncoord           !< Number of rows of coordinate array
    REAL, INTENT(in) :: lradius             !< Localization radius (the maximum radius used in this process domain) 
    INTEGER, INTENT(out) :: dim_obs_f       !< Full number of observations

! *** Local variables ***
    INTEGER :: i                            ! Counter
    REAL, ALLOCATABLE :: obs_g(:)           ! Global full observation vector (used in case of limited obs.)
    REAL, ALLOCATABLE :: ivar_obs_g(:)      ! Global full inverse variances (used in case of limited obs.)
    REAL, ALLOCATABLE :: ocoord_g(:,:)      ! Global full observation coordinates (used in case of limited obs.)
    INTEGER :: status                       ! Status flag for PDAF gather operation
    INTEGER :: localfilter                  ! Whether the filter is domain-localized


! **************************************
! *** Gather full observation arrays ***
! **************************************

    ! Check  whether the filter is domain-localized
    CALL PDAF_get_localfilter(localfilter)


    lfilter: IF (localfilter==1) THEN

       ! For domain-localized filters: gather full observations

       fullobs: IF (thisobs%use_global_obs) THEN

          ! *** Use global full observations ***

          IF (mype_filter == 0) &
               WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use global full observations'

          ! *** Initialize global dimension of observation vector ***
          CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

          ! Store full and PE-local observation dimensions in module variables
          thisobs%dim_obs_p = dim_obs_p
          thisobs%dim_obs_f = dim_obs_f

          IF (mype_filter == 0) &
               WRITE (*, '(a, 8x, a, i7)') 'PDAFomi', &
               '--- Number of full observations ', dim_obs_f

          ! *** Gather full observation vector and corresponding coordinates ***

          ! Allocate full observation arrays
          ! The arrays are deallocated in deallocate_obs in this module
          ALLOCATE(thisobs%obs_f(dim_obs_f))
          ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
          ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))

          CALL PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, thisobs%obs_f, status)
          CALL PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, ivar_obs_p, thisobs%ivar_obs_f, status)
          CALL PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, ocoord_p, thisobs%ocoord_f, ncoord, status)

       ELSE fullobs

          ! *** Use full observations limited to those relevant for a process domain ***
          ! *** This can be more efficient as in the local analysis loop less        ***
          ! *** observations have a be checked for each analysis domain              ***

          IF (mype_filter == 0) &
               WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use limited full observations'

          ! *** Initialize global dimension of observation vector ***
          CALL PDAF_gather_dim_obs_f(dim_obs_p, thisobs%dim_obs_g)

          ! *** First gather global observation vector and corresponding coordinates ***

          ! Allocate global observation arrays
          ALLOCATE(obs_g(thisobs%dim_obs_g))
          ALLOCATE(ivar_obs_g(thisobs%dim_obs_g))
          ALLOCATE(ocoord_g(ncoord, thisobs%dim_obs_g))

          CALL PDAF_gather_obs_f_flex(dim_obs_p, thisobs%dim_obs_g, obs_p, obs_g, status)
          CALL PDAF_gather_obs_f_flex(dim_obs_p, thisobs%dim_obs_g, ivar_obs_p, ivar_obs_g, status)
          CALL PDAF_gather_obs_f2_flex(dim_obs_p, thisobs%dim_obs_g, ocoord_p, ocoord_g, ncoord, status)


          ! *** Now restrict the global observation arrays to the process-relevant parts ***

          ! Get number of full observation relevant for the process domain
          ! and corresponding indices in global observation vector
     
          ALLOCATE(thisobs%id_obs_f_lim(thisobs%dim_obs_g))
          CALL PDAFomi_get_local_ids_obs_f(thisobs%dim_obs_g, lradius, ocoord_g, dim_obs_f, thisobs%id_obs_f_lim)

          ! Store full and PE-local observation dimensions in module variables
          thisobs%dim_obs_p = dim_obs_p
          thisobs%dim_obs_f = dim_obs_f

          ! Allocate global observation arrays
          ! The arrays are deallocated in deallocate_obs in this module
          ALLOCATE(thisobs%obs_f(dim_obs_f))
          ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
          ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))

          ! Get process-relevant full observation arrays
          CALL PDAFomi_limit_obs_f(thisobs, 0, obs_g, thisobs%obs_f)
          CALL PDAFomi_limit_obs_f(thisobs, 0, ivar_obs_g, thisobs%ivar_obs_f)
          DO i = 1, ncoord
             CALL PDAFomi_limit_obs_f(thisobs, 0, ocoord_g(i,:), thisobs%ocoord_f(i,:))
          END DO

          DEALLOCATE(obs_g, ivar_obs_g, ocoord_g)

       END IF fullobs

    ELSE lfilter

       ! *** For global filters use process-local observations without gathering ***

       IF (mype_filter == 0) &
            WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use process-local observations for global filters'

       ! *** Initialize global dimension of observation vector ***
       dim_obs_f = dim_obs_p

       IF (mype_filter == 0) &
            WRITE (*, '(a, 8x, a, i7)') 'PDAFomi', &
            '--- Number of full observations ', dim_obs_f

       ! *** Gather full observation vector and corresponding coordinates ***

       ! Allocate full observation arrays
       ! The arrays are deallocated in deallocate_obs in this module
       ALLOCATE(thisobs%obs_f(dim_obs_f))
       ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
       ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))

       thisobs%obs_f = obs_p
       thisobs%ivar_obs_f = ivar_obs_p
       thisobs%ocoord_f = ocoord_p

       ! Store full and PE-local observation dimensions in module variables
       thisobs%dim_obs_p = dim_obs_p
       thisobs%dim_obs_f = dim_obs_f

    END IF lfilter


  END SUBROUTINE PDAFomi_gather_obs_f



!-------------------------------------------------------------------------------
!> Gather full observational information
!!
!! This routine uses PDAF_gather_obs_f_flex from PDAF to obtain
!! a full observed state vector. The routine is usually called
!! the observation operators.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_gather_obsstate_f(thisobs, obsstate_p, obsstate_f, offset)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in) :: obsstate_p(:)      !< Vector of process-local observed state
    REAL, INTENT(inout) :: obsstate_f(:)   !< Full observed vector for all types
    INTEGER, INTENT(inout) :: offset       !< input: offset of module-type observations in obsstate_f
                                           !< output: input + number of added observations

! *** Local variables ***
    INTEGER :: status                      ! Status flag for PDAF gather operation
    INTEGER :: localfilter                 ! Whether the filter is domain-localized
    REAL, ALLOCATABLE :: obsstate_tmp(:)   ! Temporary vector of globally full observations


! **************************************
! *** Gather full observation arrays ***
! **************************************

    ! Check  whether the filter is domain-localized
    CALL PDAF_get_localfilter(localfilter)

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug gather_obsstate_f: ', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
       WRITE (*,*) '++ OMI-debug gather_obsstate_f: ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
       WRITE (*,*) '++ OMI-debug gather_obsstate_f: ', debug, 'offset', offset
       IF (.NOT.thisobs%use_global_obs) &
            WRITE (*,*) '++ OMI-debug gather_obsstate_f: ', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
       WRITE (*,*) '++ OMI-debug gather_obsstate_f: ', debug, 'obsstate_p', obsstate_p
    END IF

    lfilter: IF (localfilter==1) THEN

       ! For domain-localized filters: gather full observations

       fullobs: IF (thisobs%use_global_obs) THEN

          ! *** Gather global full observation vector ***

          CALL PDAF_gather_obs_f_flex(thisobs%dim_obs_p, thisobs%dim_obs_f, obsstate_p, &
               obsstate_f(offset+1:offset+thisobs%dim_obs_f), status)

       ELSE fullobs

          ! *** Use full observations limited to those relevant for a process domain ***
          ! *** This can be more efficient as in the local analysis loop less        ***
          ! *** observations have a be checked for each analysis domain              ***

          ALLOCATE(obsstate_tmp(thisobs%dim_obs_g))

          ! *** Gather observation vector ***
          CALL PDAF_gather_obs_f_flex(thisobs%dim_obs_p, thisobs%dim_obs_g, obsstate_p, &
               obsstate_tmp, status)

          ! Now restrict observation vector to process-relevant part
          CALL PDAFomi_limit_obs_f(thisobs, offset, obsstate_tmp, obsstate_f)

          DEALLOCATE(obsstate_tmp)

       END IF fullobs

    ELSE lfilter

       ! *** For global filters use process-local observations without gathering ***

       ! In case of a global filter store process-local observed state
       obsstate_f(offset+1:offset+thisobs%dim_obs_p) = obsstate_p(1:thisobs%dim_obs_p)

    END IF lfilter

    ! Increment offset in observaton vector
    offset = offset + thisobs%dim_obs_f

  END SUBROUTINE PDAFomi_gather_obsstate_f



!-------------------------------------------------------------------------------
!> Initialize full vector of observations
!!
!! This routine initializes the part of the full vector of
!! observations for the current observation type.
!! It has to fill the observations to obsstate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!! The routine will only be called if the adaptive forgetting
!! factor is used.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obs_f(thisobs, dim_obs_f, obsstate_f, offset)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: dim_obs_f       !< Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(:)   !< Full observation vector (dim_obs_f)
    INTEGER, INTENT(inout) :: offset       !< input: offset of module-type observations in obsstate_f
                                           !< output: input + number of added observations


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    ! Consistency check
    IF (dim_obs_f < offset+thisobs%dim_obs_f) THEN
       WRITE (*,*) 'ERROR: PDAFomi_init_obs_f - dim_obs_f is too small'
    END IF

    doassim: IF (thisobs%doassim == 1) THEN

       ! Fill part of full observation vector
       obsstate_f(offset+1 : offset+thisobs%dim_obs_f) = thisobs%obs_f(1 : thisobs%dim_obs_f)

       ! Increment offset
       offset = offset + thisobs%dim_obs_f

    END IF doassim

  END SUBROUTINE PDAFomi_init_obs_f



!-------------------------------------------------------------------------------
!> Compute mean observation error variance
!!
!! This routine will only be called, if the adaptive
!! forgetting factor feature is used. Please note that
!! this is an experimental feature.
!!
!! The routine is called in global filters (like ESTKF)
!! during the analysis or in local filters (e.g. LESTKF)
!! before the loop over local analysis domains 
!! by the routine PDAF_set_forget that estimates an 
!! adaptive forgetting factor.  The routine has to 
!! initialize the mean observation error variance.  
!! For global filters this should be the global mean,
!! while for local filters it should be the mean for the
!! PE-local  sub-domain. (init_obsvar_l_TYPE is the 
!! localized variant for local filters)
!!
!! The routine assumes a diagonal observation error
!! covariance matrix.
!!
!! If the observation counter is zero the computation
!! of the mean variance is initialized. The output is 
!! always the mean variance. If the observation counter
!! is >0 the computed was initialized for another 
!! observation. To proceed the computation, meanvar is
!! first multiplied with the observation counter to 
!! obtain the variance sum. Then the computation of the 
!! mean is continued.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obsvar_f(thisobs, meanvar, cnt_obs)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(inout) :: meanvar         !< Mean variance
    INTEGER, INTENT(inout) :: cnt_obs      !< Observation counter

! Local variables ***
    INTEGER :: i        ! Counter


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (cnt_obs==0) THEN
          ! Reset mean variance
          meanvar = 0.0
       ELSE
          ! Compute sum of variances from mean variance
          meanvar = meanvar * REAL(cnt_obs)
       END IF

       ! Add observation error variances
       DO i = 1, thisobs%dim_obs_f
          meanvar = meanvar + 1.0 / thisobs%ivar_obs_f(i)
       END DO

       ! Increment observation count
       cnt_obs = cnt_obs + thisobs%dim_obs_f

       ! Compute updated mean variance
       meanvar = meanvar / REAL(cnt_obs)

    END IF doassim

  END SUBROUTINE PDAFomi_init_obsvar_f



!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix
!!
!! The routine is called during the analysis step
!! of the global square-root filters. It has to 
!! compute the product of the inverse of the
!! process-local observation error covariance matrix
!! with the matrix of process-local observed ensemble 
!! perturbations.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_prodRinvA(thisobs, ncols, A_p, C_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs !< Data type with full observation
    INTEGER, INTENT(in) :: ncols          !< Number of columns in A_p and C_p
    REAL, INTENT(in) :: A_p(:, :)         !< Input matrix (nobs_f, ncols)
    REAL, INTENT(out)   :: C_p(:, :)      !< Output matrix (nobs_f, ncols)


! *** local variables ***
    INTEGER :: i, j       ! index of observation component
    INTEGER :: off        ! row offset in A_l and C_l
    

! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Check process-local observation dimension

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when localfilter=1
          WRITE (*,*) 'ERROR: INCONSISTENT value for DIM_OBS_P'
       END IF

       ! Initialize offset
       off = thisobs%off_obs_f

       DO j = 1, ncols
          DO i = 1, thisobs%dim_obs_f
             C_p(i+off, j) = thisobs%ivar_obs_f(i) * A_p(i+off, j)
          END DO
       END DO

    END IF doassim

  END SUBROUTINE PDAFomi_prodRinvA



!-------------------------------------------------------------------------------
!> Compute likelihood for an ensemble member
!!
!! The routine is called during the analysis step
!! of the NETF or a particle filter.
!! It has to compute the likelihood of the
!! ensemble according to the difference from the
!! observation (residual) and the error distribution
!! of the observations.
!!
!! In general this routine is similar to the routine
!! prodRinvA used for ensemble square root Kalman
!! filters. As an addition to this routine, we here have
!! to evaluate the likelihood weight according the
!! assumed observation error statistics.
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_likelihood(thisobs, nobs, obs, resid, lhood)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs   !< Data type with full observation
    INTEGER, INTENT(in) :: nobs          !< Number of observations
    REAL, INTENT(in)    :: obs(:)        ! PE-local vector of observations
    REAL, INTENT(in)    :: resid(:)      ! Input vector of residuum
    REAL, INTENT(out)   :: lhood         ! Output vector - log likelihood

! *** local variables ***
    INTEGER :: i         ! index of observation component
    REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual
    REAL :: lhood_one    ! Likelihood for this observation
    REAL :: rdummy       ! Dummy to access observation_l


    doassim: IF (thisobs%doassim == 1) THEN

! ****************************************
! *** First scale by observation error ***
! ***                   -1             ***
! ***      Rinvresid =  R  resid       ***
! ***                                  ***
! *** We assume a diagonal matrix R    ***
! ****************************************

       ! Initialize dummy to prevent compiler warning
       rdummy = obs(1)

       ALLOCATE(Rinvresid(nobs))

       DO i = 1, nobs
          Rinvresid(i) = thisobs%ivar_obs_f(i) * resid(i)
       END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

       IF (thisobs%obs_err_type==0) THEN

          ! Gaussian errors
          ! Calculate exp(-0.5*resid^T*R^-1*resid)

          ! Transform pack to log likelihood to increment its values
          IF (lhood>0.0) lhood = - LOG(lhood)

          CALL dgemv('t', nobs, 1, 0.5, resid, &
               nobs, Rinvresid, 1, 0.0, lhood_one, 1)

          lhood = EXP(-(lhood + lhood_one))

       ELSE

          ! Double-exponential errors
          ! Calculate exp(-SUM(ABS(resid)))

          ! Transform pack to log likelihood to increment its values
          IF (lhood>0.0) lhood = - LOG(lhood)

          lhood_one = 0.0
          DO i = 1, nobs
             lhood_one = lhood_one + ABS(Rinvresid(i))
          END DO

          lhood = EXP(-(lhood + lhood_one))

       END IF

       ! *** Clean up ***

       DEALLOCATE(Rinvresid)

    END IF doassim
    
  END SUBROUTINE PDAFomi_likelihood



!-------------------------------------------------------------------------------
!> Add observation error to some matrix
!!
!! The routine is called during the analysis step
!! of the stochastic EnKF. It it provided with a
!! matrix in observation space and has to add the 
!! observation error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_add_obs_error(thisobs, nobs_all, matC)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs      !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all         !< Number of observations
    REAL, INTENT(inout) :: matC(:, :)       !< Input/Output matrix (nobs_f, rank)


! *** local variables ***
    INTEGER :: i, i_all         ! index of observation component
    INTEGER :: idummy           ! Dummy to access nobs_all


    doassim: IF (thisobs%doassim == 1) THEN

      ! Initialize dummy to prevent compiler warning
       idummy = nobs_all


! *************************************************
! *** Check process-local observation dimension ***
! *************************************************

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when localfilter=1
          WRITE (*,*) 'PDAFomi ERROR: INCONSISTENT  VALUE for DIM_OBS_P'
       END IF


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

       DO i = 1, thisobs%dim_obs_f
          i_all = i + thisobs%off_obs_f
          matC(i_all, i_all) = matC(i_all, i_all) + 1.0/thisobs%ivar_obs_f(i)
       ENDDO

    END IF doassim

  END SUBROUTINE PDAFomi_add_obs_error



!-------------------------------------------------------------------------------
!> Initialize global observation error covariance matrix
!!
!! The routine is called during the analysis
!! step when an ensemble of observations is
!! generated by PDAF_enkf_obs_ensemble. 
!! It has to initialize the global observation 
!! error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obscovar(thisobs, nobs_all, covar, isdiag)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all        !< Number of observations
    REAL, INTENT(out) :: covar(:, :)       !< Input/Output matrix (nobs_f, rank)
    LOGICAL, INTENT(out) :: isdiag         !< Whether matrix R is diagonal

! *** local variables ***
    INTEGER :: i, i_all         ! index of observation component
    INTEGER :: idummy           ! Dummy to access nobs_all


    doassim: IF (thisobs%doassim == 1) THEN

       ! Initialize dummy to prevent compiler warning
       idummy = nobs_all


! *************************************
! ***   Initialize covariances      ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

       covar(:, :) = 0.0

       DO i = 1, thisobs%dim_obs_f
          i_all = i + thisobs%off_obs_f
          covar(i_all, i_all) = covar(i_all, i_all) + 1.0/thisobs%ivar_obs_f(i)
       ENDDO

       ! The matrix is diagonal
       ! This setting avoids the computation of the SVD of COVAR
       ! in PDAF_enkf_obs_ensemble
       isdiag = .TRUE.

    END IF doassim
    
  END SUBROUTINE PDAFomi_init_obscovar



!-------------------------------------------------------------------------------
!> Initialize vector of observation errors for generating synthetic obs.
!!
!! This routine initializes a vector of observation errors
!! used to perturb an observe dmodel state when generating 
!! synthetic observations.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obserr_f(thisobs, obserr_f)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs  !< Data type with full observation
    REAL, INTENT(inout) :: obserr_f(:)     !< Full vector of observation errors

! *** Local variables ***
    INTEGER :: i                           ! Counter


! *****************************************************************************
! *** Initialize vector of observation errors for generating synthetic obs. ***
! *****************************************************************************

    doassim: IF (thisobs%doassim == 1) THEN

       DO i = 1, thisobs%dim_obs_f
          IF (thisobs%ivar_obs_f(i)>0.0) THEN
             obserr_f(i + thisobs%off_obs_f) = 1.0/SQRT(thisobs%ivar_obs_f(i))
          ELSE    
             obserr_f(i + thisobs%off_obs_f) = 1.0e12
          END IF
       END DO

    END IF doassim

  END SUBROUTINE PDAFomi_init_obserr_f



!-------------------------------------------------------------------------------
!> Deallocate arrays in observation type
!!
!! This routine deallocates arrays in the data type THISOBS.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2019-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_deallocate_obs(thisobs)

    IMPLICIT NONE

! *** Arguments
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation

   ! *** Perform deallocation ***

    IF (ALLOCATED(thisobs%obs_f)) DEALLOCATE(thisobs%obs_f)
    IF (ALLOCATED(thisobs%ocoord_f)) DEALLOCATE(thisobs%ocoord_f)
    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
    IF (ALLOCATED(thisobs%ivar_obs_f)) DEALLOCATE(thisobs%ivar_obs_f)
    IF (ALLOCATED(thisobs%icoeff_p)) DEALLOCATE(thisobs%icoeff_p)
    IF (ALLOCATED(thisobs%domainsize)) DEALLOCATE(thisobs%domainsize)
    IF (ALLOCATED(thisobs%id_obs_f_lim)) DEALLOCATE(thisobs%id_obs_f_lim)

  END SUBROUTINE PDAFomi_deallocate_obs



!-------------------------------------------------------------------------------
!> Set min/max coordinate locations of a decomposed grid
!!
!! This routine sets the limiting coordinates of a 
!! process domain, i.e. the northern-, southern-,
!! eastern-, and western-most coordinate. The
!! information can be used to restrict the full
!! observations for PDAF to those that might be
!! used for the local analysis. 
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_domain_limits(lim_coords)

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in) :: lim_coords(2,2)     !< geographic coordinate array (1: longitude, 2: latitude)
                                            !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)

    ! Store domain limiting coordinates in module array
    IF (.NOT.ALLOCATED(domain_limits)) ALLOCATE(domain_limits(4))

    domain_limits(1) = lim_coords(2,1)
    domain_limits(2) = lim_coords(2,2)
    domain_limits(3) = lim_coords(1,1)
    domain_limits(4) = lim_coords(1,2)

  END SUBROUTINE PDAFomi_set_domain_limits



!-------------------------------------------------------------------------------
!> Find min/max coordinate locations in unstructured grid
!!
!! This routine finds the limiting coordinates of a 
!! process domain, i.e. the northern-, southern-,
!! eastern-, and western-most coordinate. The
!! information can be used to restrict the full
!! observations for PDAF to those that might be
!! used for the local analysis. Usually this is used
!! for unstructured grids, like in FESOM, because the
!! grid point indices do not contain information on 
!! coordinates in this case.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_domain_limits_unstr(verbose, npoints_p, coords_p)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose          !< verbosity flag 
    INTEGER, INTENT(in) :: npoints_p        !< number of process-local grid points
    REAL, INTENT(in) :: coords_p(:,:)       !< geographic coordinate array (row 1: longitude, 2: latitude)
                                            !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)

! *** Local variables ***
    INTEGER :: i                            ! Counter
    REAL :: nlimit, slimit, elimit, wlimit  ! Limiting coordinates
    REAL :: abslonmin                       ! absolute minimum longitude


! *** Determine limiting coordinates ***

    ! Initialize limiting values
    nlimit = -100.0
    slimit = 100.0
    wlimit = 100.0
    elimit = -100.0
    abslonmin = 100.0

    DO i=1, npoints_p
       ! Get North/South Limits
       IF (coords_p(2,i) < slimit) slimit = coords_p(2,i)
       IF (coords_p(2,i) > nlimit) nlimit = coords_p(2,i)

       ! Get East/West Limits
       IF (coords_p(1,i) < wlimit) wlimit = coords_p(1,i)
       IF (coords_p(1,i) > elimit) elimit = coords_p(1,i)
       IF (ABS(coords_p(1,i)) < abslonmin) THEN
          abslonmin = ABS(coords_p(1,i))
       END IF
    ENDDO

    IF (elimit*wlimit<0.0) THEN
       ! Domain crosses prime meridian or date line

       IF (wlimit<-3.1 .AND. elimit>3.1 .and. abslonmin>0.5) THEN

          ! If the domain crosses the date line, we have to search the longitudinal limits differently
          elimit = -100.0
          wlimit = 100.0
          DO i=1, npoints_p
             IF (coords_p(1,i)<0.0 .AND. coords_p(1,i)>elimit) elimit = coords_p(1,i)
             IF (coords_p(1,i)>0.0 .AND. coords_p(1,i)<wlimit) wlimit = coords_p(1,i)
          END DO
          IF (verbose==1) &
               WRITE (*,'(i3,x,a,4f10.3,a)') mype_filter, 'limit coords', nlimit, slimit, wlimit, elimit, '+++'
       ELSE
          ! In this case the domain crosses the prime meridian
          IF (verbose==1) &
               WRITE (*,'(i3,x,a,4f10.3,a)') mype_filter, 'limit coords', nlimit, slimit, wlimit, elimit, '---'
       END IF
    ELSE
       ! Standard case
       IF (verbose==1) &
            WRITE (*,'(i3,x,a,4f10.3)') mype_filter, 'limit coords', nlimit, slimit, wlimit, elimit
    END IF

    ! Store domain limiting coordinates in module array

    IF (.NOT.ALLOCATED(domain_limits)) ALLOCATE(domain_limits(4))

    domain_limits(1) = nlimit
    domain_limits(2) = slimit
    domain_limits(3) = wlimit
    domain_limits(4) = elimit

  END SUBROUTINE PDAFomi_get_domain_limits_unstr


  
!-------------------------------------------------------------------------------
!> Find observations inside or close to process domain
!!
!! This routine finds observations that lie inside the 
!! local process sub-domain or within the distance
!! LRADIUS around it. The observations are counted and
!! an index array is initialized storing the indices
!! of the process-local relevant full observations in the
!! global full observation vector.
!!
!! The routine has to be called by all filter processes
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_local_ids_obs_f(dim_obs_g, lradius, oc_f, cnt_lim, id_lim)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_g       !< Global full number of observations
    REAL, INTENT(in) :: lradius            !< Localization radius (used is a constant one here)
    REAL, INTENT(in) :: oc_f(:,:)          !< observation coordinates (radians), row 1: lon, 2: lat
                                           !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)
    INTEGER, INTENT(out) :: cnt_lim        !< Number of full observation for local process domain
    INTEGER, INTENT(out) :: id_lim(:)      !< Indices of process-local full obs. in global full vector
                                           !< It has to be allocated sufficiently large

! *** Local variables ***
    INTEGER :: i         ! Counter
    INTEGER :: flag      ! Counting flag
    REAL :: limdist      ! Limit distance normalized by r_earth
    INTEGER :: cnt_lim_max, cnt_lim_min  ! min/max number over all domains
    REAL :: maxlat       ! Highest latitude of a domain


! **********************
! *** Initialization ***
! **********************

    IF (.NOT.ALLOCATED(domain_limits)) THEN
       WRITE (*,*) 'ERROR: PDAFomi_get_local_ids_obs_f - DOMAIN_LIMITS is not initialized'
    END IF

    ! initialize index array
    id_lim = 0

    ! Limit distance around the domain
    limdist = lradius / r_earth


! ***************************************
! *** Find relevant full observations ***
! ***************************************

    cnt_lim = 0

    fullobsloop: DO i = 1, dim_obs_g

       ! Init flag for latitudinal check
       flag = 0

       ! First check in latitudinal direction
       checklat: IF (oc_f(2,i)<=domain_limits(1) .AND. oc_f(2,i)>=domain_limits(2)) THEN
          ! inside domain north-south extent
          flag=1
       ELSEIF (oc_f(2,i)>domain_limits(1)) THEN
          ! north of the domain
          IF (ABS(oc_f(2,i)-domain_limits(1)) <= limdist) flag=1
       ELSEIF (oc_f(2,i)<domain_limits(2)) THEN
          ! south of the domain
          IF (ABS(oc_f(2,i)-domain_limits(2)) <= limdist) flag=1
       END IF checklat

       ! Store highest latitude
       maxlat = MAX(ABS(domain_limits(1)), ABS(domain_limits(2)))

       ! if observation fits in the latitudinal direction check longitudinal direction
       lat_ok: IF (flag==1) THEN
          lontypes: IF (domain_limits(4)>=0.0 .OR. (domain_limits(4)<0.0 .AND. domain_limits(3)<0.0)) THEN

             IF (oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=domain_limits(4)) THEN

                ! fully inside domain extent
                cnt_lim = cnt_lim+1
                id_lim(cnt_lim) = i
             ELSEIF (oc_f(1,i)<domain_limits(3)) THEN

                ! west of the domain
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3))) <= limdist .or. &
                    (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3)))-2.0*pi) <= limdist ) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                ! east of the domain
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4))) <= limdist .or. &
                    (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4)))-2.0*pi) <= limdist ) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ENDIF
          ELSE lontypes
             IF ((oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=pi) .OR. &
                  (oc_f(1,i)<=domain_limits(4)) .AND. oc_f(1,i)>=-pi) THEN

                ! fully inside domain extent
                cnt_lim = cnt_lim+1
                id_lim(cnt_lim) = i

             ELSEIF (oc_f(1,i)<domain_limits(3) .AND. oc_f(1,i)>=0.0) THEN

                ! east of the domain
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3))) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                ENDIF
             ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                ! west of the domain
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4))) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                ENDIF
             ENDIF
          ENDIF lontypes
       ENDIF lat_ok
    END DO fullobsloop

    ! Get number of min/max process-local full observation dimensions
    CALL MPI_Allreduce (cnt_lim, cnt_lim_max, 1, MPI_INTEGER, MPI_MAX, &
         COMM_filter, MPIerr)
    CALL MPI_Allreduce (cnt_lim, cnt_lim_min, 1, MPI_INTEGER, MPI_MIN, &
         COMM_filter, MPIerr)
  
    IF (mype_filter==0) THEN
       WRITE (*,'(a,8x,a,i8)') 'PDAFomi','--- global obs. dimension', dim_obs_g
       WRITE (*,'(a,8x,a,i7,1x,i7)') 'PDAFomi','--- process-local min/max full obs. dimensions', &
            cnt_lim_min, cnt_lim_max
    END IF

  END SUBROUTINE PDAFomi_get_local_ids_obs_f


  
!-------------------------------------------------------------------------------
!> Reduce full observation vector to part relevant for local process domain
!!
!! This routine initializes a full vector of observations that only
!! contains those full observations that are relevant for a process
!! subdomain. The indices of these observations were determined
!! using get_local_ids_obs_f.
!!
!! __Revision history:__
!! * 2019-07 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_limit_obs_f(thisobs, offset, obs_f_one, obs_f_lim)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in) :: obs_f_one(:)       !< Global full observation vector (nobs_f)
    REAL, INTENT(out) :: obs_f_lim(:)      !< full observation vector for process domains (nobs_lim)
    INTEGER, INTENT(in) :: offset          !< offset of this observation in obs_f_lim

! *** Local variables ***
    INTEGER :: i         ! Counter


! ********************************************
! *** Initialize process-local full vector ***
! ********************************************

    IF (.NOT.ALLOCATED(thisobs%id_obs_f_lim)) THEN
       WRITE (*,*) 'ERROR: PDAFomi_limit_obs_f - thisobs%id_obs_f_lim is not allocated'
    END IF

    DO i = 1, thisobs%dim_obs_f
       obs_f_lim(i+offset) = obs_f_one(thisobs%id_obs_f_lim(i))
    END DO

  END SUBROUTINE PDAFomi_limit_obs_f

END MODULE PDAFomi_obs_f
