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
!! * init_obs_f \n
!!        Initialize full vector of observations for adaptive forgetting factor
!! * init_obsvar_f \n
!!        Compute mean observation error variance for adaptive forgetting factor
!! * deallocate_obs \n
!!        Deallocate arrays in observation type
!! * get_domain_limits_unstr \n
!!        Find min/max coordinate locations in unstructured grid
!! * get_local_ids_obs_f \n
!!        Find observations inside or close to process domain
!! * limit_obs_f \n
!!        Reduce full observation vector to part relevant for local process domain
!! * prodRinvA \n
!!        Multiply an intermediate matrix of the global filter analysis
!!        with the inverse of the observation error covariance matrix
!!
!! The coordinates are assumed to be in radians and are within the range 
!! -pi to +pi for longitude (- is westward) and -pi/2 to +pi/2 for latitude.
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
  REAL :: domain_limits(4)             !< Limiting coordinates (NSWE) for process domain
  REAL, PARAMETER :: r_earth=6.3675e6  !< Earth radius in meters
  REAL, PARAMETER :: pi=3.141592653589793   !< Pi

! *** Data type to define the full observations by internally shared variables of the module
  type obs_f
     ! ---- Mandatory variables to be set in init_dim_obs_f ----
     INTEGER :: doassim=0                 !< Whether to assimilate this observation type
     INTEGER :: dim_obs_p                 !< number of PE-local observations
     INTEGER :: dim_obs_f                 !< number of full observations
     INTEGER, ALLOCATABLE :: id_obs_p(:,:) !< indices of observed field in state vector
     REAL, ALLOCATABLE :: obs_f(:)        !< Full observed field
     REAL, ALLOCATABLE :: ocoord_f(:,:)   !< Coordinates of full observation vector
     REAL, ALLOCATABLE :: ivar_obs_f(:)   !< Inverse variance of full observations
     INTEGER :: disttype                  !< Type of distance computation to use for localization
     INTEGER :: ncoord                    !< Number of coordinates use for distance computation
     ! ---- Variables with predefined values - they can be changed in init_dim_obs_f  ----
     INTEGER :: obs_err_type=0            !< Type of observation error: (0) Gauss, (1) Laplace
     LOGICAL :: localfilter=.true.        !< Whether a localized filter is used
     ! ---- Optional variables - they can be set in init_dim_obs_f ----
     REAL, ALLOCATABLE :: icoeff_p(:,:)   !< Interpolation coefficients for obs. operator (optional)
     REAL, ALLOCATABLE :: domainsize(:)   !< Size of domain for periodicity (<=0 for no periodicity) (optional)
     ! ---- Mandatory variable to be set in obs_op_f ---
     INTEGER :: off_obs_f                 !< Offset of this observation in overall full obs. vector
  end type obs_f


!-------------------------------------------------------------------------------
  
CONTAINS
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
  SUBROUTINE init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs        !< Data type with full observation
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) !< Full observation vector
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in obsstate_f
                                                 !< output: input + number of added observations


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    ! Fill part of full observation vector
    obsstate_f(offset_obs+1 : offset_obs+thisobs%dim_obs_f) = thisobs%obs_f(1 : thisobs%dim_obs_f)

    ! Increment offset
    offset_obs = offset_obs + thisobs%dim_obs_f

  END SUBROUTINE init_obs_f



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
  SUBROUTINE init_obsvar_f(thisobs, meanvar, cnt_obs)

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

  END SUBROUTINE init_obsvar_f



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
  SUBROUTINE deallocate_obs(thisobs)

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

  END SUBROUTINE deallocate_obs


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
  SUBROUTINE get_domain_limits_unstr(verbose, npoints_p, coords_p)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose              !< verbosity flag 
    INTEGER, INTENT(in) :: npoints_p            !< number of process-local grid points
    REAL, INTENT(in) :: coords_p(2,npoints_p)   !< geographic coordinate array (1: longitude, 2: latitude)
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
    domain_limits(1) = nlimit
    domain_limits(2) = slimit
    domain_limits(3) = wlimit
    domain_limits(4) = elimit

  END SUBROUTINE get_domain_limits_unstr


  
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
  SUBROUTINE get_local_ids_obs_f(dim_obs_f, lradius, oc_f, cnt_lim, id_lim)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_f       !< Global full number of observations
    REAL, INTENT(in) :: lradius            !< Localization radius (used is a constant one here)
    REAL, INTENT(in) :: oc_f(2,dim_obs_f)  !< observation coordinates (radians)
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

    ! initialize index array
    id_lim = 0

    ! Limit distance around the domain
    limdist = lradius / r_earth


! ***************************************
! *** Find relevant full observations ***
! ***************************************

    cnt_lim = 0

    fullobsloop: DO i = 1, dim_obs_f

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
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3))) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                ! east of the domain
                IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4))) <= limdist) THEN
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
       WRITE (*,'(a,3x,a,i8)') 'PDAFomi','--- overall full obs. dimension', dim_obs_f
       WRITE (*,'(a,3x,a,2i6)') 'PDAFomi','--- process-local min/max full obs. dimensions', &
            cnt_lim_min, cnt_lim_max
    END IF

  END SUBROUTINE get_local_ids_obs_f


  
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
  SUBROUTINE limit_obs_f(nobs_f, nobs_f_lim, id_lim, obs_f, obs_f_lim)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs_f              !< Global full number of observations
    INTEGER, INTENT(in) :: nobs_f_lim          !< Number of full observations for process domain
    REAL, INTENT(in) :: obs_f(nobs_f)          !< Global full observation vector
    INTEGER, INTENT(in) :: id_lim(nobs_f_lim)  !< Indices of process-local full obs. in global full vector
    REAL, INTENT(out) :: obs_f_lim(nobs_f_lim) !< full observation vector for process domains

! *** Local variables ***
    INTEGER :: i         ! Counter


! ********************************************
! *** Initialize process local full vector ***
! ********************************************

    DO i = 1, nobs_f_lim
       obs_f_lim(i) = obs_f(id_lim(i))
    END DO

  END SUBROUTINE limit_obs_f



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
  SUBROUTINE prodRinvA(nobs_p, rank, ivar_obs_p, A_p, C_p)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs_p         !< Dimension of obs. vector (one or all obs. types)
    INTEGER, INTENT(in) :: rank           !< Rank of initial covariance matrix
    REAL, INTENT(in)    :: ivar_obs_p(:)  !< PE-local vector of inverse obs. variances (nobs_f)
    REAL, INTENT(in) :: A_p(:, :)         !< Input matrix (nobs_f, rank)
    REAL, INTENT(out)   :: C_p(:, :)      !< Output matrix (nobs_f, rank)


! *** local variables ***
    INTEGER :: i, j       ! index of observation component
    

! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

    DO j = 1, rank
       DO i = 1, nobs_p
          C_p(i, j) = ivar_obs_p(i) * A_p(i, j)
       END DO
    END DO

  END SUBROUTINE prodRinvA



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
  SUBROUTINE add_obs_error(nobs, ivar_obs_one, matC, offset)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs         !< Number of observations
    REAL, INTENT(in)    :: ivar_obs_one(:)  !< vector of inverse obs. variances (nobs_f)
    REAL, INTENT(inout) :: matC(:, :)   !< Input/Output matrix (nobs_f, rank)
    INTEGER, INTENT(in) :: offset       !< Offset of this observation in overall obs. vector


! *** local variables ***
    INTEGER :: i, i_all       ! index of observation component


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

  DO i = 1, nobs
     i_all = i + offset
     matC(i_all, i_all) = matC(i_all, i_all) + 1.0/ivar_obs_one(i)
  ENDDO

END SUBROUTINE add_obs_error



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
  SUBROUTINE init_obscovar(nobs, ivar_obs, offset, covar, isdiag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs         !< Number of observations
    REAL, INTENT(in)    :: ivar_obs(:)  !< vector of inverse obs. variances (nobs_f)
    REAL, INTENT(out) :: covar(:, :)    !< Input/Output matrix (nobs_f, rank)
    LOGICAL, INTENT(out) :: isdiag      !< Whether matrix R is diagonal
    INTEGER, INTENT(in) :: offset       !< Offset of this observation in overall obs. vector


! *** local variables ***
    INTEGER :: i, i_all         ! index of observation component


! *************************************
! ***   Initialize covariances      ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

    covar(:, :) = 0.0

    DO i = 1, nobs
       i_all = i + offset
       covar(i_all, i_all) = covar(i_all, i_all) + 1.0/ivar_obs(i)
    ENDDO

    ! The matrix is diagonal
    ! This setting avoids the computation of the SVD of COVAR
    ! in PDAF_enkf_obs_ensemble
    isdiag = .TRUE.
    
  END SUBROUTINE init_obscovar



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
  SUBROUTINE likelihood(nobs, obs, resid, ivar_obs, lhood, obs_err_type)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs          !< Number of observations
    REAL, INTENT(in)    :: obs(:)        ! PE-local vector of observations
    REAL, INTENT(in)    :: resid(:)      ! Input vector of residuum
    REAL, INTENT(in)    :: ivar_obs(:)   !< vector of inverse obs. variances (nobs_f)
    REAL, INTENT(out)   :: lhood         ! Output vector - log likelihood
    INTEGER, INTENT(in) :: obs_err_type  ! Type of observation error
                                         ! (0) Gaussian, (1) Laplace/double exponential

! *** local variables ***
    INTEGER :: i         ! index of observation component
    REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual
    REAL :: lhood_one    ! Likelihood for this observation


! ****************************************
! *** First scale by observation error ***
! ***                   -1             ***
! ***      Rinvresid =  R  resid       ***
! ***                                  ***
! *** We assume a diagonal matrix R    ***
! ****************************************

    ALLOCATE(Rinvresid(nobs))

    DO i = 1, nobs
       Rinvresid(i) = ivar_obs(i) * resid(i)
    END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

    IF (obs_err_type==0) THEN

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
    
  END SUBROUTINE likelihood

END MODULE PDAFomi_obs_f
