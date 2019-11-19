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
!$Id: PDAFomi_obs_f.F90 240 2019-10-23 13:30:50Z lnerger $
!BOP
!
! !MODULE:
MODULE PDAFomi_obs_f
!
! !DESCRIPTION:
! This module contains subroutines to restrict the global full vector of observations
! to those observations that are relevant for a process-local model subdomain.
!
! The coordinates are assumed to be in radians and are within the range 
! -pi to +pi for longitude (- is westward) and -pi/2 to +pi/2 for latitude.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, COMM_FILTER, MPI_INTEGER, MPI_SUM, MPIerr, MPI_MIN, MPI_MAX

  IMPLICIT NONE
  SAVE

! *** Module internal variables
  REAL :: domain_limits(4)             ! Limiting coordinates (NSWE) for process domain
  REAL, PARAMETER :: r_earth=6.3675e6  ! Earth radius in meters
  REAL, PARAMETER ::  pi=3.141592653589793   ! Pi

  ! Data type to define the full observations by internally shared variables of the module
  type obs_f
     INTEGER :: dim_obs_p                 ! number of PE-local observations
     INTEGER :: dim_obs_f                 ! number of full observations
     INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
     INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
     REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
     REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
     REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
     INTEGER :: disttype                  ! Type of distance computation to use for localization
     INTEGER :: ncoord                    ! Number of coordinates use for distance computation
  end type obs_f

! EOP  
!-------------------------------------------------------------------------------
  
CONTAINS
!BOP
!
! !ROUTINE: init_obs_f --- Initialize full vector of observations
!
! !INTERFACE:
  SUBROUTINE init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)

! !DESCRIPTION:
! This routine initializes the part of the full vector of
! observations for the current observation type.
! It has to fill the observations to obsstate_f from
! position OFFSET_OBS+1. For the return value OFFSET_OBS
! has to be incremented by the number of added observations.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    TYPE(obs_f), INTENT(inout) :: thisobs        ! Information on full observation
    INTEGER, INTENT(in) :: dim_obs_f             ! Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) ! Full observation vector
    INTEGER, INTENT(inout) :: offset_obs         ! input: offset of module-type observations in obsstate_f
                                                 ! output: input + number of added observations
!EOP


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    ! Fill part of full observation vector
    obsstate_f(offset_obs+1 : offset_obs+thisobs%dim_obs_f) = thisobs%obs_f(1 : thisobs%dim_obs_f)

    ! Increment offset
    offset_obs = offset_obs + thisobs%dim_obs_f

  END SUBROUTINE init_obs_f



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obsvar_f --- Compute mean observation error variance
!
! !INTERFACE:
  SUBROUTINE init_obsvar_f(thisobs, meanvar, cnt_obs)

! !DESCRIPTION:
! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like ESTKF)
! during the analysis or in local filters (e.g. LESTKF)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain. (init_obsvar_l_TYPE is the 
! localized variant for local filters)
!
! The implemented functionality is generic. There 
! should be no changes required as long as the 
! observation error covariance matrix is diagonal.
!
! If the observation counter is zero the computation
! of the mean variance is initialized. The output is 
! always the mean variance. If the observation counter
! is >0 first the variance sum is computed by 
! multiplying with the observation counter.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    TYPE(obs_f), INTENT(inout) :: thisobs  ! Information on full observation
    REAL, INTENT(inout) :: meanvar         ! Mean variance
    INTEGER, INTENT(inout) :: cnt_obs      ! Observation counter
!EOP

! Local variables
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
!BOP
!
! !ROUTINE: deallocate_obs --- Deallocate arrays in observation type
!
! !INTERFACE:
  SUBROUTINE deallocate_obs(thisobs)

! !DESCRIPTION:
! This routine deallocates arrays in the data type THISOBS.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    TYPE(obs_f), INTENT(inout) :: thisobs  ! Information on full observation

   ! *** Perform deallocation ***

    IF (ALLOCATED(thisobs%obs_f)) DEALLOCATE(thisobs%obs_f)
    IF (ALLOCATED(thisobs%ocoord_f)) DEALLOCATE(thisobs%ocoord_f)
    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
    IF (ALLOCATED(thisobs%ivar_obs_f)) DEALLOCATE(thisobs%ivar_obs_f)

  END SUBROUTINE deallocate_obs


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_domain_limits_unstr - find min/max coordinate locations in unstructured grid
!
! !INTERFACE:
  SUBROUTINE get_domain_limits_unstr(verbose, npoints_p, coords_p)

! !DESCRIPTION:
! This routine find the limiting coordinates of a 
! process domain, i.e. the northern-, southern-,
! eastern-, and western-most coordinate. The
! information can be used to restrict the full
! observations for PDAF to those that might be
! used for the local analysis.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: verbose              ! verbosity flag 
    INTEGER, INTENT(in) :: npoints_p            ! number of process-local grid points
    REAL, INTENT(in) :: coords_p(2,npoints_p)   ! geographic coordinate array (1: longitude, 2: latitude)
                                                ! ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)
!EOP

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
!BOP
!
! !ROUTINE: get_local_ids_obs_f - find observations inside or close to process domain
!
! !INTERFACE:
  SUBROUTINE get_local_ids_obs_f(dim_obs_f, lradius, oc_f, cnt_lim, id_lim)


! !DESCRIPTION:
! This routine finds observations that lie inside the 
! local process sub-domain or within the distance
! LRADIUS around it. The observations are counted and
! an index array is initialized storing the indices
! of the process-local relevant full observations in the
! global full observation vector.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS
    INTEGER, INTENT(in) :: dim_obs_f       ! Global full number of observations
    REAL, INTENT(in) :: lradius            ! Localization radius (used is a constant one here)
    REAL, INTENT(in) :: oc_f(2,dim_obs_f)  ! observation coordinates (radians)
                                           ! ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)
    INTEGER, INTENT(out) :: cnt_lim        ! Number of full observation for local process domain
    INTEGER, INTENT(out) :: id_lim(:)      ! Indices of process-local full obs. in global full vector
!EOP

! Local variables
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
       WRITE (*,'(a,3x,a,i8)') 'PDAF-USER','--- overall full obs. dimension', dim_obs_f
       WRITE (*,'(a,3x,a,2i6)') 'PDAF-USER','--- process-local min/max full obs. dimensions', &
            cnt_lim_min, cnt_lim_max
    END IF

  END SUBROUTINE get_local_ids_obs_f


  
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: limit_obs_f - reduce full observation vector to part relevant for local process domain
!
! !INTERFACE:
  SUBROUTINE limit_obs_f(nobs_f, nobs_f_lim, id_lim, obs_f, obs_f_lim)


! !DESCRIPTION:
! This routine initializes a full vector of observations that only
! contains those full observations that are relevant for a process
! subdomain. The indices of these observations were determined
! using GET_LOCAL_IDS_OBS_F.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS
    INTEGER, INTENT(in) :: nobs_f              ! Global full number of observations
    INTEGER, INTENT(in) :: nobs_f_lim          ! Number of full observations for process domain
    REAL, INTENT(in) :: obs_f(nobs_f)          ! Global full observation vector
    INTEGER, INTENT(in) :: id_lim(nobs_f_lim)  ! Indices of process-local full obs. in global full vector
    REAL, INTENT(out) :: obs_f_lim(nobs_f_lim) ! full observation vector for process domains
!EOP

! Local variables
    INTEGER :: i         ! Counter


! *** Initialize process local full vector ***

    DO i = 1, nobs_f_lim
       obs_f_lim(i) = obs_f(id_lim(i))
    END DO

  END SUBROUTINE limit_obs_f

END MODULE PDAFomi_obs_f
