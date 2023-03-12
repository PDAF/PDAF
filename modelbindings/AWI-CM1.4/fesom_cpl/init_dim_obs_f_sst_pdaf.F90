!$Id: init_dim_obs_f_pdaf.F90 1981 2019-03-12 14:51:18Z qtang $
!BOP
!
! !ROUTINE: init_dim_obs_f_sst_pdaf --- Set full dimension of SST observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_sst_pdaf(step, dim_obs_f)

! !DESCRIPTION:
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2019-03 - Lars Nerger - Initial code from splitting init_dim_obs_f_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, COMM_filter
  USE mod_assim_pdaf, &
       ONLY: delt_obs_ocn, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       obs_sst, obs_sst_error, ivariance_obs, id_observed_sst_p, ocoord_n2d, &
       rms_obs, dim_obs_p, filtertype, &
       mean_ice_p, mean_sst_p, local_range, loc_radius, &
       sst_exclude_ice, sst_exclude_diff
  USE g_parfe, &
       ONLY: mydim_nod2d, MPI_INTEGER, MPI_SUM, MPIerr
  USE o_mesh, &
       ONLY: nod2d, coord_nod2d
  USE o_param, &
       ONLY: r_earth, rad
  USE g_rotate_grid, &
       ONLY: r2g

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in)    :: step      ! Current time step
  INTEGER, INTENT(inout) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
! Called by: PDAF_letkf_update   (as U_init_dim_obs)
!EOP

! Local variables
  INTEGER :: i, iter_file, i_obs           ! Counters
  INTEGER :: fileid                        ! ID for NetCDF file
  INTEGER :: id_state                      ! ID for state
  INTEGER :: stat(100)                     ! Status for NetCDF functions
  INTEGER :: startv(2),countv(2)           ! Vectors for reading fields
  CHARACTER(len=4)   :: mype_string        ! String for process rank
  REAL, ALLOCATABLE :: ocoord_n2d_p(:,:)   ! PE-local coordinates of observed SST/profile field
  CHARACTER(len=100) :: sst_file = ''      ! Complete name of SST observation file without path
  REAL(4), ALLOCATABLE :: all_sst_p(:)     ! PE-local complete SST field read from file
  REAL, ALLOCATABLE :: obs_sst_p(:)        ! PE-local observed SST field
  REAL, ALLOCATABLE :: obs_sst_error_p(:)  ! PE-local observed SST error
  REAL :: distance_error                   ! Distance between the obs and the equation 
  REAL, ALLOCATABLE :: weight_error(:)     ! weight accounting for the observation error
  REAL, ALLOCATABLE :: ivariance_obs_p(:)  ! PE-local inverse observation error variance
  INTEGER :: status                        ! Status flag for PDAF gather operation
  INTEGER :: cnt_ex_ice_p, cnt_ex_diff_p   ! PE-local counts of excluded points due to ice and difference
  INTEGER :: cnt_ex_ice, cnt_ex_diff       ! Global counts of excluded points due to ice and difference


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  IF (mype_filter==0) &
       WRITE (*,'(a,5x,a)') 'FESOM-PDAF', 'Assimilate SST observations'

  ! set localization radius
  loc_radius(:) = local_range


! **********************************
! *** Read PE-local observations ***
! **********************************

  ! Initialize complete file name
  WRITE(mype_string,'(i4.4)') mype_filter

  sst_file=TRIM(file_sst_prefix)//TRIM(mype_string)//TRIM(file_sst_suffix)

  ! Allocate array
  ALLOCATE(all_sst_p(myDim_nod2D))

  ! Position to read from file
  iter_file = step / delt_obs_ocn

  ! Read SST observation and std
     
  stat(1) = NF_OPEN(TRIM(path_obs_sst)//sst_file, NF_NOWRITE, fileid)        
  stat(2) = NF_INQ_VARID(fileid, 'obs', id_state)
 
  ! *** Read state estimate ***
  startv(2) = iter_file
  countv(2) = 1
  startv(1) = 1
  countv(1) = myDim_nod2D 

  stat(3) = NF_GET_VARA_REAL (fileid, id_state, startv, countv, all_sst_p)

  ! *** close file  ***
  stat(4) = NF_CLOSE(fileid)
 
  ! check status flag
  DO i=1,4
     IF (stat(i).NE.NF_NOERR) WRITE(*,*) &
          'NetCDF error in reading full SST, no.',i, &
          ' file ',sst_file
     stat(i)=0
  END DO

  ! *** Exclude observations if mean_ice is not zero and SST>0 ***

  IF (sst_exclude_ice) THEN
     cnt_ex_ice_p = 0
     DO i = 1, myDim_nod2D
        IF (mean_ice_p(i) > 0.0 .AND.  all_sst_p(i)<999.0 &
             .AND. (all_sst_p(i)>=0.0 .OR. mean_sst_p(i)>=0.0)) THEN
           all_sst_p(i) = 1.0e6
           cnt_ex_ice_p = cnt_ex_ice_p + 1
        END IF
     END DO

     CALL MPI_Allreduce (cnt_ex_ice_p, cnt_ex_ice, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)

     IF (mype_filter==0) &
          WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
          '--- Observations excluded because of ice', cnt_ex_ice
  END IF

  ! *** Set localization radius to zero for grid points with ice ***

  IF (mype_filter == 0) &
       WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
       '--- Set localization radius to zero for points with ice'

  DO i = 1, myDim_nod2D
     IF (mean_ice_p (i) > 0.0) THEN
        loc_radius (i) = 0.0
     END IF
  END DO

  ! *** Exclude observations if difference from ensemble mean is beyond limit SST_EXCLUDE_DIFF ***

  IF (sst_exclude_diff > 0.0) THEN
     cnt_ex_diff_p = 0
     DO i = 1, myDim_nod2D
        IF (ABS(mean_sst_p(i) - all_sst_p(i)) > sst_exclude_diff .AND. all_sst_p(i)<=999.0) THEN
           all_sst_p(i) = 1.0e6
           cnt_ex_diff_p = cnt_ex_diff_p+1
        END IF
     END DO

     CALL MPI_Allreduce (cnt_ex_diff_p, cnt_ex_diff, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)

     IF (mype_filter==0) &
          WRITE (*,'(a,5x,a,f6.2,a,i7)') 'FESOM-PDAF', &
          '--- Observations excluded due to difference >',sst_exclude_diff,'degC:', cnt_ex_diff
  END IF


! ***************************************************
! *** Count available observations                ***
! *** and initialize index and coordinate arrays. ***
! ***************************************************

  ! *** Count PE-local number of observations ***
  dim_obs_p = 0
  DO i = 1, myDim_nod2d
     IF (ABS(all_sst_p(i)) < 999.0) dim_obs_p=dim_obs_p+1
  ENDDO

  ! *** Initialize index vector of observed surface nodes ***
  IF (ALLOCATED(id_observed_sst_p)) THEN
     DEALLOCATE(id_observed_sst_p)
  END IF
  ALLOCATE(id_observed_sst_p(dim_obs_p))

  i_obs=0
  DO i = 1, myDim_nod2d
     IF (ABS(all_sst_p(i)) < 999.0) THEN
        i_obs = i_obs + 1
        id_observed_sst_p(i_obs) = i
     END IF
  ENDDO

  ! *** Initialize PE-local vector of observations ***
  ALLOCATE(obs_sst_p(dim_obs_p))
  DO i = 1, dim_obs_p
     obs_sst_p(i) = REAL(all_sst_p(id_observed_sst_p(i)), 8)
  ENDDO

  ! *** Initialize coordinate arrays for PE-local observations
  ALLOCATE(ocoord_n2d_p(2, dim_obs_p))
  DO i = 1, dim_obs_p
     ! Rotate to geographic coordinates and store
     CALL r2g(ocoord_n2d_p(1, i), ocoord_n2d_p(2, i), &
          coord_nod2d(1, id_observed_sst_p(i)), coord_nod2d(2, id_observed_sst_p(i)))
  ENDDO


! ***************************************
! *** Define local observation errors ***
! ***************************************

  ALLOCATE(obs_sst_error_p(dim_obs_p))

  ! *** Set constant error *** 
  obs_sst_error_p(:) = rms_obs 

  ! Set inverse observation error variance
  ALLOCATE(ivariance_obs_p(dim_obs_p))
  DO i = 1, dim_obs_p
     ivariance_obs_p(i) = 1.0 / rms_obs**2
  END DO


! *************************************************
! *** Gather global information on observations ***
! *************************************************

  ! *** Initialize global dimension of observation vector ***
  CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

  ! *** Gather full observation vectors ***

  ! Allocate global observation arrays
  ! The arrays are deallocated in prepoststep
  ALLOCATE(obs_sst(dim_obs_f))
  ALLOCATE(obs_sst_error(dim_obs_f))
  ALLOCATE(ivariance_obs(dim_obs_f))
 
  CALL PDAF_gather_obs_f(obs_sst_p, obs_sst, status)

  CALL PDAF_gather_obs_f(obs_sst_error_p, obs_sst_error, status)

  CALL PDAF_gather_obs_f(ivariance_obs_p, ivariance_obs, status)
 

  ! *** Gather global coordinate information
  ! The array is deallocated in prepoststep
  ALLOCATE(ocoord_n2d(2, dim_obs_f))

  CALL PDAF_gather_obs_f2(ocoord_n2d_p, ocoord_n2d, 2, status)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(all_sst_p, obs_sst_p)
  DEALLOCATE(obs_sst_error_p)
  DEALLOCATE(ocoord_n2d_p, ivariance_obs_p)

END SUBROUTINE init_dim_obs_f_sst_pdaf

