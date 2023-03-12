!$Id: init_ens_rnd.F90 1136 2011-08-25 08:23:35Z lnerger $
!BOP
!
! !ROUTINE: init_ens_rnd --- Initialize ensemble by random drawing from trajectory
!
! !INTERFACE:
SUBROUTINE init_ens_rnd(dim, dim_ens, state, ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK):
!
! The routine is called by init_seik. It 
! initializes an ensemble of dim\_ens states
! by random drawing from the long trajectory
! representing the true state.
!
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2010-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_model, &
       ONLY: step_null
  USE output_netcdf, &
       ONLY: file_state
  USE mod_assimilation, &
       ONLY: seedset

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: init_seik
! Calls: memcount
!EOP

! *** local variables ***
  INTEGER :: i, s                     ! counters
  INTEGER :: dim_file                 ! State dimension in file
  INTEGER :: nsteps_file              ! Number of time steps stored in file
  INTEGER :: stat(1000)               ! Array for status flag
  INTEGER :: fileid                   ! ID for NetCDF file
  INTEGER :: id_state                 ! ID for state
  INTEGER :: id_dim                   ! ID for dimension
  INTEGER, ALLOCATABLE :: rnd_indx(:) ! Array of random indices for reading from file
  REAL, ALLOCATABLE :: rnd_array(:)   ! Array of random values
  INTEGER :: max_indx                 ! Maximum number of indices considered in file
  INTEGER :: offset_indx              ! Index-offset in file (=step_null)
  INTEGER :: iseed(4)                 ! Seed-array for DLARNV
  INTEGER :: pos(2)                   ! Position index for reading from file
  INTEGER :: cnt(2)                   ! Count index for reading from file


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(9x, a)') '--- generate ensemble'
  WRITE (*, '(9x, a)') &
       '--- Use random drawing from true trajectory'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
  
  WRITE(*,'(9x,a,a)') '--- Reading true states from ', TRIM(file_state)

  s = 1
  stat(s) = NF_OPEN(file_state, NF_NOWRITE, fileid)

  ! Read size of state vector
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_file)

  ! Read number of time steps stored in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'timesteps', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, nsteps_file)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from init file, no.', i
  END DO


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

! *** Determine index interval to scale random values

  ! Maximum number of indices
  max_indx = nsteps_file - step_null

  ! Index offset in file
  offset_indx = step_null
  WRITE (*, '(9x,a,i7,a,i5)') &
       '--- Select from ',max_indx,' time steps starting at step ',offset_indx
  WRITE (*, '(9x,a,i3)') &
       '--- Use seed set no. ', seedset

! *** Initialize array of random indices

  allocate(rnd_indx(dim_ens))
  allocate(rnd_array(dim_ens))
  ! count allocated memory
  CALL memcount(2, 'r', dim_ens)
  CALL memcount(2, 'i', dim_ens)

  ! Initialize seed
  IF (seedset == 2) THEN
     iseed(1)=1
     iseed(2)=5
     iseed(3)=7
     iseed(4)=9
  ELSE IF (seedset == 3) THEN
     iseed(1)=2
     iseed(2)=5
     iseed(3)=7
     iseed(4)=9
  ELSE IF (seedset == 4) THEN
     iseed(1)=1
     iseed(2)=6
     iseed(3)=7
     iseed(4)=9
  ELSE IF (seedset == 5) THEN
     iseed(1)=1
     iseed(2)=5
     iseed(3)=8
     iseed(4)=9
  ELSE IF (seedset == 6) THEN
     iseed(1)=2
     iseed(2)=5
     iseed(3)=8
     iseed(4)=9
  ELSE IF (seedset == 7) THEN
     iseed(1)=2
     iseed(2)=6
     iseed(3)=8
     iseed(4)=9
  ELSE IF (seedset == 8) THEN
     iseed(1)=2
     iseed(2)=6
     iseed(3)=8
     iseed(4)=11
  ELSE IF (seedset == 9) THEN
     iseed(1)=3
     iseed(2)=6
     iseed(3)=8
     iseed(4)=11
  ELSE IF (seedset == 10) THEN
     iseed(1)=3
     iseed(2)=7
     iseed(3)=8
     iseed(4)=11
  ELSE
     iseed(1)=2*300+1
     iseed(2)=2*150+5
     iseed(3)=2*13+7
     iseed(4)=2*31+9
  END IF

  ! Generate uniformly distributed random values
  CALL dlarnv(1, iseed, dim_ens, rnd_array)

  ! Scale and Transform random values for state selection
  DO i = 1, dim_ens
     rnd_indx(i) = INT(rnd_array(i) * max_indx)
  END DO

! *** Read states according to rnd_indx to initialize ensemble

  s = 1
  stat(s) = NF_INQ_VARID(fileid, 'state', id_state)

  DO i = 1, dim_ens
     
     pos(2) = rnd_indx(i)
     cnt(2) = 1
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_state, pos, cnt, ens(:, i))

  END DO
  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading states from init file, no.', i
  END DO


! ****************
! *** clean up ***
! ****************

  stat(s) = nf_close(fileid)

  deallocate(rnd_indx, rnd_array)

END SUBROUTINE init_ens_rnd
