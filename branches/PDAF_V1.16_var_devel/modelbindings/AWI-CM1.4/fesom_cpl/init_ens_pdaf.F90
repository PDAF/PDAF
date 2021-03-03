!$Id: init_ens_pdaf.F90 2271 2020-04-08 13:04:09Z lnerger $
!>  Routine to initialize ensemble for PDAF
!!
!! User-supplied routine for PDAF.
!!
!! If only a single filter algorithm is used, the 
!! ensemble initialization can be performed directly
!! in this routine. If a single filter is implemented,
!! one can perform the initialization directly here.
!!
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: file_init, path_init, read_inistate, file_inistate, &
       varscale, off_fields_p
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_filter, COMM_filter_fesom, abort_parallel
  USE g_parfe, &
       ONLY: MPI_DOUBLE_PRECISION, MPIerr

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                     !< Process-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)            !< Process-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)     !< Process-local state ensemble
  INTEGER, INTENT(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  INTEGER :: i, member, s, row, col    ! Counters
  INTEGER :: rank                      ! Rank stored in init file
  INTEGER :: dim_p_file                ! Local state dimension in init file
  INTEGER :: fileid                    ! ID for NetCDF file
  INTEGER :: id_dim                    ! ID for dimension
  INTEGER :: id_state,id_svals,id_eof  ! IDs for fields
  INTEGER :: startv(2),countv(2)       ! Vectors for reading fields
  INTEGER :: stat(7)                   ! Status flag for NetCDF commands
  REAL :: fac                          ! Square-root of dim_ens or dim_ens-1
  REAL,ALLOCATABLE :: eof_p(:,:)       ! Matrix of eigenvectors of covariance matrix
  REAL,ALLOCATABLE :: svals(:)         ! Singular values
  REAL,ALLOCATABLE :: omega(:,:)       ! Transformation matrix Omega
  CHARACTER(len=4)   :: mype_string    ! String for process rank
  CHARACTER(len=150) :: file           ! File holding initial state estimate
  INTEGER :: dim_p_read                ! state dimension to be read from file


! **********************
! *** INITIALIZATION ***
! **********************

  dim_p_read = dim_p

  mype0: IF (mype_filter == 0) THEN
    WRITE (*, '(/a, 8x,a)') 'FESOM-PDAF', 'Generate state ensemble from covariance matrix'
    WRITE (*, '(a, 8x,a)') &
         'FESOM-PDAF', '--- use 2nd order exact sampling (SEIK type)'
    WRITE (*, '(a, 8x,a,i5)') 'FESOM-PDAF', '--- number of EOFs:',dim_ens-1
  END IF mype0

  ! allocate memory for temporary fields
  ALLOCATE(eof_p(dim_p, dim_ens-1))
  ALLOCATE(svals(dim_ens-1))
  ALLOCATE(omega(dim_ens, dim_ens-1))
  

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
  
  IF (mype_filter == 0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- Read initial state'

  write(mype_string,'(i4.4)') mype_filter

  file=Trim(path_init)//Trim(file_init)//TRIM(mype_string)//'.nc'

  s = 1
  stat(s) = NF_OPEN(file, NF_NOWRITE, fileid)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) THEN
        WRITE(*, *) 'NetCDF error in opening initialization file, no.', i
        STOP
     END IF
  END DO

  ! Read size of state vector
  s = 1
  stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_file)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from file, no.', i
  END DO


  checkdim: IF (dim_p_read == dim_p_file .AND. rank >= dim_ens-1) THEN

     IF (mype_filter == 0) WRITE (*,'(8x,a)') '--- Read covariance matrix'

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF_INQ_VARID(fileid, 'running_meanstate', id_state)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

     ! Read initialization information
     s = s + 1
     stat(s) = NF_GET_VAR_DOUBLE(fileid, id_state, state_p(1:dim_p_read))

     startv(2) = 1
     countv(2) = dim_ens-1
     startv(1) = 1
     countv(1) = dim_p_read
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_eof, startv, countv, eof_p)

     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_svals, 1, dim_ens-1, svals)

     s = s + 1
     stat(s) = nf_close(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
     END DO

     IF (mype_filter==0) THEN
        WRITE(*,*) 'svals', svals
     END IF


! **************************************************
! *** Initialize initial state from model fields ***
! **************************************************

     CALL collect_state_pdaf(dim_p, state_p)


! ********************************
! *** Initialize initial state ***
! ********************************

     readinistate: IF (read_inistate) THEN
        IF (mype_filter == 0) WRITE (*,'(a,8x,a)') 'FESOM-PDAF','--- Read separate initial state'

        write(mype_string,'(i4.4)') mype_filter

        file=Trim(path_init)//Trim(file_inistate)//TRIM(mype_string)//'.nc'

        s = 1
        stat(s) = NF_OPEN(file, NF_NOWRITE, fileid)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) THEN
              WRITE(*, *) 'NetCDF error in opening initialization file, no.', i
              STOP
           END IF
        END DO

        ! Read size of state vector
        s = 1
        stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
        s = s + 1
        stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_file)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading dimensions from file, no.', i
        END DO

        ! Inquire IDs for mean state
        s = 1
        stat(s) = NF_INQ_VARID(fileid, 'meanstate', id_state)

        ! Read initialization information
        s = s + 1
        stat(s) = NF_GET_VAR_DOUBLE(fileid, id_state, state_p(1:dim_p_read))

        s = s + 1
        stat(s) = nf_close(fileid)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
        END DO

     END IF readinistate


! ******************************************
! *** Set variance of ice fields to zero ***
! ******************************************

     if (mype_filter==0) write (*,*) 'RESET VARIANCE OF ICE FIELDS TO ZERO!'

     s = 0
     DO member = 1, dim_ens-1
        DO i = 1+off_fields_p(7), dim_p
           eof_p(i, member) = 0.0
           s = s+1
        END DO
     END DO


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

     IF (dim_ens>1) THEN
        ! Only initialize Omega if ensemble size > 0

        IF (mype_filter==0) THEN

           WRITE (*,'(a,8x,a)') 'FESOM-PDAF','--- generate state ensemble'

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_omega(dim_ens-1, Omega, 1, 1)

            ! ***      Generate ensemble of states         ***
            ! *** x_i = x + sqrt(FAC) eofV (Omega C^(-1))t ***

            ! A = Omega C^(-1)
           DO col = 1, dim_ens-1
              DO row = 1, dim_ens
                 Omega(row, col) = Omega(row,col) * svals(col)
              END DO
           END DO
        END IF
        CALL MPI_Bcast(Omega, dim_ens*(dim_ens-1), MPI_DOUBLE_PRECISION, 0, &
             COMM_filter_fesom, MPIerr)
     END IF

     ! state_ens = state+ sqrt(dim_ens-1) eofV A^T

     DO col = 1,dim_ens
        ens_p(1:dim_p,col) = state_p(1:dim_p)
     END DO

     IF (dim_ens>1) THEN
        ! Only add perturbations if ensemble size > 0

        fac = varscale * SQRT(REAL(dim_ens-1))

        CALL DGEMM('n', 't', dim_p, dim_ens, dim_ens-1, &
             fac, eof_p, dim_p, Omega, dim_ens, 1.0, ens_p, dim_p)
     END IF

  ELSE checkdim

      ! *** Rank stored in file is smaller than requested EOF rank ***
     WRITE(*,*) 'FESOM-PDAF: ','ERROR: Rank stored in file is smaller than requested EOF rank'
     CALL abort_parallel()

  END IF checkdim


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eof_p, omega)

END SUBROUTINE init_ens_pdaf
  
