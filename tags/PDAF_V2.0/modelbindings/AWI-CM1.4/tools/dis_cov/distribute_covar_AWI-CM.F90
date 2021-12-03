! $Id: distribute_covar.F90 1604 2016-05-30 06:42:16Z lnerger $
!BOP
!Modified by QT 2018-02-27 ---- for AWI-CM 

! !Program: distribute_covar --- Generate distributed covariance matrix files
!
! !INTERFACE:
PROGRAM distribute_covar

! !DESCRIPTION:
! This program generates files holding distributed covariance matrix
! information according to the distribution information of the
! mesh. The global covariance matrix is read from a NetCDF file.
! The input file contains values in single precision. For compatibility
! with FESOM the distributed output files will be in double precision.

! !USES:
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
!EOP

  ! Local variables
  INTEGER :: i, s, iter, pe, irank
  CHARACTER(len=120) :: inpath, outpath, distpath
  CHARACTER(len=120) :: infile, outfile, partfile, mylistfile
  CHARACTER(len=150) :: ncfile_in, ncfile_out
  CHARACTER(len=10) :: mype_string
  CHARACTER(len=150) :: title
  CHARACTER(len=150) :: fieldsstr
  INTEGER :: ncid_in, ncid_out
  INTEGER :: id_dim
  INTEGER :: id_sigma, id_list
  INTEGER :: id_mean, id_svec
  INTEGER :: dimid_rank, dimid_one, dimid_n2d, dimid_n3d, dimid_nfields, dimid_state
  INTEGER :: id_mu, id_mv, id_mw, id_mz, id_ms, id_mt
  INTEGER :: id_maice, id_mhice, id_mhsnow, id_muice, id_mvice
  INTEGER :: id_svdu, id_svdv, id_svdw, id_svdz, id_svds, id_svdt
  INTEGER :: id_svdaice, id_svdhice, id_svdhsnow, id_svduice, id_svdvice
  INTEGER :: id_u, id_v, id_w, id_z, id_t, id_s, id_out
  INTEGER :: id_aice, id_hice, id_hsnow, id_uice, id_vice
  INTEGER :: startpos
  INTEGER :: n2d, n3d, rank, dim_state, nfields
  INTEGER :: edim_n2d, edim_n3d
  INTEGER :: offsets(6)     ! Field offsets in state vector
  INTEGER :: idummy
  INTEGER :: sum_n2d, sum_n3d
  INTEGER :: stat(100)
  INTEGER :: dim_state_l
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)
  INTEGER :: npes
  INTEGER, ALLOCATABLE :: n2d_dist(:), n3d_dist(:)
  REAL(kind=4), ALLOCATABLE :: svals(:)
  REAL(kind=4), ALLOCATABLE :: field_n2d(:), field_n3d(:)
  REAL(kind=4), ALLOCATABLE :: field_myn2d(:), field_myn3d(:)
  REAL(kind=8), ALLOCATABLE :: state_l(:)
  INTEGER, ALLOCATABLE :: mylist_n2d(:), mylist_n3d(:)


! ************************************************
! *** Configuration                            ***
! ************************************************

  ! Path to and name of file holding global covariance matrix
  inpath = '/gfs1/work/hbkqtang/output/cov_6/'
  infile = 'cov.nc'

  ! Path to mesh
  distpath = '/gfs1/work/hbkqtang/input/CORE2_final/dist/384p/'
  partfile = 'rpart.out'
  mylistfile = 'my_list'

  ! Path to and name stub of output files
  outpath = '/gfs2/work/hbkqtang/output/dis_cov/'
  outfile = 'covar'


! ************************************************
! *** Init                                     ***
! ************************************************

  WRITE (*,'(3x,a)') '**************************************************************'
  WRITE (*,'(3x,a)') '*** Generate distributed covariance matrix files for AWi-CM***'
  WRITE (*,'(3x,a/)') '**************************************************************'

  ncfile_in = TRIM(inpath)//TRIM(infile)
  WRITE (*,*) 'Read fields from file: ',TRIM(ncfile_in)

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  WRITE (*,*) 'Write distributed fields to files: ',TRIM(ncfile_out),'_XXXX.nc'
  WRITE (*,*) 'Read partitioning information from directory: ',TRIM(distpath)


! *******************************************
! *** Open model file and read dimensions ***
! *******************************************

  s = 1
  stat(s) = NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in)

  ! Get dimensions
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, rank)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nodes_2D', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, n2d)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nodes_3D', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, n3d)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_state)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nfields', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, nfields)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO

  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  WRITE (*,'(10x,1x,a3,i12)') 'n2d', n2d
  WRITE (*,'(10x,1x,a3,i12)') 'n3d', n3d
  WRITE (*,'(10x,a4,i12)') 'rank', rank
  WRITE (*,'(5x,a,i12)') 'dim_state', dim_state
  WRITE (*,'(5x,2x,a,i12)') 'nfields', nfields


! ***********************************
! *** Read mesh partitioning data ***
! ***********************************
! 
  OPEN(unit=10,file=TRIM(distpath)//TRIM(partfile), status='old', form='formatted')

  READ(10,*) npes

  ALLOCATE(n2d_dist(npes),n3d_dist(npes))
  READ(10,*) n2d_dist
  READ(10,*) n3d_dist
  CLOSE(10)

  sum_n2d = 0
  sum_n3d = 0
  DO i= 1, npes
     sum_n2d = sum_n2d + n2d_dist(i)
     sum_n3d = sum_n3d + n3d_dist(i)
  END DO

  IF (sum_n2d /= n2d .OR. sum_n3d /= n3d) THEN
     WRITE (*,*) 'Partitioned mesh not consistent with global mesh'
     STOP
  END IF

  WRITE (*,'(/1x,a)') 'Mesh distribution:'
  WRITE (*,'(5x,a,i6)') 'Number of PEs:', npes
  WRITE (*,'(5x,a,4x,a)') 'Local nodes:   nod2d','nod3d'
  DO i = 1, npes
     WRITE (*,'(18x,i7,i10)') n2d_dist(i), n3d_dist(i)
  END DO


! ************************************
! *** Initialize distributed files ***
! ************************************

  WRITE (*,'(/1x,a/)') '------- Initialize distributed files -------------'

  ! Read singular values from input file

  ALLOCATE(svals(rank))

  s = 1
  stat(s) = NF_INQ_VARID(ncid_in, 'sigma', id_sigma)
  s = s + 1
  stat(s) = NF_GET_VAR_REAL(ncid_in, id_sigma, svals)
  s = s + 1
  title = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'title', title)
  s = s + 1
  fieldsstr = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'state_fields', fieldsstr)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading from input file, no.', i
  END DO

  ! Initialize one file for each PE
  initoutfiles: DO pe = 0, npes - 1

     dim_state_l = n2d_dist(pe + 1) + 5 * n3d_dist(pe + 1)

     WRITE(mype_string,'(i4.4)') pe
     WRITE (*,*) 'pe and dim_state_l', pe, dim_state_l

     s = 1
     !stat(s) = NF_CREATE(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', IOR(NF_NOCLOBBER,NF_64BIT_OFFSET), ncid_out) 
     stat(s) = NF_CREATE(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', 0, ncid_out)
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(title), &
          TRIM(title)) 
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'state_fields', LEN_TRIM(fieldsstr), &
          TRIM(fieldsstr)) 

     ! Define dimensions
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'rank',  rank, dimid_rank)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'nodes_2D', n2d_dist(pe + 1), dimid_n2d)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'nodes_3D', n3d_dist(pe + 1), dimid_n3d)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'dim_state', dim_state_l, dimid_state)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'one', 1, dimid_one)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'nfields', nfields, dimid_nfields)

     ! Define variables

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'sigma', NF_DOUBLE, 1, dimid_rank, Id_sigma)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_n2d', NF_INT, 1, dimid_n2d, id_list)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_n3d', NF_INT, 1, dimid_n3d, id_list)

     ! running mean state (for the last snap shot)
     dimids(1) = DimId_state
     dimids(2) = dimid_one

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'running_meanstate', NF_DOUBLE, 2, dimids, Id_mean)

     ! singular state vectors
     dimids(1) = DimId_state
     dimids(2) = dimid_rank

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'V', NF_DOUBLE, 2, dimids, Id_svec)

     ! End define mode
     s = s + 1
     stat(s)=NF_ENDDEF(ncid_out)

     ! Write singular values
     s = s + 1
     stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, id_sigma, REAL(svals(1:rank),8))

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in init of output file, no.', i
     END DO

  END DO initoutfiles


! ****************************************************
! *** Generate and write distributed state vectors ***
! ****************************************************

  WRITE (*,'(1x,a/)') '------- Generate and write distributed files -------------'

  ! Allocate global 2d and 3d fields
  ALLOCATE(field_n2d(n2d))
  ALLOCATE(field_n3d(n3d))

  ! *** Read, distribute and write mean state

  peloop: DO pe = 0, npes - 1

     WRITE (*,'(5x, a, i5)') '-- Process ', pe

     ! Size of local state vector
     dim_state_l = n2d_dist(pe + 1) + 5 * n3d_dist(pe + 1)

     ! allocate fields
     ALLOCATE(field_myn2d(n2d_dist(pe+1)))
     ALLOCATE(field_myn3d(n3d_dist(pe+1)))
     ALLOCATE(state_l(dim_state_l))

     ! define offsets
     offsets (1) = 0                                       ! offset of field SSH
     offsets (2) = n2d_dist(pe + 1)                        ! offset of field u
     offsets (3) = n2d_dist(pe + 1) + n3d_dist(pe + 1)     ! offset of field v
     offsets (4) = n2d_dist(pe + 1) + 2 * n3d_dist(pe + 1) ! offset of field w
     offsets (5) = n2d_dist(pe + 1) + 3 * n3d_dist(pe + 1) ! offset of field t
     offsets (6) = n2d_dist(pe + 1) + 4 * n3d_dist(pe + 1) ! offset of field s

     ! Read mylist
     WRITE(mype_string,'(i4.4)') pe

     OPEN(unit=10,file=TRIM(distpath)//TRIM(mylistfile)//TRIM(mype_string)//'.out', &
          status='old', form='formatted')
     READ(10, *) idummy
     READ(10, *) idummy
     READ(10, *) edim_n2d
     ALLOCATE(mylist_n2d(n2d_dist(pe + 1) + edim_n2d))
     READ(10, *) mylist_n2d
     READ(10, *) idummy
     READ(10, *) edim_n3d
     ALLOCATE(mylist_n3d(n3d_dist(pe + 1) + edim_n3d))
     READ(10, *) mylist_n3d
     CLOSE(10)

     ! Open output file
     s = 1
     stat(s) = NF_OPEN(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_WRITE, ncid_out)

     ! Write MyLists
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_n2d', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, mylist_n2d)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_n3d', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, mylist_n3d)

     ! Inquire IDs for mean state and singular vectors
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'running_meanstate', id_mean)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'V', id_svec)

     ! Get IDs for global fields from input file
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'ssh_mean', id_mz)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'u_mean', id_mu)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'v_mean', id_mv)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'wpot_mean', id_mw)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'temp_mean', id_mt)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'salt_mean', id_ms)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'ssh_svd', id_svdz)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'u_svd', id_svdu)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'v_svd', id_svdv)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'w_svd', id_svdw)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'temp_svd', id_svdt)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'salt_svd', id_svds)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in inquire from output file, no.', i, ' pe ',pe
     END DO

     loopfields: DO irank = 0, rank

        IF (irank == 0) THEN
           ! Treat mean state

           ! Ids for input
           id_z = id_mz
           id_u = id_mu
           id_v = id_mv
           id_w = id_mw
           id_t = id_mt
           id_s = id_ms

           ! Id for output
           id_out = id_mean
           ! Column of matrix in file
           startpos = 1
        ELSE
           ! Treat singular vectors

           ! Ids for input
           id_z = id_svdz
           id_u = id_svdu
           id_v = id_svdv
           id_w = id_svdw
           id_t = id_svdt
           id_s = id_svds

           ! Id for output
           id_out = id_svec
           ! Column of matrix in file
           startpos = irank
        END IF

        ! *** Generate local state vector

        ! Read global SSH
        startv(2) = startpos
        countv(2) = 1
        startv(1) = 1
        countv(1) = n2d  
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_z, startv, countv, field_n2d)

        ! Initialize local nodes of local state vector
        DO i = 1, n2d_dist(pe+1)
           state_l(i + offsets (1)) = field_n2d(mylist_n2d(i))
        END DO

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading from input file, no.', i, &
                ' pe, irank',pe, irank
        END DO
        ! Read global U
        countv(1) = n3d  
        s = 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_u, startv, countv, field_n3d)

        ! Initialize local nodes of local state vector
        DO i = 1, n3d_dist(pe+1)
           state_l(i + offsets (2)) = REAL(field_n3d(mylist_n3d(i)),8)
        END DO

        ! Read global V
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_v, startv, countv, field_n3d)

        ! Initialize local nodes of local state vector
        DO i = 1, n3d_dist(pe+1)
           state_l(i + offsets (3)) = REAL(field_n3d(mylist_n3d(i)),8)
        END DO

        ! Read global W
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_w, startv, countv, field_n3d)

        ! Initialize local nodes of local state vector
        DO i = 1, n3d_dist(pe+1)
           state_l(i + offsets (4)) = REAL(field_n3d(mylist_n3d(i)),8)
        END DO

        ! Read global Temperature
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_t, startv, countv, field_n3d)

        ! Initialize local nodes of local state vector
        DO i = 1, n3d_dist(pe+1)
           state_l(i + offsets (5)) = REAL(field_n3d(mylist_n3d(i)),8)
        END DO

        ! Read global Salinity
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_s, startv, countv, field_n3d)

        ! Initialize local nodes of local state vector
        DO i = 1, n3d_dist(pe+1)
           state_l(i + offsets (6)) = REAL(field_n3d(mylist_n3d(i)),8)
        END DO

        ! Write local state vector
        startv(2) = startpos
        countv(2) = 1
        startv(1) = 1
        countv(1) = dim_state_l
        s = 1
        stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_out, startv, countv, state_l)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in writing to output file, no.', i, &
                ' pe, irank',pe, irank
        END DO

     END DO loopfields

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in closing output file, no.', i, ' pe ',pe
     END DO

     DEALLOCATE(field_myn2d, field_myn3d)
     DEALLOCATE(state_l, mylist_n2d, mylist_n3d)

  END DO peloop

  stat(1) = nf_close(ncid_in)

  WRITE (*,'(/1x,a/)') '------- END -------------'

END PROGRAM distribute_covar
