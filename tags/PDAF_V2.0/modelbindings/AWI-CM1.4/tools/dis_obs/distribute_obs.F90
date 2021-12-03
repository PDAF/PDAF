! $Id: distribute_obs.F90 644 2012-06-20 09:27:59Z lnerger $
!BOP
!
! !Program: distribute_obs --- Generate distributed observation files
!
! !INTERFACE:
PROGRAM distribute_obs

! !DESCRIPTION:
! This program generates files holding distributed observation
! information according to the distribution information of the
! mesh. The observation infromation is read from a NetCDF file 
! holding the global observation information.

! !USES:
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
!EOP

! Local variables
  INTEGER :: i, s, iter, pe
  LOGICAL :: twin_data
  CHARACTER(len=120) :: inpath, outpath, distpath
  CHARACTER(len=120) :: infile, outfile, partfile, mylistfile
  CHARACTER(len=150) :: ncfile_in, ncfile_out
  CHARACTER(len=150) :: attstr
  CHARACTER(len=10) :: mype_string
  INTEGER :: ncid_in, ncid_out
  INTEGER :: id_dim, id_stderr, id_time, id_list
  INTEGER :: id_sst, id_my_sst
  INTEGER :: dimid_iter, dimid_one, dimid_n2d
  INTEGER :: varid_time
  INTEGER :: n2d, n3d, steps
  INTEGER :: idummy
  INTEGER :: sum_n2d, sum_n3d
  INTEGER :: stat(100)
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)
  INTEGER :: npes
  INTEGER, ALLOCATABLE :: n2d_dist(:), n3d_dist(:)
  INTEGER, ALLOCATABLE :: mylist_n2d(:)
  REAL, ALLOCATABLE :: times(:)
  REAL, ALLOCATABLE :: sst(:), my_sst(:)
  REAL :: stderr_sst
  INTEGER, ALLOCATABLE :: cntsst(:), cntsstday(:)


! ************************************************
! *** Configuration                            ***
! ************************************************

  ! Path to and name of file holding global model fields
  inpath = 'fesom_obs/'
  infile = 'obs_SST.nc'

  ! Path to mesh
  distpath = '../input/CORE2_final/dist/384p/'
  partfile = 'rpart.out'
  mylistfile = 'my_list'

  ! Path to and name stub of output files
  outpath = 'fesom_obs/dist384/'
  outfile = 'obs_SST_dis384'

  ! Whether the observation file hold twin data 
  ! only in this case the estimated stderr_sst and time array are handled
  twin_data = .FALSE.


! ************************************************
! *** Init                                     ***
! ************************************************

  WRITE (*,'(3x,a)') '********************************************************'
  WRITE (*,'(3x,a)') '*** Generate distributed observation files for FESOM ***'
  WRITE (*,'(3x,a/)') '********************************************************'

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
  s = s + 1

  ! Get dimensions
  stat(s) = NF_INQ_DIMID(ncid_in, 'nodes_2D', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, n2d)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'time', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, steps)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO

  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  WRITE (*,'(5x,1x,a3,i10)') 'n2d', n2d
  WRITE (*,'(5x,a4,i10)') 'iter', steps


! ***********************************
! *** Read mesh partitioning data ***
! ***********************************

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
  
  IF (sum_n2d /= n2d) THEN
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

  ! Read stderr and time information from input file

  ALLOCATE(times(steps))

  IF (twin_data) THEN
    s = 1
    stat(s) = NF_INQ_VARID(ncid_in, 'time', id_time)
    s = s + 1
    stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_time, times)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'std', id_stderr)
    s = s + 1
    stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_stderr, stderr_sst)

    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in reading from input file, no.', i
    END DO
  ENDIF

  ! Initialize one file for each PE
  initoutfiles: DO pe = 0, npes - 1

     WRITE(mype_string,'(i4.4)') pe

     s = 1
     stat(s) = NF_CREATE(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', 0, ncid_out) 

     attstr  = 'daily SST observations from year 2016'
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
          TRIM(attstr)) 

     ! Define dimensions
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'iter', steps, dimid_iter)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'one', 1, dimid_one)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'nodes_2D', n2d_dist(pe + 1), dimid_n2d)

     ! Define variables
     IF (twin_data) THEN
        s = s + 1
        stat(s) = NF_DEF_VAR(ncid_out, 'stderr_sst', NF_DOUBLE, 1, dimid_one, Id_stderr) 
        s = s + 1
        stat(s) = NF_DEF_VAR(ncid_out, 'time', NF_DOUBLE, 1, dimid_iter, Id_time) 
     ENDIF
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'nodlist', NF_INT, 1, dimid_n2d, id_list)

     dimids(1) = DimId_n2d
     dimids(2) = dimid_iter

     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'obs', NF_DOUBLE, 2, dimids(1:2), id_my_sst)
     s = s + 1
     stat(s) = NF_ENDDEF(ncid_out) 

     IF (twin_data) THEN
        ! Write std error
        s = s + 1
        stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, Id_stderr, stderr_sst)

        ! Write times
        s = s + 1
        stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, Id_time, times)
     ENDIF

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in init of output file, no.', i
     END DO

  END DO initoutfiles


! ***********************************
! *** Generate distributed fields ***
! ***********************************

  WRITE (*,'(/1x,a/)') '------- Generate distributed files -------------'

  ! Allocate global sst field
  ALLOCATE(sst(n2d))
  ALLOCATE(cntsst(steps))
  ALLOCATE(cntsstday(npes))
  cntsst = 0
  cntsstday = 0

  peloop: DO pe = 0, npes - 1

     ! allocate fields
     ALLOCATE(mylist_n2d(n2d_dist(pe+1)))
     ALLOCATE(my_sst(n2d_dist(pe+1)))
     
     ! Read mylist
     WRITE(mype_string,'(i4.4)') pe

     OPEN(unit=10,file=TRIM(distpath)//TRIM(mylistfile)//TRIM(mype_string)//'.out', &
          status='old', form='formatted')

     READ(10,*) idummy
     READ(10,*) idummy
     READ(10,*) idummy
     READ(10,*) mylist_n2d
     CLOSE(10)

     ! Open output file
     s = 1
     stat(s) = NF_OPEN(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_WRITE, ncid_out)

     ! Write MyList
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, mylist_n2d)

     ! Get ID for global sst field
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'sst', id_sst)

     ! Get ID for local sst field
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'obs', id_my_sst)

     ! Loop through time slices and generate and write local fields
     stepsloop: DO iter = 1, steps

        ! Read global sst field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = n2d  
        stat(1) = NF_GET_VARA_DOUBLE(ncid_in, id_sst, startv, countv, sst)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading SSH'

        ! Select local nodes
        DO i = 1, n2d_dist(pe+1)
           my_sst(i) = sst(mylist_n2d(i))
           IF (my_sst(i)>-100.0 .AND. my_sst(i)<100.0) cntsst(iter) = cntsst(iter)+1
           IF (iter==1) THEN
              IF (my_sst(i)>-100.0 .AND. my_sst(i)<100.0) cntsstday(pe+1) = cntsstday(pe+1)+1
           END IF
        END DO

        ! Write local sst field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = n2d_dist(pe+1)
        stat(1) = NF_PUT_VARA_DOUBLE(ncid_out, id_my_sst, startv, countv, my_sst)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in writing local SSH'

     END DO stepsloop

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in writing to output file, no.', i
     END DO

     DEALLOCATE(mylist_n2d, my_sst)

  END DO peloop

  stat(1) = nf_close(ncid_in)

  DO iter=1,steps
     WRITE (*,*) 'step, cntsst', iter, cntsst(iter)
  END DO
  DO pe=1,npes
     write (*,*) 'pe, cntsstday1', pe-1, cntsstday(pe)
  END DO

  WRITE (*,'(1x,a/)') '------- END -------------'

  DEALLOCATE(cntsst, cntsstday)

END PROGRAM distribute_obs
