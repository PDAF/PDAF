! !Program: obs_interp --- Interpolate original SST observations into fesom mesh grid
!
! !INTERFACE:

PROGRAM obs_interp

! This program handles SST observations from Copernicus.
! (Data product: SST_GLO_SST_L3S_NRT_OBSERVATIONS_010_010)
! The observations have a higher resolution than the model grid.
! Here a subper-obbing of the observations is performed by
! averaging them onto the unstructured model grid. The
! prepared observations are stored in a netCDF file.
! Processed are both the SST observation values and the standard 
! deviation.
 
  IMPLICIT NONE
#include "netcdf.inc"

  INTEGER :: i, n, j, q, ind, s, month, iday, cnt_day, i_obs, j_obs      ! Counters
  INTEGER :: nod2D, elem2D, node, rmax
  REAL, ALLOCATABLE :: coord_nod2D (:,:)         
  INTEGER, ALLOCATABLE :: elem2D_nodes (:,:)
  REAL :: x, y, ocean_area  
  CHARACTER(len=120) :: meshpath, obspath, obsfile, ncfile_in
  CHARACTER(len=120) :: outpath, outfile, ncfile_out
  CHARACTER(len=150) :: attstr
  REAL, ALLOCATABLE :: elem_area(:)
  REAL, ALLOCATABLE :: area(:)
  REAL, ALLOCATABLE :: mesh_resolution (:)
  REAL :: radius
  INTEGER, ALLOCATABLE :: nod_in_elem2D_num (:)
  INTEGER, ALLOCATABLE :: nod_in_elem2D (:,:)
  INTEGER :: ncid_in, id_dim, dim_x, dim_y, id_SST, id_std, id_time, id_lon, id_lat
  INTEGER :: ncid_out, dimid_n2d, dimid_time
  INTEGER :: id_interp_SST, id_interp_std
  INTEGER :: elnodes(3), ed(2), elem, nz
  REAL :: a(2), b(2), ax, ay, lon, lat, vol
  REAL, ALLOCATABLE :: work_array(:)
  LOGICAL :: cartesian
  REAL :: cyclic_length     ! [radians]
  REAL, PARAMETER :: r_earth=6367500.0
  REAL, PARAMETER :: pi=3.14159265358979
  INTEGER(2), ALLOCATABLE :: obs_file (:,:)
  INTEGER(1), ALLOCATABLE :: std_file (:,:)
  INTEGER :: cnt 
  INTEGER :: stat(100)
  INTEGER :: countv(3), startv(3), startv_out(2), countv_out(2)
  REAL :: radius_coord (2)
  REAL :: xmin, xmax, ymin, ymax, coord_x, coord_y
  REAL :: mean_obs, mean_std      ! temporay store mean value for each time step
  REAL, ALLOCATABLE :: average_obs (:), average_std (:)        ! store mean value over all 2D nodes
  REAL, ALLOCATABLE :: robs_file (:,:)
  REAL, ALLOCATABLE :: rstd_file (:,:)
  REAL(4), ALLOCATABLE :: obs_lon(:)
  REAL(4), ALLOCATABLE :: obs_lat(:)
  INTEGER :: month_day (12)
  CHARACTER(len=4) :: fileprefix
  CHARACTER(len=102) :: filesuffix
  INTEGER :: dimids(2)
  INTEGER :: day, jmin, jmax, imin, imax 
  REAL, ALLOCATABLE :: sec_in_year(:)
  REAL :: fillvalue
  REAL :: al, be, ga, det
  REAL, PARAMETER :: alphaEuler=50.0   ! Rotation parameter for rotated FESOM mesh
  REAL, PARAMETER :: betaEuler=15.0    ! Rotation parameter for rotated FESOM mesh
  REAL, PARAMETER :: gammaEuler=-90.0  ! Rotation parameter for rotated FESOM mesh
  REAL :: rotate_matrix(3,3)           ! Rotation matrix
  REAL, PARAMETER :: rad=pi/180.0
  REAL :: xr, yr, zr, xg, yg, zg
  INTEGER :: nobs_day
  REAL :: avg_res_factor
  INTEGER :: allcnt_day, mincnt_day, maxcnt_day
  REAL :: obs_min  ! minimum SST value for interpolated observation
  LOGICAL :: write_area


! *********************
! *** Configuration ***
! *********************

  ! Path to the fesom mesh
  meshpath = '../input/CORE2_final/'

  ! Path and name of SST observation file
  obspath = '../observation/SST_CMEMS/'

  ! Path to and name of output file holding observations
  outpath = 'fesom_obs/'
  outfile = 'obs_SST.nc'

  ! Specify observation file
  fileprefix = '2016'
  filesuffix = '-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.nc'

  ! Configuration for the year (2016 is a leap year)
  day = 366
  month_day = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

  ! Resolution factor over which the observations are averaged
  ! used is mesh_resolution(n) * avg_res_factor
  avg_res_factor = 0.5

  ! Define minimum SST which is allowed to use
  obs_min = -3.0

  ! Whether to write the areas into a file
  write_area = .false.


! ********************
! *** Preparations ***
! ********************

  ! Compute the rotation matrix
  al=alphaEuler*rad
  be=betaEuler*rad
  ga=gammaEuler*rad

  ! rotation matrix
  rotate_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
  rotate_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
  rotate_matrix(1,3)=sin(ga)*sin(be)
  rotate_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
  rotate_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
  rotate_matrix(2,3)=cos(ga)*sin(be)
  rotate_matrix(3,1)=sin(be)*sin(al)
  rotate_matrix(3,2)=-sin(be)*cos(al)
  rotate_matrix(3,3)=cos(be)


! ********************************
! *** Read 2D mesh information ***
! ********************************

  ! Read 2D nodes and elements  
  open (20,file=trim(MeshPath)//'nod2d.out',  status='old')
  open (21,file=trim(MeshPath)//'elem2d.out', status='old')
  write(*,*) '2D mesh is opened'

  read(20,*) nod2D 
  allocate(coord_nod2D(2, nod2D))
  do i=1, nod2D
     read(20,*) n, x, y, ind

     ! Store coordinates radians
     coord_nod2D(1, i)=x*rad
     coord_nod2D(2, i)=y*rad
  end do
  close(20)

  read(21,*)  elem2D      
  allocate(elem2D_nodes(3,elem2D))  
  do n=1, elem2D
     read(21,*) elem2D_nodes(:,n)     
  end do
  close(21)
  
  write(*,*) 'read_2Dmesh: DONE'


! *****************************
! *** Compute the mesh area ***
! *****************************

  cartesian = .false.
  cyclic_length = 2.0*pi

  ALLOCATE(elem_area(elem2D))
  ALLOCATE(area(nod2D))   !! Extra size just for simplicity
                                             !! in some further routines
  ALLOCATE(mesh_resolution(nod2D))

  ! Get nodes indiced belonging to element
  ALLOCATE(nod_in_elem2D_num(nod2D))
  nod_in_elem2D_num=0
  DO n=1,elem2D
     DO j=1,3
        node=elem2D_nodes(j,n)
        IF (node>nod2D) CYCLE
        nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
     END DO
  END DO

  ! =========================
  ! The areas of triangles:
  ! =========================

  DO n=1, elem2D
     elnodes=elem2D_nodes(:,n)
     ay=SUM(coord_nod2D(2,elnodes))/3.0
     ay=COS(ay)
     IF (cartesian) ay=1.0
     a=coord_nod2D(:,elnodes(2))-coord_nod2D(:,elnodes(1))
     b=coord_nod2D(:,elnodes(3))-coord_nod2D(:,elnodes(1))
     IF(a(1)>cyclic_length/2.) a(1)=a(1)-cyclic_length
     IF(a(1)<-cyclic_length/2.) a(1)=a(1)+cyclic_length
     IF(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
     IF(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
     a(1)=a(1)*ay
     b(1)=b(1)*ay
     elem_area(n)=0.5*ABS(a(1)*b(2)-b(1)*a(2))
  END DO

  ! =====================================================
  ! Scalar element 
  ! areas at different levels (there can be partly land)
  ! =====================================================

  rmax=maxval(nod_in_elem2D_num(1:nod2D))
  allocate(nod_in_elem2D(rmax,nod2D))
  nod_in_elem2D = 0
  nod_in_elem2D_num=0 
  do n=1,elem2D   
     do j=1,3
        node=elem2D_nodes(j,n)
        if (node>nod2D) cycle 
        nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
        nod_in_elem2D(nod_in_elem2D_num(node),node)=n
     end do
  end do
 
  area=0.0
  DO n=1, nod2D
     DO j=1, nod_in_elem2D_num(n)
        elem=nod_in_elem2D(j,n)
        area(n)=area(n)+elem_area(elem)/3.0
     END DO
  END DO

  IF (write_area) THEN
     open (23,file='area_nod2d.out',  status='replace')
     DO n=1, nod2D
        write (23,*) area(n)
     END DO
  END IF
  close(23)

  ! Only areas through which there is exchange are counted

  ! ===========================
  ! Update to proper dimension
  ! ===========================

  elem_area=elem_area*r_earth*r_earth
  area=area*r_earth*r_earth

  ! coordinates are in radians, edge_dxdy are in meters,
  ! and areas are now in m^2
 
  ALLOCATE(work_array(nod2D))
  mesh_resolution=SQRT(area(:)/pi)*2.0
  DO q=1, 3 !apply mass matrix N times to smooth the field
     DO n=1, nod2D
        vol=0.0
        work_array(n)=0.0
        DO j=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(j, n)
           elnodes=elem2D_nodes(:,elem)
           work_array(n)=work_array(n)+SUM(mesh_resolution(elnodes))/3.0*elem_area(elem)
           vol=vol+elem_area(elem)
        END DO
        work_array(n)=work_array(n)/vol
     END DO
     DO n=1,nod2D
        mesh_resolution(n)=work_array(n)
     ENDDO
  END DO
  DEALLOCATE(work_array)

  vol=0.0
  DO n=1, nod2D
     vol=vol+area(n)
  END DO
  ocean_area=vol

  WRITE(*,*)  'Mesh statistics:'
  WRITE(*,*)  'maxArea: ', MAXVAL(elem_area), '   MinArea: ', MINVAL(elem_area)
  WRITE(*,*)  'maxScArea: ', MAXVAL(area(:)), '   MinScArea: ', MINVAL(area(:))
  WRITE(*,*) 'Total ocean area is: ', ocean_area, ' m^2'

  WRITE(*,*) 'mesh_areas finished'

  DEALLOCATE(nod_in_elem2D_num)


! *******************************************
! *** Read SST data from observation file ***
! *******************************************

  fillvalue = 1.0E7 

  ALLOCATE(average_obs(nod2D))
  ALLOCATE(average_std(nod2D))

  ! *** Initialize output file

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  WRITE (*,*) 'Write interpolated average SST observations to file: ',TRIM(ncfile_out)

  s = 1
  stat(s) = NF_CREATE(ncfile_out, 0, ncid_out)

  attstr  = 'interpolated SST observations'
  s = s + 1
  stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
       TRIM(attstr))
  s = s + 1
  stat(s) = NF_PUT_ATT_REAL(ncid_out, NF_GLOBAL, '_FillValue', NF_REAL, 1, REAL(fillvalue,4)) 

  ! Define dimensions
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'nodes_2D', nod2D, dimid_n2d)
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'time', day, dimid_time)

  ! Define variables
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'time', NF_REAL, 1, dimid_time, id_time)

  dimids(1) = dimid_n2D
  dimids(2) = dimid_time
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'sst', NF_REAL, 2, dimids, id_interp_SST)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'std', NF_REAL, 2, dimids, id_interp_std)
  s = s + 1
  stat(s) = NF_ENDDEF(ncid_out)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in init of output file, no.', i
  END DO

  ! Write time
  ALLOCATE (sec_in_year (day))
  DO i = 1, day
    sec_in_year (i) = 86400.0 * REAL(i)
  END DO

  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_time, REAL(sec_in_year, 4))

 
  ! *** Loop over one year ***
  cnt_day = 0

  monthloop: DO month = 1, 12
     dayloop: DO iday = 1, month_day(month)
        WRITE (obsfile,'(A,I2.2,I2.2,A)') TRIM(fileprefix), month, iday, TRIM(filesuffix)
        cnt_day = cnt_day + 1  

        ! Read netCDF file holding original observations
        ncfile_in = TRIM(obspath)//TRIM(obsfile)
        WRITE (*,'(8x,a,a)') 'Read observations from file: ', TRIM(ncfile_in)

        s = 1
        stat(s) = NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in)
        s = s + 1

        ! Get dimensions
        stat(s) = NF_INQ_DIMID(ncid_in, 'lon', id_dim)
        s = s + 1
        stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_x)
        s = s + 1
        stat(s) = NF_INQ_DIMID(ncid_in, 'lat', id_dim)
        s = s + 1
        stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_y)
  
        ! Get state variable IDs
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'adjusted_sea_surface_temperature', id_SST)
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'adjusted_standard_deviation_error', id_std)
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'lon', id_lon)
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'lat', id_lat) 
   
        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
        END DO

        ! *** Read observations and error information ***

        ! Allocate fields for file reading
        ALLOCATE(obs_file(dim_x, dim_y))
        ALLOCATE(robs_file(dim_x, dim_y))
        ALLOCATE(obs_lon(dim_x))
        ALLOCATE(obs_lat(dim_y))

        ! Read coordinate
        stat(1) = NF_GET_VAR_REAL(ncid_in, id_lon, obs_lon)
        stat(2) = NF_GET_VAR_REAL(ncid_in, id_lat, obs_lat)

        ! Read SST data

        startv(3) = 1
        countv(3) = 1
        startv(2) = 1
        countv(2) = dim_y 
        startv(1) = 1
        countv(1) = dim_x
        stat(1) = NF_GET_VARA_INT2(ncid_in, id_SST, startv(1:3), countv(1:3), obs_file)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading observation_SST'
        DO j = 1, dim_y
           DO i = 1, dim_x
              ! transfer to degree celcius
              robs_file(i,j) = REAL(obs_file(i,j))*0.01
           END DO
        END DO

        ! Read std data
        WRITE (*,*) '-- Read std data ---'
        ALLOCATE (std_file(dim_x, dim_y))
        ALLOCATE (rstd_file(dim_x, dim_y))     

        startv(3) = 1
        countv(3) = 1
        startv(2) = 1
        countv(2) = dim_y
        startv(1) = 1
        countv(1) = dim_x
        stat(1) = NF_GET_VARA_INT1(ncid_in, id_std, startv(1:3), countv(1:3), std_file)
        IF (stat(1) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading observation_SST'
        DO j = 1, dim_y
           DO i = 1, dim_x
              ! transform to degree celcius
              rstd_file(i,j) = REAL(std_file(i,j))*0.01 + 1
           END DO
        END DO
        WRITE (*,*) '--- Done ---'


!*************************************************************************
!*** Average original observations onto FESOM mesh and write into file ***
!*************************************************************************

        nobs_day = 0
        allcnt_day = 0
        maxcnt_day = 0
        mincnt_day = 10000

        DO n = 1, nod2D

           ! search radius in meters
           radius  = mesh_resolution(n) * avg_res_factor

           ! convert rotated coordinates to geographical ones

           ! Rotated Cartesian coordinates:
           xr=cos(coord_nod2D (2, n))*cos(coord_nod2D (1, n))
           yr=cos(coord_nod2D (2, n))*sin(coord_nod2D (1, n))
           zr=sin(coord_nod2D (2, n))

           ! Geographical Cartesian coordinates:
           xg=rotate_matrix(1,1)*xr + rotate_matrix(2,1)*yr + rotate_matrix(3,1)*zr
           yg=rotate_matrix(1,2)*xr + rotate_matrix(2,2)*yr + rotate_matrix(3,2)*zr
           zg=rotate_matrix(1,3)*xr + rotate_matrix(2,3)*yr + rotate_matrix(3,3)*zr

           ! Geographical coordinates:
           lat=asin(zg)
           if(yg==0. .and. xg==0.) then
              lon=0.0     ! exactly at the poles
           else
              lon=atan2(yg,xg)
           end if

           ! Convert distance to longitude and latitude  ---  approximately (in radians)
           radius_coord (1) = ABS (radius / (r_earth * COS(lat)))
           radius_coord (2) = radius / r_earth

           ! Define the search box and count the number of observations within the search box
           ! This is in radians
           xmin = lon - radius_coord (1)
           xmax = lon + radius_coord (1)
           ymin = lat - radius_coord (2)
           ymax = lat + radius_coord (2)

           ! Convert to degrees
           xmin = xmin / rad
           xmax = xmax / rad 
           ymin = ymin / rad
           ymax = ymax / rad 

           ! Count all observations within the search box and compute the average obs for each 2D node
           cnt  = 0
           mean_obs = 0
           mean_std = 0

           ! compute the start and end j
           jmin = INT((ymin + 79.95) * 10) + 1
           jmax = INT((ymax + 79.95) * 10) + 1
           imin = INT((xmin + 179.95) * 10) +1
           imax = INT((xmax + 179.95) * 10) + 1
           IF (jmin < 1) jmin = 1
           IF (jmax > dim_y) jmax = dim_y
           IF (imin < 1) imin = 1
           IF (imax > dim_x) imax = dim_x

           ! Averaging
           DO j_obs = jmin, jmax
              coord_y = 0.1 * REAL(j_obs) - 80.05
              IF (coord_y > ymin .AND. coord_y < ymax )  THEN
                 DO i_obs = imin, imax
                    coord_x = 0.1 * REAL(i_obs) - 180.05
                    IF (coord_x > xmin .AND. coord_x < xmax  &
                         .AND. robs_file (i_obs,j_obs) > obs_min ) THEN
                       cnt  = cnt + 1
                       mean_obs = mean_obs + robs_file(i_obs, j_obs)
                       mean_std = mean_std + rstd_file(i_obs, j_obs)
                    END IF
                 END DO
              END IF
           END DO

           IF (cnt /= 0) THEN
              mean_obs = mean_obs / cnt
              mean_std = mean_std / cnt
              nobs_day = nobs_day+1
           ELSE
              mean_obs = fillvalue
              mean_std = fillvalue
           END IF
           average_obs(n) = mean_obs
           average_std(n) = mean_std

           allcnt_day = allcnt_day + cnt

           if (cnt<mincnt_day) mincnt_day = cnt
           if (cnt>maxcnt_day) maxcnt_day = cnt
        END DO ! n=1, nod2D

        ! Write number of obs. per day
        write (*,*) 'Number of observations for this day: ', nobs_day
        write (*,'(a, f10.2, 2i7)') 'mean/min/max averaged number of raw observations: ', &
             REAL(allcnt_day) / REAL(nobs_day), mincnt_day, maxcnt_day

        !  Write SST

        startv_out(2) = cnt_day
        countv_out(2) = 1
        startv_out(1) = 1
        countv_out(1) = nod2D

        s = s + 1
        stat(s) = NF_PUT_VARA_REAL(ncid_out, id_interp_SST, startv_out, countv_out, REAL(average_obs, 4))

        !  Write std
        s = s + 1
        stat(s) = NF_PUT_VARA_REAL(ncid_out, id_interp_std, startv_out, countv_out, REAL(average_std, 4))
 
        DEALLOCATE(obs_file, robs_file, std_file, rstd_file,obs_lon,obs_lat)
     END DO dayloop
  END DO monthloop

  ! Close file
  s = 1
  stat(s) = NF_CLOSE(ncid_out)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in write of output file, no.', i
  END DO
  WRITE (*,*) '--- Done ---'

END PROGRAM obs_interp
