!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!mod_read_obs.F90: Module for reading of observation files
!-------------------------------------------------------------------------------------------

module mod_read_obs
  use iso_C_binding

  implicit none
  integer, allocatable :: idx_obs_nc(:), x_idx_obs_nc(:), y_idx_obs_nc(:), z_idx_obs_nc(:), var_id_obs_nc(:,:)
  integer(c_int), allocatable,target :: idx_obs_pf(:), x_idx_obs_pf(:), y_idx_obs_pf(:), z_idx_obs_pf(:), ind_obs_pf(:)
  type(c_ptr),bind(C,name="tidx_obs") :: ptr_tidx_obs
  type(c_ptr),bind(C,name="xidx_obs") :: ptr_xidx_obs
  type(c_ptr),bind(C,name="yidx_obs") :: ptr_yidx_obs
  type(c_ptr),bind(C,name="zidx_obs") :: ptr_zidx_obs
  type(c_ptr),bind(C,name="ind_obs")  :: ptr_ind_obs

  !kuw: obs variables for clm
  real, allocatable :: clmobs_lon(:), clmobs_lat(:)
  integer, allocatable :: clmobs_layer(:)
  !integer :: clmobs_layer
  real, allocatable :: clmobs_dr(:)
  real, allocatable :: clm_obs(:)
  real, allocatable :: clm_obserr(:)
  !kuw end

  real, allocatable :: pressure_obs(:)
  real, allocatable :: pressure_obserr(:)
  integer :: multierr=0, dim_nx, dim_ny
contains
  subroutine read_obs_nc()
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    use netcdf
    implicit none
    integer :: ncid, pres_varid, idx_varid,  x_idx_varid,  y_idx_varid,  z_idx_varid
    ! This is the name of the data file we will read.
    character (len = *), parameter :: dim_name = "dim_obs"
    character (len = *), parameter :: pres_name = "obs_pf"
    character (len = *), parameter :: idx_name = "idx"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status

    call check( nf90_open(obs_filename, nf90_nowrite, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    !print *, "dimid is ", dimid

    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    !print *, "name is ", recorddimname, ", len is ", dim_obs

    allocate(idx_obs_nc(dim_obs))
    allocate(pressure_obs(dim_obs))

    call check( nf90_inq_varid(ncid, pres_name, pres_varid) )
    call check( nf90_inq_varid(ncid, idx_name, idx_varid) )

    call check(nf90_get_var(ncid, pres_varid, pressure_obs))
    status =  nf90_get_var(ncid, idx_varid, idx_obs_nc)

    ! Read the surface pressure and idxerature data from the file.
    ! Since we know the contents of the file we know that the data
    ! arrays in this program are the correct size to hold all the data.

    allocate(x_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, X_IDX_NAME, x_idx_varid) )

    !print *, "pressure is ", pressure_obs
    !print *, "idx is ", idx_obs_nc
    call check( nf90_get_var(ncid, x_idx_varid, x_idx_obs_nc) )

    !print *, "x_idx is ", x_idx_obs_nc

    allocate(y_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Y_IDX_NAME, y_idx_varid) )

    call check( nf90_get_var(ncid, y_idx_varid, y_idx_obs_nc) )

    !print *, "y_idx is ", y_idx_obs_nc

    allocate(z_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Z_IDX_NAME, z_idx_varid) )

    call check( nf90_get_var(ncid, z_idx_varid, z_idx_obs_nc) )

    !print *, "z_idx is ", z_idx_obs_nc

    call check( nf90_close(ncid) )
    !print *,"*** SUCCESS reading example file ", obs_filename, "! "

  end subroutine read_obs_nc

  subroutine read_obs_nc_multi(current_observation_filename)
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    use netcdf
    implicit none
    integer :: ncid, pres_varid,presserr_varid, idx_varid,  x_idx_varid,  &
         y_idx_varid,  z_idx_varid
    ! This is the name of the data file we will read.
    character (len = *), parameter :: dim_name = "dim_obs"
    character (len = *), parameter :: pres_name = "obs_pf"
    character (len = *), parameter :: presserr_name = "obserr_pf"
    character (len = *), parameter :: idx_name = "idx"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status
    integer :: haserr
    character (len = *), intent(in) :: current_observation_filename

    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    !print *, "dimid is "!, dimid

    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    !print *, "name is ", recorddimname, ", len is ", dim_obs," dim_nx ", dim_nx
 
    allocate(idx_obs_nc(dim_obs))
    allocate(pressure_obs(dim_obs))

    call check( nf90_inq_varid(ncid, pres_name, pres_varid) )
    call check( nf90_inq_varid(ncid, idx_name, idx_varid) )

    call check(nf90_get_var(ncid, pres_varid, pressure_obs))
    status =  nf90_get_var(ncid, idx_varid, idx_obs_nc)

    ! Read the surface pressure and idxerature data from the file.
    ! Since we know the contents of the file we know that the data
    ! arrays in this program are the correct size to hold all the data.
    !check, if observation errors are present in observation file
    haserr = nf90_inq_varid(ncid, presserr_name, presserr_varid) 
    if(haserr == nf90_noerr) then
      multierr = 1
      if(.not.allocated(pressure_obserr)) allocate(pressure_obserr(dim_obs))
      call check(nf90_get_var(ncid, presserr_varid, pressure_obserr))
    endif

    if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
    allocate(x_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, X_IDX_NAME, x_idx_varid) )

    !print *, "pressure is ", pressure_obs
    !print *, "idx is ", idx_obs_nc
    call check( nf90_get_var(ncid, x_idx_varid, x_idx_obs_nc) )

    !print *, "x_idx is ", x_idx_obs_nc

    if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
    allocate(y_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Y_IDX_NAME, y_idx_varid) )

    call check( nf90_get_var(ncid, y_idx_varid, y_idx_obs_nc) )

    !print *, "y_idx is ", y_idx_obs_nc

    if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
    allocate(z_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Z_IDX_NAME, z_idx_varid) )

    call check( nf90_get_var(ncid, z_idx_varid, z_idx_obs_nc) )

    !print *, "z_idx is ", z_idx_obs_nc

    call check( nf90_close(ncid) )
    !print *,"*** SUCCESS reading example file ", current_observation_filename, "! "

  end subroutine read_obs_nc_multi

  subroutine read_obs_nc_multiscalar()
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    USE mod_parallel_pdaf, &
         ONLY: comm_filter
    USE mod_parallel_model, &
         ONLY: mpi_info_null
    USE netcdf
    implicit none
    integer :: ncid, pres_varid, idx_varid,  x_idx_varid,  y_idx_varid,  z_idx_varid, &
               var_id_varid
    integer :: comm, omode, info
    ! This is the name of the data file we will read.
    character (len = *), parameter :: dim_name = "dim_obs"
    character (len = *), parameter :: pres_name = "obs_pf"
    character (len = *), parameter :: idx_name = "idx"
    character (len = *), parameter :: var_id_name = "var_id"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character (len = *), parameter :: dim_nx_name = "dim_nx"
    character (len = *), parameter :: dim_ny_name = "dim_ny"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status

     call check( nf90_open(obs_filename, nf90_nowrite, ncid) )
    !call check( nf90_open_par(obs_filename, nf90_nowrite, comm_filter, mpi_info_null, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    !print *, "dimid is ", dimid

    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    !print *, "name is ", recorddimname, ", len is ", dim_obs
    call check(nf90_inq_dimid(ncid, dim_nx_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_nx))
    call check(nf90_inq_dimid(ncid, dim_ny_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_ny))

    allocate(idx_obs_nc(dim_obs))
    allocate(pressure_obs(dim_obs))
    allocate(var_id_obs_nc(dim_ny, dim_nx))

    call check( nf90_inq_varid(ncid, pres_name, pres_varid) )
    call check( nf90_inq_varid(ncid, idx_name, idx_varid) )

    call check(nf90_get_var(ncid, pres_varid, pressure_obs))
    status =  nf90_get_var(ncid, idx_varid, idx_obs_nc)
    call check(nf90_inq_varid(ncid, var_id_name, var_id_varid))
    call check(nf90_get_var(ncid, var_id_varid, var_id_obs_nc))

    ! Read the surface pressure and idxerature data from the file.
    ! Since we know the contents of the file we know that the data
    ! arrays in this program are the correct size to hold all the data.

    allocate(x_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, X_IDX_NAME, x_idx_varid) )

    !print *, "pressure is ", pressure_obs
    !print *, "idx is ", idx_obs_nc
    call check( nf90_get_var(ncid, x_idx_varid, x_idx_obs_nc) )

    !print *, "x_idx is ", x_idx_obs_nc

    allocate(y_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Y_IDX_NAME, y_idx_varid) )

    call check( nf90_get_var(ncid, y_idx_varid, y_idx_obs_nc) )

    !print *, "y_idx is ", y_idx_obs_nc

    allocate(z_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Z_IDX_NAME, z_idx_varid) )

    call check( nf90_get_var(ncid, z_idx_varid, z_idx_obs_nc) )

    !print *, "z_idx is ", z_idx_obs_nc

    call check( nf90_close(ncid) )
    !print *,"*** SUCCESS reading example file ", obs_filename, "! "

  end subroutine read_obs_nc_multiscalar 

subroutine read_obs_nc_multiscalar_files(current_observation_filename)
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    use netcdf
    implicit none
    integer :: ncid, pres_varid,presserr_varid, idx_varid,  x_idx_varid,  y_idx_varid,  &
               z_idx_varid, var_id_varid
    ! This is the name of the data file we will read.
    character (len = *), parameter :: dim_name = "dim_obs"
    character (len = *), parameter :: pres_name = "obs_pf"
    character (len = *), parameter :: presserr_name = "obserr_pf"
    character (len = *), parameter :: idx_name = "idx"
    character (len = *), parameter :: var_id_name = "var_id"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character (len = *), parameter :: dim_nx_name = "dim_nx"
    character (len = *), parameter :: dim_ny_name = "dim_ny"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status
    integer :: haserr
    character (len = *), intent(in) :: current_observation_filename

    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    !print *, "dimid is ", dimid

    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    call check(nf90_inq_dimid(ncid, dim_nx_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_nx))
    call check(nf90_inq_dimid(ncid, dim_ny_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_ny))
    !print *, "name is ", recorddimname, ", len is ", dim_obs
    !print *, "dim_nx dim_nx ", dim_nx, dim_ny

    allocate(idx_obs_nc(dim_obs))
    allocate(pressure_obs(dim_obs))
    allocate(var_id_obs_nc(dim_ny, dim_nx))

    call check( nf90_inq_varid(ncid, pres_name, pres_varid) )
    call check( nf90_inq_varid(ncid, idx_name, idx_varid) )
    call check(nf90_get_var(ncid, pres_varid, pressure_obs))
    status =  nf90_get_var(ncid, idx_varid, idx_obs_nc)
    call check(nf90_inq_varid(ncid, var_id_name, var_id_varid))
    call check(nf90_get_var(ncid, var_id_varid, var_id_obs_nc))

    ! Read the surface pressure and idxerature data from the file.
    ! Since we know the contents of the file we know that the data
    ! arrays in this program are the correct size to hold all the data.
    !check, if observation errors are present in observation file
    haserr = nf90_inq_varid(ncid, presserr_name, presserr_varid) 
    if(haserr == nf90_noerr) then
      multierr = 1
      if(.not.allocated(pressure_obserr)) allocate(pressure_obserr(dim_obs))
      call check(nf90_get_var(ncid, presserr_varid, pressure_obserr))
    endif

    if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
    allocate(x_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, X_IDX_NAME, x_idx_varid) )

    !print *, "pressure is ", pressure_obs
    !print *, "idx is ", idx_obs_nc
    call check( nf90_get_var(ncid, x_idx_varid, x_idx_obs_nc) )

    !print *, "x_idx is ", x_idx_obs_nc

    if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
    allocate(y_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Y_IDX_NAME, y_idx_varid) )

    call check( nf90_get_var(ncid, y_idx_varid, y_idx_obs_nc) )

    !print *, "y_idx is ", y_idx_obs_nc


    if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
    allocate(z_idx_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, Z_IDX_NAME, z_idx_varid) )

    call check( nf90_get_var(ncid, z_idx_varid, z_idx_obs_nc) )

    !print *, "z_idx is ", z_idx_obs_nc

    call check( nf90_close(ncid) )
    !print *,"*** SUCCESS reading example file ", current_observation_filename, "! "

  end subroutine read_obs_nc_multiscalar_files

  subroutine get_obsindex_currentobsfile(no_obs) bind(c,name='get_obsindex_currentobsfile')
    USE mod_assimilation, only: obs_filename
    use netcdf
    use mod_parallel_model, only:tcycle

    implicit none
    integer, intent(out) :: no_obs
    character (len = 110) :: filename
    integer :: ncid, varid
    character (len = *), parameter :: dim_name   = "dim_obs"
    character (len = *), parameter :: idx_name   = "idx"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character (len = *), parameter :: gwind_name = "gw_indicator"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status
    integer :: haserr

    !write(filename, '(a, i5.5)') trim(obs_filename)//'.', tcycle

    if(allocated(idx_obs_pf))   deallocate(idx_obs_pf)
    if(allocated(x_idx_obs_pf)) deallocate(x_idx_obs_pf)
    if(allocated(y_idx_obs_pf)) deallocate(y_idx_obs_pf)
    if(allocated(z_idx_obs_pf)) deallocate(z_idx_obs_pf)
    if(allocated(ind_obs_pf))   deallocate(ind_obs_pf)

    call check( nf90_open(filename, nf90_nowrite, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, no_obs))

    if(.not.allocated(idx_obs_pf)) allocate(idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, idx_name, varid) )
    status =  nf90_get_var(ncid, varid, idx_obs_pf)

    if(.not.allocated(x_idx_obs_pf)) allocate(x_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, X_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, x_idx_obs_pf) )

    if(.not.allocated(y_idx_obs_pf)) allocate(y_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, Y_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, y_idx_obs_pf) )

    if(.not.allocated(z_idx_obs_pf)) allocate(z_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, Z_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, z_idx_obs_pf) )

    if(.not.allocated(ind_obs_pf)) allocate(ind_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, gwind_name, varid) )
    call check( nf90_get_var(ncid, varid, ind_obs_pf) )

    ptr_tidx_obs = c_loc(idx_obs_pf)
    ptr_xidx_obs = c_loc(x_idx_obs_pf)
    ptr_yidx_obs = c_loc(y_idx_obs_pf)
    ptr_zidx_obs = c_loc(z_idx_obs_pf)
    ptr_ind_obs  = c_loc(ind_obs_pf)

    call check( nf90_close(ncid) )

  end subroutine get_obsindex_currentobsfile

  !kuw: routine to read clm soil moisture observations
  subroutine read_obs_nc_multi_clm(current_observation_filename)
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    use netcdf
    implicit none
    integer :: ncid, clmobs_varid, dr_varid,  clmobs_lon_varid,  clmobs_lat_varid,  &
         clmobs_layer_varid, clmobserr_varid
    character (len = *), parameter :: dim_name   = "dim_obs"
    character (len = *), parameter :: obs_name   = "obs_clm"
    character (len = *), parameter :: dr_name    = "dr"
    character (len = *), parameter :: lon_name   = "lon"
    character (len = *), parameter :: lat_name   = "lat"
    character (len = *), parameter :: layer_name = "layer"
    character (len = *), parameter :: obserr_name   = "obserr_clm"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status, haserr
    character (len = *), intent(in) :: current_observation_filename

    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )
    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))

    allocate(clmobs_lon(dim_obs))
    allocate(clmobs_lat(dim_obs))
    allocate(clm_obs(dim_obs))
    allocate(clmobs_layer(dim_obs))
    allocate(clmobs_dr(2))

    call check( nf90_inq_varid(ncid, obs_name, clmobs_varid) )
    call check(nf90_get_var(ncid, clmobs_varid, clm_obs))

    !call check( nf90_inq_varid(ncid, dr_name, dr_varid) )

    !check, if observation errors are present in observation file
    haserr = nf90_inq_varid(ncid, obserr_name, clmobserr_varid) 
    if(haserr == nf90_noerr) then
      multierr = 1
      if(.not.allocated(clm_obserr)) allocate(clm_obserr(dim_obs))
      call check(nf90_get_var(ncid, clmobserr_varid, clm_obserr))
    endif


    call check( nf90_inq_varid(ncid, lon_name, clmobs_lon_varid) )
    call check( nf90_get_var(ncid, clmobs_lon_varid, clmobs_lon) )

    call check( nf90_inq_varid(ncid, lat_name, clmobs_lat_varid) )
    call check( nf90_get_var(ncid, clmobs_lat_varid, clmobs_lat) )

    call check( nf90_inq_varid(ncid, layer_name, clmobs_layer_varid) )
    call check( nf90_get_var(ncid, clmobs_layer_varid, clmobs_layer) )

    call check( nf90_inq_varid(ncid, dr_name, dr_varid) )
    call check( nf90_get_var(ncid, dr_varid, clmobs_dr) )

    call check( nf90_close(ncid) )

  end subroutine read_obs_nc_multi_clm
    !kuw end
  subroutine read_obs_nc_multiscalar_clm_files(current_observation_filename)
    USE mod_assimilation, &
         ONLY: obs_p, obs_index_p, dim_obs, obs_filename
    use netcdf
    implicit none
    integer :: ncid, clmobs_varid, dr_varid,  clmobs_lon_varid,  clmobs_lat_varid,  &
               clmobs_layer_varid, clmobserr_varid, var_id_varid, x, y
    character (len = *), parameter :: dim_name   = "dim_obs"
    character (len = *), parameter :: dim_nx_name = "dim_nx"
    character (len = *), parameter :: dim_ny_name = "dim_ny"
    character (len = *), parameter :: obs_name   = "obs_clm"
    character (len = *), parameter :: var_id_name = "var_id"
    character (len = *), parameter :: dr_name    = "dr"
    character (len = *), parameter :: lon_name   = "lon"
    character (len = *), parameter :: lat_name   = "lat"
    character (len = *), parameter :: layer_name = "layer"
    character (len = *), parameter :: obserr_name   = "obserr_clm"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status, haserr
    character (len = *), intent(in) :: current_observation_filename

    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )
    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    call check(nf90_inq_dimid(ncid, dim_nx_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_nx))
    call check(nf90_inq_dimid(ncid, dim_ny_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_ny))

    allocate(clmobs_lon(dim_obs))
    allocate(clmobs_lat(dim_obs))
    allocate(clm_obs(dim_obs)) 
    allocate(clmobs_layer(dim_obs))
    allocate(clmobs_dr(2))
    allocate(var_id_obs_nc(dim_ny, dim_nx))
    !allocate(var_id_obs_nc(dim_obs))

    call check( nf90_inq_varid(ncid, obs_name, clmobs_varid) )
    call check(nf90_get_var(ncid, clmobs_varid, clm_obs))
    call check(nf90_inq_varid(ncid, var_id_name, var_id_varid))
    call check(nf90_get_var(ncid, var_id_varid, var_id_obs_nc))
    call check( nf90_inq_varid(ncid, dr_name, dr_varid) )

    !check, if observation errors are present in observation file
    haserr = nf90_inq_varid(ncid, obserr_name, clmobserr_varid) 
    if(haserr == nf90_noerr) then
      multierr = 1
      if(.not.allocated(clm_obserr)) allocate(clm_obserr(dim_obs))
      call check(nf90_get_var(ncid, clmobserr_varid, clm_obserr))
    endif

    call check( nf90_inq_varid(ncid, lon_name, clmobs_lon_varid) )
    call check( nf90_get_var(ncid, clmobs_lon_varid, clmobs_lon) )

    call check( nf90_inq_varid(ncid, lat_name, clmobs_lat_varid) )
    call check( nf90_get_var(ncid, clmobs_lat_varid, clmobs_lat) )

    call check( nf90_inq_varid(ncid, layer_name, clmobs_layer_varid) )
    call check( nf90_get_var(ncid, clmobs_layer_varid, clmobs_layer) )

    call check( nf90_inq_varid(ncid, dr_name, dr_varid) )
    call check( nf90_get_var(ncid, dr_varid, clmobs_dr) )

    call check( nf90_close(ncid) )

  end subroutine read_obs_nc_multiscalar_clm_files

  subroutine clean_obs_nc()
    implicit none
   ! if(allocated(idx_obs_nc))deallocate(idx_obs_nc)
    if(allocated(pressure_obs))deallocate(pressure_obs)
    !if(allocated(pressure_obserr))deallocate(pressure_obserr)
    !if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
    !if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
    !if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
    !kuw: clean clm observations
    if(allocated(clmobs_lon))deallocate(clmobs_lon)
    if(allocated(clmobs_lat))deallocate(clmobs_lat)
    if(allocated(clm_obs))deallocate(clm_obs)
    if(allocated(clmobs_layer))deallocate(clmobs_layer)
    if(allocated(clmobs_dr))deallocate(clmobs_dr)
    !if(allocated(clm_obserr))deallocate(clm_obserr)
    !kuw end
  end subroutine clean_obs_nc

  subroutine clean_obs_pf() bind(c,name='clean_obs_pf')
    implicit none
    if(allocated(idx_obs_pf))deallocate(idx_obs_pf)
    if(allocated(x_idx_obs_pf))deallocate(x_idx_obs_pf)
    if(allocated(y_idx_obs_pf))deallocate(y_idx_obs_pf)
    if(allocated(z_idx_obs_pf))deallocate(z_idx_obs_pf)
  end subroutine clean_obs_pf

  subroutine check(status)

    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check


end module mod_read_obs
