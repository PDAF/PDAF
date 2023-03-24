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
!init_dim_obs_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                       'init_dim_obs_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_obs_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
  !
  ! The routine is called at the beginning of each
  ! analysis step.  It has to initialize the size of
  ! the observation vector according to the current
  ! time step for the PE-local domain.
  !
  ! !REVISION HISTORY:
  ! 2013-02 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  !   USE mod_assimilation, &
  !        ONLY : nx, ny, local_dims, obs_p, obs_index_p, &
  !        coords_obs, local_dims_obs
  !   USE mod_parallel_pdaf, &
  !        ONLY: mype_filter, npes_filter, COMM_filter, MPI_INTEGER, &
  !        MPIerr, MPIstatus
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, comm_filter, npes_filter
  use mod_parallel_model, &
       only: mpi_integer, model, mpi_double_precision, mpi_in_place, mpi_sum, &
       mype_world
  USE mod_assimilation, &
       ONLY: obs_p, obs_index_p, dim_obs, obs_filename, &
       obs, &
       pressure_obserr_p, clm_obserr_p, &
       obs_nc2pdaf, &
       local_dims_obs,
       ! dim_obs_p, &
       dim_obs_f, &
       obs_id_p, &
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
!hcp 
!CLMSA needs the physical  coordinates of the elements of state vector 
!and observation array.        
       longxy, latixy, longxy_obs, latixy_obs, &
!hcp end
#endif
#endif
       var_id_obs, maxlon, minlon, maxlat, &
       minlat, maxix, minix, maxiy, miniy, lon_var_id, ix_var_id, lat_var_id, iy_var_id, &
       screen
  Use mod_read_obs, &
       only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, &
       read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, &
       z_idx_obs_nc, clm_obs, &
       var_id_obs_nc, dim_nx, dim_ny, &
       clmobs_lon, clmobs_lat, clmobs_layer, clmobs_dr, clm_obserr
  use mod_tsmp, &
      only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
      xcoord, ycoord, zcoord, xcoord_fortran, ycoord_fortran, &
      zcoord_fortran, &
#endif
#endif
      tag_model_clm, point_obs

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  !kuw
  use shr_kind_mod, only: r8 => shr_kind_r8
  USE clmtype,                  ONLY : clm3
  use decompMod , only : get_proc_bounds, get_proc_global
  !kuw end
  !hcp
  !use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
  !latixy, longxy_obs, latixy_obs
  USE enkf_clm_mod, only: domain_def_clm
  !hcp end
#endif
#endif

  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  ! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector
  ! !CALLING SEQUENCE:
  ! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
  ! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
  ! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
  ! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
  ! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
  !EOP

  ! *** Local variables
  integer :: ierror
  INTEGER :: max_var_id
  INTEGER :: tmp_dim_obs_f
  INTEGER :: i,j,k,count  ! Counters
  INTEGER :: m,l          ! Counters
  logical :: is_multi_observation_files
  character (len = 110) :: current_observation_filename
  integer,allocatable :: local_dis(:),local_dim(:)

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)
  ! pft: "plant functional type"
  integer :: begp, endp   ! per-proc beginning and ending pft indices
  integer :: begc, endc   ! per-proc beginning and ending column indices
  integer :: begl, endl   ! per-proc beginning and ending landunit indices
  integer :: begg, endg   ! per-proc gridcell ending gridcell indices
  integer :: numg         ! total number of gridcells across all processors
  integer :: numl         ! total number of landunits across all processors
  integer :: numc         ! total number of columns across all processors
  integer :: nump         ! total number of pfts across all processors
  real    :: deltax, deltay
  !real    :: deltaxy, y1 , x1, z1, x2, y2, z2, R, dist, deltaxy_max
  logical :: is_use_dr
#endif
#endif

  ! ****************************************
  ! *** Initialize observation dimension ***
  ! ****************************************

  ! Read observation file
  ! ---------------------

  !  if I'm root in filter, read the nc file
  is_multi_observation_files = .true.
  if (is_multi_observation_files) then
      ! Set name of current NetCDF observation file
      write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step
  else
      ! Single NetCDF observation file (currently NOT used)
      write(current_observation_filename, '(a, i5.5)') trim(obs_filename)
  end if

  if (mype_filter .eq. 0) then
      ! Read current NetCDF observation file
      call read_obs_nc(current_observation_filename)
  end if

  ! Broadcast first variables
  ! -------------------------
  ! Dimension of observation vector
  call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! Switch for vector of observation errors
  call mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! broadcast dim_ny and dim_nx
  if(point_obs.eq.0) then
     call mpi_bcast(dim_nx, 1, MPI_INTEGER, 0, comm_filter, ierror)
     call mpi_bcast(dim_ny, 1, MPI_INTEGER, 0, comm_filter, ierror)
  endif

  ! Allocate observation arrays for non-root procs
  ! ----------------------------------------------
  if (mype_filter .ne. 0) then ! for all non-master proc
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
      ! if exist ParFlow-type obs
     !if(model == tag_model_parflow) then
        if(allocated(pressure_obs)) deallocate(pressure_obs)
        allocate(pressure_obs(dim_obs))
        if (multierr.eq.1) then
             if (allocated(pressure_obserr)) deallocate(pressure_obserr)
             allocate(pressure_obserr(dim_obs))
        endif
        if(allocated(idx_obs_nc)) deallocate(idx_obs_nc)
        allocate(idx_obs_nc(dim_obs))
        if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
        allocate(x_idx_obs_nc(dim_obs))
        if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
        allocate(y_idx_obs_nc(dim_obs))
        if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
        allocate(z_idx_obs_nc(dim_obs))
        if(point_obs.eq.0) then
           if(allocated(var_id_obs_nc))deallocate(var_id_obs_nc)
           allocate(var_id_obs_nc(dim_ny, dim_nx))
        endif
     !end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
      ! if exist CLM-type obs
!     if(model == tag_model_clm) then
        if(allocated(clm_obs)) deallocate(clm_obs)
        allocate(clm_obs(dim_obs))
        if(allocated(clmobs_lon)) deallocate(clmobs_lon)
        allocate(clmobs_lon(dim_obs))
        if(allocated(clmobs_lat)) deallocate(clmobs_lat)
        allocate(clmobs_lat(dim_obs))
        if(allocated(clmobs_dr)) deallocate(clmobs_dr)
        allocate(clmobs_dr(2))
        if(allocated(clmobs_layer)) deallocate(clmobs_layer)
        allocate(clmobs_layer(dim_obs))
        if(point_obs.eq.0) then
            if(allocated(var_id_obs_nc)) deallocate(var_id_obs_nc)
            allocate(var_id_obs_nc(dim_ny, dim_nx))
        endif
        if(multierr.eq.1) then 
            if(allocated(clm_obserr)) deallocate(clm_obserr)
            allocate(clm_obserr(dim_obs))
        end if
!     end if
#endif
#endif
  end if

  ! Broadcast the idx and pressure
  ! ------------------------------
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  !if(model == tag_model_parflow) then
      ! if exist ParFlow-type obs
     call mpi_bcast(pressure_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     if(multierr.eq.1) call mpi_bcast(pressure_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     ! broadcast xyz indices
     call mpi_bcast(x_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     call mpi_bcast(y_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     call mpi_bcast(z_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     if(point_obs.eq.0) call mpi_bcast(var_id_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  !end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  !if(model == tag_model_clm) then
      ! if exist CLM-type obs
     call mpi_bcast(clm_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     if(multierr.eq.1) call mpi_bcast(clm_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_lon, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_lat, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_dr,  2, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_layer, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     if(point_obs.eq.0) call mpi_bcast(var_id_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  !end if
#endif
#endif

  ! CLM grid information
  ! --------------------
  ! Results used only in `localize_covar_pdaf` for LEnKF
  ! Calling could be restricted to LEnKF
!hcp
!use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
!latixy, longxy_obs, latixy_obs
! Index arrays of longitudes and latitudes
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
     ! if exist CLM-type obs
  if(model .eq. tag_model_clm) then
      ! Generate CLM index arrays from lon/lat values
      call domain_def_clm(clmobs_lon, clmobs_lat, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

      ! Obtain general CLM index information
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
      call get_proc_global(numg, numl, numc, nump)
  end if
#endif
#endif
!hcp end

  ! Number of observations in process-local domain
  ! ----------------------------------------------
  ! Additionally `obs_id_p` is set (the NetCDF index of the
  ! observation corresponding to the state index in the local domain)
  dim_obs_p = 0

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  if (model .eq. tag_model_parflow) then

     if(allocated(obs_id_p)) deallocate(obs_id_p)
     allocate(obs_id_p(enkf_subvecsize))
     obs_id_p(:) = 0

     do i = 1, dim_obs
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              dim_obs_p = dim_obs_p + 1
              obs_id_p(j) = i
           end if
        end do
     end do
  end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  ! Switch for how to check index of CLM observations
  ! True: Use snapping distance between long/lat on CLM grid
  ! False: Use index arrays from `domain_def_clm`
  is_use_dr = .true.
  
  if(model .eq. tag_model_clm) then

     if(allocated(obs_id_p)) deallocate(obs_id_p)
     allocate(obs_id_p(endg-begg+1))
     obs_id_p(:) = 0

     do i = 1, dim_obs
        count = 1
        do j = begg, endg
            if(is_use_dr) then
                deltax = abs(lon(j)-clmobs_lon(i))
                deltay = abs(lat(j)-clmobs_lat(i))
            end if
            if(((is_use_dr).and.(deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2)))
                .or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(count)) .and. (latixy_obs(i) == latixy(count)))) then
                dim_obs_p = dim_obs_p + 1
                obs_id_p(count) = i
                ! One observation can be associated with multiple grid
                ! points If it is clear that each observation is only
                ! associated to a single grid point, then the
                ! following line can be uncommented
                ! exit
            end if
            count = count + 1
        end do
    end do
  end if
#endif
#endif

  if (screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dim_obs_p=", dim_obs_p
  end if

  ! add and broadcast size of local observation dimensions using mpi_allreduce 
  call mpi_allreduce(dim_obs_p, tmp_dim_obs_f, 1, MPI_INTEGER, MPI_SUM, &
       comm_filter, ierror) 
  ! Set dimension of full observation vector
  dim_obs_f = tmp_dim_obs_f

  if (screen > 2) then
      if (mype_filter==0) then
          print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_f_pdaf: dim_obs_f=", dim_obs_f
      end if
  end if

  allocate(local_dis(npes_filter))
  allocate(local_dim(npes_filter))
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dim, 1, MPI_INTEGER, comm_filter, ierror)
  local_dis(1) = 0
  do i = 2, npes_filter
     local_dis(i) = local_dis(i-1) + local_dim(i-1)
  end do
  deallocate(local_dim)

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: local_dis=", local_dis
  end if

  ! Write process-local observation arrays
  ! --------------------------------------
  ! allocate index for mapping between observations in nc input and
  ! sorted by pdaf: obs_nc2pdaf

  ! Non-trivial example: The second observation in the NetCDF file
  ! (`i=2`) is the only observation in the subgrid (`count = 1`) of
  ! the first PE (`mype_filter = 0`):
  !
  ! i = 2
  ! count = 1
  ! mype_filter = 0
  ! 
  ! obs_nc2pdaf(local_dis(mype_filter+1)+count) = i
  !-> obs_nc2pdaf(local_dis(1)+1) = 2
  !-> obs_nc2pdaf(1) = 2

  IF (ALLOCATED(obs)) DEALLOCATE(obs)
  ALLOCATE(obs(dim_obs))
  !IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
  !ALLOCATE(obs_index(dim_obs))
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_p(dim_obs_p))
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  if(point_obs.eq.0) then
      IF (ALLOCATED(var_id_obs)) DEALLOCATE(var_id_obs)
      ALLOCATE(var_id_obs(dim_obs_p))
  end if

  if (allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
  allocate(obs_nc2pdaf(dim_obs))
  obs_nc2pdaf = 0

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  if (model .eq. tag_model_parflow) then
     ! allocate pressure_obserr_p observation error for parflow run at PE-local domain 
!     if((multierr.eq.1) .and. (.not.allocated(pressure_obserr_p))) allocate(pressure_obserr_p(dim_obs_p))
     !hcp pressure_obserr_p must be reallocated because the numbers of obs are
     !not necessary the same for all observation files.
     if(multierr.eq.1) then 
        if (allocated(pressure_obserr_p)) deallocate(pressure_obserr_p)
        allocate(pressure_obserr_p(dim_obs_p))
     endif
     !hcp fin

  if (point_obs.eq.0) then
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     if(allocated(ix_var_id)) deallocate(ix_var_id) 
     allocate(ix_var_id(max_var_id))
     if(allocated(iy_var_id)) deallocate(iy_var_id)
     allocate(iy_var_id(max_var_id))
     if(allocated(maxix)) deallocate(maxix)
     allocate(maxix(max_var_id))
     if(allocated(minix)) deallocate(minix)
     allocate(minix(max_var_id))
     if(allocated(maxiy)) deallocate(maxiy)
     allocate(maxiy(max_var_id))
     if(allocated(miniy)) deallocate(miniy)
     allocate(miniy(max_var_id))
     
     ix_var_id(:) = 0
     iy_var_id(:) = 0
     maxix = -999
     minix = 9999999
     maxiy = -999
     miniy = 9999999
     do j = 1, max_var_id
        do m = 1, dim_nx
           do k = 1, dim_ny 
              i = (m-1)* dim_ny + k  
              if (var_id_obs_nc(k,m) == j) then      
                 maxix(j) = MAX(x_idx_obs_nc(i),maxix(j))
                 minix(j) = MIN(x_idx_obs_nc(i),minix(j))
                 maxiy(j) = MAX(y_idx_obs_nc(i),maxiy(j))
                 miniy(j) = MIN(y_idx_obs_nc(i),miniy(j))
              end if
           end do
        end do
        ix_var_id(j) = (maxix(j) + minix(j))/2.0
        iy_var_id(j) = (maxiy(j) + miniy(j))/2.0
     end do

     count = 1
     do m = 1, dim_nx
        do k = 1, dim_ny
           i = (m-1)* dim_ny + k    
           obs(i) = pressure_obs(i)  
           ! coords_obs(1, i) = idx_obs_nc(i)
           do j = 1, enkf_subvecsize
              if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
                 obs_index_p(count) = j
                 obs_p(count) = pressure_obs(i)
                 var_id_obs(count) = var_id_obs_nc(k,m)
                 if(multierr.eq.1) pressure_obserr_p(count) = pressure_obserr(i)
                 count = count + 1
              end if
           end do
        end do
     end do
  else if (point_obs.eq.1) then
     count = 1
     do i = 1, dim_obs
        obs(i) = pressure_obs(i)  
        ! coords_obs(1, i) = idx_obs_nc(i)
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              !print *, j
              !obs_index(count) = j
              !obs(count) = pressure_obs(i)
              obs_index_p(count) = j
              obs_p(count) = pressure_obs(i)
              if(multierr.eq.1) pressure_obserr_p(count) = pressure_obserr(i)
              obs_nc2pdaf(local_dis(mype_filter+1)+count) = i
              count = count + 1
           end if
        end do
     end do
  end if
  end if
  call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  if(model .eq. tag_model_clm) then
     ! allocate clm_obserr_p observation error for clm run at PE-local domain
!     if((multierr.eq.1) .and. (.not.allocated(clm_obserr_p))) allocate(clm_obserr_p(dim_obs_p))
     if(multierr.eq.1) then
         if (allocated(clm_obserr_p)) deallocate(clm_obserr_p)
         allocate(clm_obserr_p(dim_obs_p))
     endif
  if(point_obs.eq.0) then
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     if(allocated(lon_var_id)) deallocate(lon_var_id)
     allocate(lon_var_id(max_var_id))
     if(allocated(lat_var_id)) deallocate(lat_var_id)
     allocate(lat_var_id(max_var_id))
     if(allocated(maxlon)) deallocate(maxlon)
     allocate(maxlon(max_var_id))
     if(allocated(minlon)) deallocate(minlon)
     allocate(minlon(max_var_id))
     if(allocated(maxlat)) deallocate(maxlat)
     allocate(maxlat(max_var_id))
     if(allocated(minlat)) deallocate(minlat)
     allocate(minlat(max_var_id))

     lon_var_id(:) = 0
     lat_var_id(:) = 0
     maxlon = -999
     minlon = 9999999
     maxlat = -999
     minlat = 9999999
     do j = 1, max_var_id
        do m = 1, dim_nx
           do k = 1, dim_ny
              i = (m-1)* dim_ny + k    
              if (var_id_obs_nc(k,m) == j) then      
                 maxlon(j) = MAX(longxy_obs(i),maxlon(j))
                 minlon(j) = MIN(longxy_obs(i),minlon(j))
                 maxlat(j) = MAX(latixy_obs(i),maxlat(j))
                 minlat(j) = MIN(latixy_obs(i),minlat(j))
              end if
           end do
           lon_var_id(j) = (maxlon(j) + minlon(j))/2.0 
           lat_var_id(j) = (maxlat(j) + minlat(j))/2.0
           !print *, 'j  lon_var_id  lat_var_id ', j, lon_var_id(j), lat_var_id(j)
        enddo  ! allocate clm_obserr_p observation error for clm run at PE-local domain
     enddo

     count = 1
     do m = 1, dim_nx
        do l = 1, dim_ny
           i = (m-1)* dim_ny + l        
           obs(i) = clm_obs(i) 
           k = 1
           do j = begg,endg
              if((longxy_obs(i) == longxy(k)) .and. (latixy_obs(i) == latixy(k))) then
                 obs_index_p(count) = k 
                 obs_p(count) = clm_obs(i)
                 var_id_obs(count) = var_id_obs_nc(l,m)
                 if(multierr.eq.1) clm_obserr_p(count) = clm_obserr(i)
                 count = count + 1
              endif
              k = k + 1
           end do
        end do
     end do
  else if(point_obs.eq.1) then

     count = 1
     do i = 1, dim_obs
        obs(i) = clm_obs(i) 
        k = 1
       do j = begg,endg
            if(is_use_dr) then
                deltax = abs(lon(j)-clmobs_lon(i))
                deltay = abs(lat(j)-clmobs_lat(i))
            end if
            if(((is_use_dr).and.(deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2)))
                .or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(k)) .and. (latixy_obs(i) == latixy(k)))) then
              !obs_index_p(count) = j + (size(lon) * (clmobs_layer(i)-1))
              !obs_index_p(count) = j + ((endg-begg+1) * (clmobs_layer(i)-1))
              !obs_index_p(count) = j-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
              obs_index_p(count) = k + ((endg-begg+1) * (clmobs_layer(i)-1))
              !write(*,*) 'obs_index_p(',count,') is',obs_index_p(count)
              obs_p(count) = clm_obs(i)
              if(multierr.eq.1) clm_obserr_p(count) = clm_obserr(i)
              obs_nc2pdaf(local_dis(mype_filter+1)+count) = i
              count = count + 1
           end if
           k = k + 1
        end do
     end do
  end if
  end if
  call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
#endif
#endif

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: obs_nc2pdaf=", obs_nc2pdaf
  end if

  ! allocate array of local observation dimensions with total PEs
  IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
  ALLOCATE(local_dims_obs(npes_filter))

  ! Gather array of local observation dimensions 
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
       comm_filter, ierror)

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
!!#if (defined PARFLOW_STAND_ALONE || defined COUP_OAS_PFL)
  IF (model == tag_model_parflow) THEN
     !print *, "Parflow: converting xcoord to fortran"
     call C_F_POINTER(xcoord, xcoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(ycoord, ycoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(zcoord, zcoord_fortran, [enkf_subvecsize])
  ENDIF
#endif
#endif

  !  clean up the temp data from nc file
  ! ------------------------------------
  deallocate(local_dis)
  call clean_obs_nc()

END SUBROUTINE init_dim_obs_pdaf

