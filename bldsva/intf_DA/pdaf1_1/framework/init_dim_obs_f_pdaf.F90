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
!init_dim_obs_f_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                        'init_dim_obs_f_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_obs_f_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: LSEIK/LETKF/LESTKF
  !
  ! The routine is called in PDAF\_lseik\_update 
  ! at the beginning of the analysis step before 
  ! the loop through all local analysis domains. 
  ! It has to determine the dimension of the 
  ! observation vector according to the current 
  ! time step for all observations required for 
  ! the analyses in the loop over all local 
  ! analysis domains on the PE-local state domain.
  !
  ! !REVISION HISTORY:
  ! 2013-02 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  !   USE mod_assimilation, &
  !        ONLY : nx, ny, local_dims, &
  !        obs_p, obs_index_p, coords_obs, local_dims_obs
  !   USE mod_parallel_pdaf, &
  !        ONLY: mype_filter, npes_filter, COMM_filter, MPI_INTEGER, &
  !        MPIerr, MPIstatus
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter
  use mod_parallel_model, &
       only: mpi_integer, model, mpi_double_precision, mpi_double, mpi_sum
  USE mod_assimilation, &
       ONLY: obs, obs_index_p, dim_obs, obs_filename, pressure_obserr_p, &
       clm_obserr_p, local_dims_obs, obs_p, global_to_local, dim_obs_p, obs_id_p, &
       longxy, latixy, longxy_obs, latixy_obs, var_id_obs, maxlon, minlon, maxlat, &
       minlat, maxix, minix, maxiy, miniy, lon_var_id, ix_var_id, lat_var_id, iy_var_id 
  Use mod_read_obs, &
       only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, read_obs_nc_multiscalar_clm_files, &
       read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, read_obs_nc_multiscalar_files, &
       z_idx_obs_nc, read_obs_nc_multi, read_obs_nc_multi_clm, clm_obs, read_obs_nc_multiscalar, &
       clmobs_lon, clmobs_lat, clmobs_layer, clmobs_dr, clm_obserr, var_id_obs_nc, dim_nx, dim_ny
  use mod_tsmp, &
#if defined CLMSA
  only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, point_obs
#else
  only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, xcoord, ycoord, zcoord, xcoord_fortran, ycoord_fortran, &
       zcoord_fortran, point_obs
#endif

#if defined CLMSA
  !kuw
  use shr_kind_mod    , only : r8 => shr_kind_r8 
  USE clmtype,    ONLY : clm3
  use decompMod , only : get_proc_bounds, get_proc_global
  USE enkf_clm_mod, only: domain_def_clm
  !kuw end
#endif

  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_lseik_update   (as U_init_dim_obs)
  ! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
  ! Called by: PDAF_letkf_update   (as U_init_dim_obs)
  !EOP

  ! *** Local variables
  INTEGER :: ierror, max_var_id
  INTEGER :: i,j,k,m,l,count  ! Counters
  logical :: is_multi_observation_files
  character (len = 110) :: current_observation_filename
  INTEGER, ALLOCATABLE :: displ(:), recv_counts(:), recv(:)
#if defined CLMSA
  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)
  integer :: begp, endp   ! per-proc beginning and ending pft indices
  integer :: begc, endc   ! per-proc beginning and ending column indices
  integer :: begl, endl   ! per-proc beginning and ending landunit indices
  integer :: begg, endg   ! per-proc gridcell ending gridcell indices
  integer :: numg         ! total number of gridcells across all processors
  integer :: numl         ! total number of landunits across all processors
  integer :: numc         ! total number of columns across all processors
  integer :: nump         ! total number of pfts across all processors
  real    :: deltax, deltay
  !real    :: deltaxy, y1 , x1, z1, x2, y2, z2, R, dist
#endif

#if defined CLMSA
  lon   => clm3%g%londeg
  lat   => clm3%g%latdeg
  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)
#endif

  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

  !  if I'm root in filter, read the nc file
  is_multi_observation_files = .true.
  if (is_multi_observation_files) then
     write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step
#if defined CLMSA
     if (mype_filter .eq. 0) then
        if(model == tag_model_parflow .and. point_obs.eq.1) then
           call read_obs_nc_multi(current_observation_filename)
        else if(model == tag_model_parflow .and. point_obs.eq.0) then
           call read_obs_nc_multiscalar_files(current_observation_filename)   
        end if
        if(model == tag_model_clm .and. point_obs.eq.1)  then
           call read_obs_nc_multi_clm(current_observation_filename)
        else if(model == tag_model_clm .and. point_obs.eq.0)  then
           call read_obs_nc_multiscalar_clm_files(current_observation_filename)
        end if
     end if
#else
     if (mype_filter.eq.0 .and. point_obs.eq.1) then 
        call read_obs_nc_multi(current_observation_filename)
     else if (mype_filter.eq.0 .and. point_obs.eq.0) then    
        call read_obs_nc_multiscalar_files(current_observation_filename)
     endif
#endif
  else
     if (mype_filter.eq.0 .and. point_obs.eq.1) then
        call read_obs_nc()
     else if (mype_filter.eq.0 .and. point_obs.eq.0) then   
        call read_obs_nc_multiscalar()
     endif
  end if

  ! broadcast dim_obs
  call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! broadcast multierr
  call mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! broadcast dim_ny and dim_nx
  if(point_obs.eq.0) then
     call mpi_bcast(dim_nx, 1, MPI_INTEGER, 0, comm_filter, ierror)
     call mpi_bcast(dim_ny, 1, MPI_INTEGER, 0, comm_filter, ierror)
  endif   

  ! allocate for non-root procs
  if (mype_filter .ne. 0) then ! for all non-master proc
#ifndef CLMSA
     if(model == tag_model_parflow) then
        allocate(idx_obs_nc(dim_obs))
        allocate(pressure_obs(dim_obs))
        if((multierr.eq.1) .and. (.not.allocated(pressure_obserr))) allocate(pressure_obserr(dim_obs))
        if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
        allocate(x_idx_obs_nc(dim_obs))
        if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
        allocate(y_idx_obs_nc(dim_obs))
        if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
        allocate(z_idx_obs_nc(dim_obs))
        if(point_obs.eq.0) allocate(var_id_obs_nc(dim_ny, dim_nx))
     end if
#endif

#if defined CLMSA
     if(model == tag_model_clm) then
        allocate(clm_obs(dim_obs))
        allocate(clmobs_lon(dim_obs))
        allocate(clmobs_lat(dim_obs))
        allocate(clmobs_dr(2))
        allocate(clmobs_layer(dim_obs))
        if(point_obs.eq.0) allocate(var_id_obs_nc(dim_ny, dim_nx))
        if(multierr.eq.1) allocate(clm_obserr(dim_obs))
     end if
#endif
  end if

#ifndef CLMSA
  ! boardcast the idx and pressure
  !if(model == tag_model_parflow) then ! for all non-master proc
  call mpi_bcast(pressure_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  if(multierr.eq.1) call mpi_bcast(pressure_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  if(point_obs.eq.0) call mpi_bcast(var_id_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  call mpi_bcast(idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  call mpi_bcast(x_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  call mpi_bcast(y_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  call mpi_bcast(z_idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
  !end if
#endif

#if defined CLMSA
  if(model == tag_model_clm) then
     call mpi_bcast(clm_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_lon, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_lat, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     !call mpi_bcast(clmobs_dr,  dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_dr,  2, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(clmobs_layer, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     if(point_obs.eq.0) call mpi_bcast(var_id_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     if(multierr.eq.1) call mpi_bcast(clm_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  end if
#endif

  ! select the obs in my domain
  dim_obs_p = 0
#ifndef CLMSA
  if(allocated(obs_id_p)) deallocate(obs_id_p)
  allocate(obs_id_p(enkf_subvecsize))
  obs_id_p(:) = 0
  if (model .eq. tag_model_parflow) then
     do i = 1, dim_obs
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              dim_obs_p = dim_obs_p + 1
              obs_id_p(j) = i
           end if
        end do
     end do
     ! Set dimension of full observation vector
     dim_obs_f = dim_obs
  end if
#endif

#if defined CLMSA
  call domain_def_clm(clmobs_lon, clmobs_lat, dim_obs, longxy, latixy, longxy_obs, latixy_obs)
  if(allocated(obs_id_p)) deallocate(obs_id_p)
  allocate(obs_id_p(endg-begg+1))
  obs_id_p(:) = 0
  if(model .eq. tag_model_clm) then
     do i = 1, dim_obs
        count = 1
        do j = begg, endg
           if((longxy_obs(i) == longxy(count)) .and. (latixy_obs(i) == latixy(count))) then
              dim_obs_p = dim_obs_p + 1
              obs_id_p(count) = i
              EXIT
           endif
           count = count + 1
        end do
     end do
     ! Set dimension of full observation vector
     dim_obs_f = dim_obs
  end if
#endif

  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  IF (ALLOCATED(obs)) DEALLOCATE(obs)
  IF (ALLOCATED(var_id_obs)) DEALLOCATE(var_id_obs)  
  ALLOCATE(obs(dim_obs))
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))
  IF(point_obs.eq.0) ALLOCATE(var_id_obs(dim_obs_p))

#ifndef CLMSA
  if (model .eq. tag_model_parflow .and. point_obs.eq.0) then
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     allocate(ix_var_id(max_var_id))
     allocate(iy_var_id(max_var_id))
     allocate(maxix(max_var_id))
     allocate(minix(max_var_id))
     allocate(maxiy(max_var_id))
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

     ! allocate pressure_obserr_p observation error for parflow run at PE-local domain 
     if((multierr.eq.1) .and. (.not.allocated(pressure_obserr_p))) allocate(pressure_obserr_p(dim_obs_p))
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
  else if (model .eq. tag_model_parflow .and. point_obs.eq.1) then
     ! allocate pressure_obserr_p observation error for parflow run at PE-local domain 
     if((multierr.eq.1) .and. (.not.allocated(pressure_obserr_p))) allocate(pressure_obserr_p(dim_obs_p))
     count = 1
     do i = 1, dim_obs
        obs(i) = pressure_obs(i)  
        ! coords_obs(1, i) = idx_obs_nc(i)
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              obs_index_p(count) = j
              obs_p(count) = pressure_obs(i)
              if(multierr.eq.1) pressure_obserr_p(count) = pressure_obserr(i)
              count = count + 1
           end if
        end do
     end do
  end if
#endif

#if defined CLMSA
  if(model .eq. tag_model_clm .and. point_obs.eq.0) then
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     allocate(lon_var_id(max_var_id))
     allocate(lat_var_id(max_var_id))
     allocate(maxlon(max_var_id))
     allocate(minlon(max_var_id))
     allocate(maxlat(max_var_id))
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

   if(multierr.eq.1) allocate(clm_obserr_p(dim_obs_p))
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
  else if(model .eq. tag_model_clm .and. point_obs.eq.1) then
     ! allocate clm_obserr_p observation error for clm run at PE-local domain 
     if(multierr.eq.1) allocate(clm_obserr_p(dim_obs_p))
     count = 1
     do i = 1, dim_obs
        obs(i) = clm_obs(i) 
        k = 1
        do j = begg,endg
           if((longxy_obs(i) == longxy(k)) .and. (latixy_obs(i) == latixy(k))) then
              obs_index_p(count) = k 
              obs_p(count) = clm_obs(i)
              if(multierr.eq.1) clm_obserr_p(count) = clm_obserr(i)
              count = count + 1
           endif
           k = k + 1
        end do
     end do
  end if
#endif

  ! allocate array of local observation dimensions with total PEs
  IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
  ALLOCATE(local_dims_obs(npes_filter))

  ! Gather array of local observation dimensions 
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
       comm_filter, ierror)

#ifndef CLMSA
  IF (model == tag_model_parflow) THEN
     !print *, "Parflow: converting xcoord to fortran"
     call C_F_POINTER(xcoord, xcoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(ycoord, ycoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(zcoord, zcoord_fortran, [enkf_subvecsize])
  ENDIF
#endif

  !  clean up the temp data from nc file
  call clean_obs_nc() 

END SUBROUTINE init_dim_obs_f_pdaf

