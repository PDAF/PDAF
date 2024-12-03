!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!init_dim_obs_pdaf.F90: TSMP-PDAF implementation of routine
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
  ! Used in the filters: SEIK/EnKF/ETKF/ESTKF
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
       ONLY: mype_filter, comm_filter, npes_filter, abort_parallel, &
       mpi_integer, mpi_double_precision, mpi_in_place, mpi_sum, &
       mype_world
  USE mod_assimilation, &
       ONLY: obs_p, obs_index_p, dim_obs, obs_filename, &
       obs, &
       obs_interp_indices_p, &
       obs_interp_weights_p, &
       pressure_obserr_p, clm_obserr_p, &
       obs_pdaf2nc, &
       local_dims_obs, &
       local_disp_obs, &
       ! dim_obs_p, &
       obs_id_p, &
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
!hcp 
!CLMSA needs the physical  coordinates of the elements of state vector 
!and observation array.        
       longxy, latixy, longxy_obs, latixy_obs, &
       longxy_obs_floor, latixy_obs_floor, &
!hcp end
#endif
#endif
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
       sc_p, idx_obs_nc_p, &
#endif
#endif
       var_id_obs, maxlon, minlon, maxlat, &
       minlat, maxix, minix, maxiy, miniy, lon_var_id, ix_var_id, lat_var_id, iy_var_id, &
       screen
  USE mod_assimilation, ONLY: obs_nc2pdaf
  Use mod_read_obs, &
       only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, &
       read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, &
       z_idx_obs_nc, &
       x_idx_interp_d_obs_nc, y_idx_interp_d_obs_nc, &
       clm_obs, &
       var_id_obs_nc, dim_nx, dim_ny, &
       clmobs_lon, clmobs_lat, clmobs_layer, clmobs_dr, clm_obserr
  use mod_read_obs, only: dampfac_state_time_dependent_in
  use mod_read_obs, only: dampfac_param_time_dependent_in
  use mod_tsmp, &
      only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize
  use mod_tsmp, &
      only: nx_glob, ny_glob, nz_glob, crns_flag
  use mod_tsmp, only: da_print_obs_index
  use mod_tsmp, only: tag_model_clm
  use mod_tsmp, only: point_obs
  use mod_tsmp, only: obs_interp_switch
  use mod_tsmp, &
      only: is_dampfac_state_time_dependent, &
      dampfac_state_time_dependent, is_dampfac_param_time_dependent, dampfac_param_time_dependent
  use mod_tsmp, only: model

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  !kuw
  use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef CLMFIVE
  use GridcellType, only: grc
  use ColumnType, only : col
  ! use GetGlobalValuesMod, only: GetGlobalWrite
  ! use clm_varcon, only: nameg
  use enkf_clm_mod, only: state_clm2pdaf_p
  use enkf_clm_mod, only: clmstatevec_only_active
  use enkf_clm_mod, only: clmstatevec_max_layer
#else  
  USE clmtype,                  ONLY : clm3
#endif  
  use decompMod , only : get_proc_bounds, get_proc_global
  !kuw end
  !hcp
  !use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
  !latixy, longxy_obs, latixy_obs
  USE enkf_clm_mod, only: domain_def_clm
  USE enkf_clm_mod, only: get_interp_idx
  use enkf_clm_mod, only: clmstatevec_allcol
  !hcp end
#endif
#endif

  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  ! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector
  ! !CALLING SEQUENCE:
  ! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT    (as U_init_dim_obs)
  ! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
  ! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
  ! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
  !EOP

  ! *** Local variables
  integer :: ierror
  INTEGER :: max_var_id         ! Multi-scale DA
  INTEGER :: sum_dim_obs_p
  INTEGER :: c                ! CLM Column index
  INTEGER :: g                ! CLM Gridcell index
  INTEGER :: cg
  INTEGER :: i,j,k        ! Counters
  INTEGER :: cnt          ! Counters
  INTEGER :: cnt_interp   ! Counter for interpolation grid cells
  INTEGER :: m,l          ! Counters
  logical :: is_multi_observation_files
  character (len = 110) :: current_observation_filename
  integer :: k_cnt !,nsc !hcp
  real    :: sum_interp_weights

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)
  integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays
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
  logical :: obs_snapped     !Switch for checking multiple observation counts
  logical :: newgridcell
#endif
#endif

  character (len = 27) :: fn    !TSMP-PDAF: function name for obs_index_p output

  ! ****************************************
  ! *** Initialize observation dimension ***
  ! ****************************************

  ! Read observation file
  ! ---------------------

  ! Default: no local damping factors
  is_dampfac_state_time_dependent = 0
  is_dampfac_param_time_dependent = 0

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
  !! broadcast crns_flag
  !call mpi_bcast(crns_flag, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! broadcast dim_ny and dim_nx
  if(point_obs.eq.0) then
     call mpi_bcast(dim_nx, 1, MPI_INTEGER, 0, comm_filter, ierror)
     call mpi_bcast(dim_ny, 1, MPI_INTEGER, 0, comm_filter, ierror)
  endif
  ! broadcast damping factor flags
  call mpi_bcast(is_dampfac_state_time_dependent, 1, MPI_INTEGER, 0, comm_filter, ierror)
  call mpi_bcast(is_dampfac_param_time_dependent, 1, MPI_INTEGER, 0, comm_filter, ierror)

  ! broadcast dampfac_state_time_dependent_in
  if(is_dampfac_state_time_dependent.eq.1) then

     if (mype_filter .ne. 0) then ! for all non-master proc
       if(allocated(dampfac_state_time_dependent_in)) deallocate(dampfac_state_time_dependent_in)
       allocate(dampfac_state_time_dependent_in(1))
     end if

     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": Before setting dampfac_state_time_dependent"
     end if

     call mpi_bcast(dampfac_state_time_dependent_in, 1, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_state_time_dependent_in=", dampfac_state_time_dependent_in
     end if

     ! Set C-version of dampfac_state_time_dependent with value read from obsfile
     dampfac_state_time_dependent = dampfac_state_time_dependent_in(1)

     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_state_time_dependent=", dampfac_state_time_dependent
     end if

  end if

  ! broadcast dampfac_param_time_dependent_in
  if(is_dampfac_param_time_dependent.eq.1) then

     if (mype_filter .ne. 0) then ! for all non-master proc
       if(allocated(dampfac_param_time_dependent_in)) deallocate(dampfac_param_time_dependent_in)
       allocate(dampfac_param_time_dependent_in(1))
     end if

     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": Before setting dampfac_param_time_dependent"
     end if

     call mpi_bcast(dampfac_param_time_dependent_in, 1, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_param_time_dependent_in=", dampfac_param_time_dependent_in
     end if

     ! Set C-version of dampfac_param_time_dependent with value read from obsfile
     dampfac_param_time_dependent = dampfac_param_time_dependent_in(1)

     if (screen > 2) then
       print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_param_time_dependent=", dampfac_param_time_dependent
     end if

  end if
  
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
        if(obs_interp_switch .eq. 1) then
            if(allocated(x_idx_interp_d_obs_nc))deallocate(x_idx_interp_d_obs_nc)
            allocate(x_idx_interp_d_obs_nc(dim_obs))
            if(allocated(y_idx_interp_d_obs_nc))deallocate(y_idx_interp_d_obs_nc)
            allocate(y_idx_interp_d_obs_nc(dim_obs))
        end if
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
     if(obs_interp_switch .eq. 1) then
         call mpi_bcast(x_idx_interp_d_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
         call mpi_bcast(y_idx_interp_d_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     end if
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

      ! Interpolation of measured states: Save the indices of the
      ! nearest grid points
      if (obs_interp_switch .eq. 1) then
         ! Get the floored values for latitudes and longitudes
         call get_interp_idx(clmobs_lon, clmobs_lat, dim_obs, longxy_obs_floor, latixy_obs_floor)
      end if

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else      
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif

      ! Obtain CLM index information
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
        cnt = 1
        obs_snapped = .false.
        do g = begg, endg
            if(is_use_dr) then
                deltax = abs(lon(g)-clmobs_lon(i))
                deltay = abs(lat(g)-clmobs_lat(i))
            end if
            ! Assigning observations to grid cells according to
            ! snapping distance or index arrays
            if(((is_use_dr).and.(deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(cnt)) .and. (latixy_obs(i) == latixy(cnt)))) then
                dim_obs_p = dim_obs_p + 1
                obs_id_p(cnt) = i

                ! if (is_use_dr) then
                !   call GetGlobalWrite(g,nameg)
                ! end if

                ! Check if observation has already been snapped.
                ! Comment out if multiple grids per observation are wanted.
                if (obs_snapped) then
                  print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR Observation snapped at multiple grid cells."
                  print *, "i=", i
                  if (is_use_dr) then
                    print *, "clmobs_lon(i)=", clmobs_lon(i)
                    print *, "clmobs_lat(i)=", clmobs_lat(i)
                  end if
                  call abort_parallel()
                end if

                ! Set observation as counted
                obs_snapped = .true.
            end if
            cnt = cnt + 1
        end do
    end do
  end if
#endif
#endif

  if (screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dim_obs_p=", dim_obs_p
  end if

  ! Dimension of full observation vector
  ! ------------------------------------

  ! add and broadcast size of local observation dimensions using mpi_allreduce 
  call mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
       comm_filter, ierror) 

  ! Check sum of dimensions of PE-local observation vectors against
  ! dimension of full observation vector
  if (.not. sum_dim_obs_p == dim_obs) then
    print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR Sum of local observation dimensions"
    print *, "sum_dim_obs_p=", sum_dim_obs_p
    print *, "dim_obs=", dim_obs
    call abort_parallel()
  end if

  !  Gather local observation dimensions and displacements in arrays
  ! ----------------------------------------------------------------

  ! Allocate array of local observation dimensions
  IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
  ALLOCATE(local_dims_obs(npes_filter))

  ! Gather array of local observation dimensions
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
       comm_filter, ierror)

  ! Allocate observation displacement array local_disp_obs
  IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)
  ALLOCATE(local_disp_obs(npes_filter))

  ! Set observation displacement array local_disp_obs
  local_disp_obs(1) = 0
  do i = 2, npes_filter
     local_disp_obs(i) = local_disp_obs(i-1) + local_dims_obs(i-1)
  end do

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: local_disp_obs=", local_disp_obs
  end if

  ! Write index mapping array NetCDF->PDAF
  ! --------------------------------------
  ! Set index mapping `obs_pdaf2nc` between observation order in
  ! NetCDF input and observation order in pdaf as determined by domain
  ! decomposition.

  ! Use-case: Correct index order in loops over NetCDF-observation
  ! file input arrays.

  ! Trivial example: The order in the NetCDF file corresponds exactly
  ! to the order in the domain decomposition in PDAF, e.g. for a
  ! single PE per component model run.

  ! Non-trivial example: The first observation in the NetCDF file is
  ! not located in the domain/subgrid of the first PE. Rather, the
  ! second observation in the NetCDF file (`i=2`) is the only
  ! observation (`cnt = 1`) in the subgrid of the first PE
  ! (`mype_filter = 0`). This leads to a non-trivial index mapping,
  ! e.g. `obs_pdaf2nc(1)==2`:
  !
  ! i = 2
  ! cnt = 1
  ! mype_filter = 0
  ! 
  ! obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
  !-> obs_pdaf2nc(local_disp_obs(1)+1) = 2
  !-> obs_pdaf2nc(1) = 2

  if (allocated(obs_pdaf2nc)) deallocate(obs_pdaf2nc)
  allocate(obs_pdaf2nc(dim_obs))
  obs_pdaf2nc = 0
  if (allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
  allocate(obs_nc2pdaf(dim_obs))
  obs_nc2pdaf = 0

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  if (model .eq. tag_model_parflow) then
  if (point_obs.eq.1) then

    cnt = 1
    do i = 1, dim_obs
      do j = 1, enkf_subvecsize
        if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
          obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
          obs_nc2pdaf(i) = local_disp_obs(mype_filter+1)+cnt
          cnt = cnt + 1
        end if
      end do
    end do

  end if
  end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  if(model .eq. tag_model_clm) then
  if (point_obs.eq.1) then

    cnt = 1
    do i = 1, dim_obs
      do g = begg,endg
        newgridcell = .true.
        do c = begc,endc
          cg =   mycgridcell(c)
          if(cg .eq. g) then
            if(newgridcell) then

              if(is_use_dr) then
                deltax = abs(lon(g)-clmobs_lon(i))
                deltay = abs(lat(g)-clmobs_lat(i))
              end if

              if(((is_use_dr).and.(deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then
                obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
                obs_nc2pdaf(i) = local_disp_obs(mype_filter+1)+cnt
                cnt = cnt + 1
              end if

              newgridcell = .false.

            end if
          end if
        end do
      end do
    end do

  end if
  end if
#endif
#endif

  ! collect values from all PEs, by adding all PE-local arrays (works
  ! since only the subsection belonging to a specific PE is non-zero)
  call mpi_allreduce(MPI_IN_PLACE,obs_pdaf2nc,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
  call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: obs_pdaf2nc=", obs_pdaf2nc
  end if


  ! Write process-local observation arrays
  ! --------------------------------------
  IF (ALLOCATED(obs)) DEALLOCATE(obs)
  ALLOCATE(obs(dim_obs))
  !IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
  !ALLOCATE(obs_index(dim_obs))
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_p(dim_obs_p))
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  if(obs_interp_switch .eq. 1) then
      ! Array for storing indices from states that are interpolated to observation locations
      IF (ALLOCATED(obs_interp_indices_p)) DEALLOCATE(obs_interp_indices_p)
      ALLOCATE(obs_interp_indices_p(dim_obs_p, 4)) ! Later 8 for 3D / ParFlow
      IF (ALLOCATED(obs_interp_weights_p)) DEALLOCATE(obs_interp_weights_p)
      ALLOCATE(obs_interp_weights_p(dim_obs_p, 4)) ! Later 8 for 3D / ParFlow
  end if
  if(point_obs.eq.0) then
      IF (ALLOCATED(var_id_obs)) DEALLOCATE(var_id_obs)
      ALLOCATE(var_id_obs(dim_obs_p))
  end if

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

     if(crns_flag.eq.1) then 
        if (allocated(sc_p)) deallocate(sc_p)
        allocate(sc_p(nz_glob, dim_obs_p))
        if (allocated(idx_obs_nc_p)) deallocate(idx_obs_nc_p)
        allocate(idx_obs_nc_p(dim_obs_p))
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

     cnt = 1
     do m = 1, dim_nx
        do k = 1, dim_ny
           i = (m-1)* dim_ny + k    
           obs(i) = pressure_obs(i)  
           ! coords_obs(1, i) = idx_obs_nc(i)
           do j = 1, enkf_subvecsize
              if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
                 obs_index_p(cnt) = j
                 obs_p(cnt) = pressure_obs(i)
                 var_id_obs(cnt) = var_id_obs_nc(k,m)
                 if(multierr.eq.1) pressure_obserr_p(cnt) = pressure_obserr(i)
                 cnt = cnt + 1
              end if
           end do
        end do
     end do
  else if (point_obs.eq.1) then

     !hcp
     if(crns_flag.eq.1) then
         idx_obs_nc(:)=nx_glob*(y_idx_obs_nc(:)-1)+x_idx_obs_nc(:)
     endif
     !hcp fin
     cnt = 1
     do i = 1, dim_obs
        obs(i) = pressure_obs(i)  
        ! coords_obs(1, i) = idx_obs_nc(i)
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              !print *, j
              !obs_index(cnt) = j
              !obs(cnt) = pressure_obs(i)
              obs_index_p(cnt) = j
              obs_p(cnt) = pressure_obs(i)
              if(multierr.eq.1) pressure_obserr_p(cnt) = pressure_obserr(i)
              if(crns_flag.eq.1) then
                  idx_obs_nc_p(cnt)=idx_obs_nc(i)
                  !Allocate(sc_p(cnt)%scol_obs_in(nz_glob))       
              endif
              cnt = cnt + 1
           end if
        end do
     end do
     do i = 1, dim_obs_p
      if(crns_flag.eq.1) then 
        do k = 1, nz_glob
          k_cnt=idx_obs_nc_p(i)+(k-1)*nx_glob*ny_glob
          do j = 1, enkf_subvecsize
             if (k_cnt .eq. idx_map_subvec2state_fortran(j)) sc_p(nz_glob-k+1,i)=j
          enddo
        enddo
      endif
     enddo

     if(obs_interp_switch.eq.1) then
         ! loop over all obs and save the indices of the nearest grid
         ! points to array obs_interp_indices_p and save the distance
         ! weights to array obs_interp_weights_p (later normalized)
         cnt = 1
         do i = 1, dim_obs
             cnt_interp = 0
             do j = 1, enkf_subvecsize
                 ! First: ix and iy smaller than observation location
                 if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
                     obs_interp_indices_p(cnt, 1) = j
                     obs_interp_weights_p(cnt, 1) = sqrt(abs(x_idx_interp_d_obs_nc(i)) * abs(x_idx_interp_d_obs_nc(i)) + abs(y_idx_interp_d_obs_nc(i)) * abs(y_idx_interp_d_obs_nc(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Second: ix larger than observation location, iy smaller
                 if (idx_obs_nc(i) + 1 .eq. idx_map_subvec2state_fortran(j)) then
                     obs_interp_indices_p(cnt, 2) = j
                     obs_interp_weights_p(cnt, 2) = sqrt(abs(1.0-x_idx_interp_d_obs_nc(i)) * abs(1.0-x_idx_interp_d_obs_nc(i)) + abs(y_idx_interp_d_obs_nc(i)) * abs(y_idx_interp_d_obs_nc(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Third: ix smaller than observation location, iy larger
                 if (idx_obs_nc(i) + nx_glob .eq. idx_map_subvec2state_fortran(j)) then
                     obs_interp_indices_p(cnt, 3) = j
                     obs_interp_weights_p(cnt, 3) = sqrt(abs(x_idx_interp_d_obs_nc(i)) * abs(x_idx_interp_d_obs_nc(i)) + abs(1.0-y_idx_interp_d_obs_nc(i)) * abs(1.0-y_idx_interp_d_obs_nc(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Fourth: ix and iy larger than observation location
                 if (idx_obs_nc(i) + nx_glob + 1 .eq. idx_map_subvec2state_fortran(j)) then
                     obs_interp_indices_p(cnt, 4) = j
                     obs_interp_weights_p(cnt, 4) = sqrt(abs(1.0-x_idx_interp_d_obs_nc(i)) * abs(1.0-x_idx_interp_d_obs_nc(i)) + abs(1.0-y_idx_interp_d_obs_nc(i)) * abs(1.0-y_idx_interp_d_obs_nc(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Check if all four corners are found
                 if(cnt_interp == 4) then
                     cnt = cnt + 1
                     ! exit
                 end if
             end do
         end do

         do i = 1, dim_obs

             ! Sum of distance weights
             sum_interp_weights = sum(obs_interp_weights_p(i, :))

             do j = 1, 4
                 ! Normalize distance weights
                 obs_interp_weights_p(i, j) = obs_interp_weights_p(i, j) / sum_interp_weights
             end do
         end do

     end if

  end if
  end if
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

     cnt = 1
     do m = 1, dim_nx
        do l = 1, dim_ny
           i = (m-1)* dim_ny + l        
           obs(i) = clm_obs(i) 
           do g = begg,endg
              if((longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1))) then
                 obs_index_p(cnt) = g-begg+1
                 obs_p(cnt) = clm_obs(i)
                 var_id_obs(cnt) = var_id_obs_nc(l,m)
                 if(multierr.eq.1) clm_obserr_p(cnt) = clm_obserr(i)
                 cnt = cnt + 1
              endif
           end do
        end do
     end do
  else if(point_obs.eq.1) then

     cnt = 1

     do i = 1, dim_obs
        obs(i) = clm_obs(i)

       do g = begg,endg
         newgridcell = .true.

         do c = begc,endc

           cg =   mycgridcell(c)

           if(cg .eq. g) then

             if(newgridcell) then

               if(is_use_dr) then
                 deltax = abs(lon(g)-clmobs_lon(i))
                 deltay = abs(lat(g)-clmobs_lat(i))
               end if

               if(((is_use_dr).and.(deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then

                 ! Different settings of observation-location-index in
                 ! state vector depending on the method of state
                 ! vector assembling.
                 if(clmstatevec_allcol.eq.1) then
#ifdef CLMFIVE
                   if(clmstatevec_only_active.eq.1) then

                     ! Error if observation deeper than clmstatevec_max_layer
                     if(clmobs_layer(i) > min(clmstatevec_max_layer, col%nbedrock(c))) then
                       print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR observation layer deeper than clmstatevec_max_layer or bedrock."
                       print *, "i=", i
                       print *, "c=", c
                       print *, "clmobs_layer(i)=", clmobs_layer(i)
                       print *, "col%nbedrock(c)=", col%nbedrock(c)
                       print *, "clmstatevec_max_layer=", clmstatevec_max_layer
                       call abort_parallel()
                     end if
                     obs_index_p(cnt) = state_clm2pdaf_p(c,clmobs_layer(i))
                   else
#endif
                     obs_index_p(cnt) = c-begc+1 + ((endc-begc+1) * (clmobs_layer(i)-1))
#ifdef CLMFIVE
                   end if
#endif
                 else
                   obs_index_p(cnt) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
                 end if

                 !write(*,*) 'obs_index_p(',cnt,') is',obs_index_p(cnt)
                 obs_p(cnt) = clm_obs(i)
                 if(multierr.eq.1) clm_obserr_p(cnt) = clm_obserr(i)
                 cnt = cnt + 1
               end if

               newgridcell = .false.

             end if

           end if

         end do
       end do
     end do

     if(obs_interp_switch.eq.1) then
         ! loop over all obs and save the indices of the nearest grid
         ! points to array obs_interp_indices_p and save the distance
         ! weights to array obs_interp_weights_p (later normalized)
         cnt = 1
         do i = 1, dim_obs
             cnt_interp = 0
             do g = begg,endg
                 ! First: latitude and longitude smaller than observation location
                 if((longxy_obs_floor(i) == longxy(g-begg+1)) .and. (latixy_obs_floor(i) == latixy(g-begg+1))) then

                     obs_interp_indices_p(cnt, 1) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
                     obs_interp_weights_p(cnt, 1) = sqrt(abs(lon(g)-clmobs_lon(i)) * abs(lon(g)-clmobs_lon(i)) + abs(lat(g)-clmobs_lat(i)) * abs(lat(g)-clmobs_lat(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Second: latitude larger than observation location, longitude smaller than observation location
                 if((longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs_floor(i) == latixy(g-begg+1))) then
                     obs_interp_indices_p(cnt, 2) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
                     obs_interp_weights_p(cnt, 2) =sqrt(abs(lon(g)-clmobs_lon(i)) * abs(lon(g)-clmobs_lon(i)) + abs(lat(g)-clmobs_lat(i)) * abs(lat(g)-clmobs_lat(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Third: latitude smaller than observation location, longitude larger than observation location
                 if((longxy_obs_floor(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1))) then
                     obs_interp_indices_p(cnt, 3) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
                     obs_interp_weights_p(cnt, 3) = sqrt(abs(lon(g)-clmobs_lon(i)) * abs(lon(g)-clmobs_lon(i)) + abs(lat(g)-clmobs_lat(i)) * abs(lat(g)-clmobs_lat(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Fourth: latitude and longitude larger than observation location
                 if((longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1))) then
                     obs_interp_indices_p(cnt, 4) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
                     obs_interp_weights_p(cnt, 4) = sqrt(abs(lon(g)-clmobs_lon(i)) * abs(lon(g)-clmobs_lon(i)) + abs(lat(g)-clmobs_lat(i)) * abs(lat(g)-clmobs_lat(i)))
                     cnt_interp = cnt_interp + 1
                 end if
                 ! Check if all four corners are found
                 if(cnt_interp == 4) then
                     cnt = cnt + 1
                     ! exit
                 end if
             end do
         end do

         do i = 1, dim_obs

             ! Sum of distance weights
             sum_interp_weights = sum(obs_interp_weights_p(i, :))

             do j = 1, 4
                 ! Normalize distance weights
                  obs_interp_weights_p(i, j) = obs_interp_weights_p(i, j) / sum_interp_weights
              end do
         end do

     end if

  end if
  end if
#endif
#endif

#ifdef PDAF_DEBUG
  IF (da_print_obs_index > 0) THEN
    ! TSMP-PDAF: For debug runs, output the state vector in files
    WRITE(fn, "(a,i5.5,a,i5.5,a)") "obs_index_p_", mype_world, ".", step, ".txt"
    OPEN(unit=71, file=fn, action="write")
    DO i = 1, dim_obs_p
      WRITE (71,"(i10)") obs_index_p(i)
    END DO
    CLOSE(71)
  END IF
#endif


  !  clean up the temp data from nc file
  ! ------------------------------------
  call clean_obs_nc()

END SUBROUTINE init_dim_obs_pdaf

