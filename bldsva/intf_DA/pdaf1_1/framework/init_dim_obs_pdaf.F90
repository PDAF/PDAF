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
  !        ONLY : nx, ny, local_dims, obs_p, obs_index_p
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, comm_filter, npes_filter
  use mod_parallel_model, &
       only: mpi_integer, model, mpi_double_precision, mpi_in_place, mpi_sum, &
       mype_world
  USE mod_assimilation, &
       ONLY: obs_p, obs_index_p, dim_obs, obs_filename, dim_state_p, &
       pressure_obserr_p, clm_obserr_p, obs_nc2pdaf, &
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
!hcp 
!CLMSA needs the physical  coordinates of the elements of state vector 
!and observation array.        
       longxy, latixy, longxy_obs, latixy_obs, &
!hcp end
#endif
#endif
       screen
  Use mod_read_obs, &
       only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, &
       read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, &
       z_idx_obs_nc, clm_obs, &
       clmobs_lon, clmobs_lat,clmobs_layer, clmobs_dr, clm_obserr
  use mod_tsmp, &
      only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
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
  INTEGER :: i,j,k,count   ! Counters
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
  !end if
#endif
#endif

  ! Generate CLM index arrays from lon/lat values
  ! ---------------------------------------------
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
      call domain_def_clm(clmobs_lon, clmobs_lat, dim_obs, longxy, latixy, longxy_obs, latixy_obs)
  end if
#endif
#endif
!hcp end

  ! Number of observations in process-local domain
  ! ----------------------------------------------
  dim_obs_p = 0

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  if (model .eq. tag_model_parflow) then
     do i = 1, dim_obs
        do j = 1, enkf_subvecsize
           if (idx_obs_nc(i) .eq. idx_map_subvec2state_fortran(j)) then
              dim_obs_p = dim_obs_p + 1
           end if
        end do
     end do
     ! saving the size of local observation vector to variable dim_state_p
     !dim_state_p = dim_obs_p
  end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  if(model .eq. tag_model_clm) then

     lon   => clm3%g%londeg
     lat   => clm3%g%latdeg
     call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
     call get_proc_global(numg, numl, numc, nump)

     do i = 1, dim_obs
       do j = begg, endg
           deltax = abs(lon(j)-clmobs_lon(i))
           deltay = abs(lat(j)-clmobs_lat(i))
           if((deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))) then
              dim_obs_p = dim_obs_p + 1
           end if
        end do
     end do
     ! saving the size of local observation vector to variable dim_state_p
     !dim_state_p = dim_obs_p
  end if
#endif
#endif

  if (screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dim_obs_p=", dim_obs_p
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
  
  !IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
  !IF (ALLOCATED(obs)) DEALLOCATE(obs)
  !ALLOCATE(obs_index(dim_obs_p))
  !ALLOCATE(obs(dim_obs_p))
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))

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
     count = 1
     do i = 1, dim_obs
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

     count = 1
     do i = 1, dim_obs
       do j = begg,endg
           deltax = abs(lon(j)-clmobs_lon(i))
           deltay = abs(lat(j)-clmobs_lat(i))
           if((deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))) then
              !obs_index_p(count) = j + (size(lon) * (clmobs_layer(i)-1))
              !obs_index_p(count) = j + ((endg-begg+1) * (clmobs_layer(i)-1))
              obs_index_p(count) = j-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
              !write(*,*) 'obs_index_p(',count,') is',obs_index_p(count)
              obs_p(count) = clm_obs(i)
              if(multierr.eq.1) clm_obserr_p(count) = clm_obserr(i)
              obs_nc2pdaf(local_dis(mype_filter+1)+count) = i
              count = count + 1
           end if
        end do
     end do
  end if
  call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
#endif
#endif

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: obs_nc2pdaf=", obs_nc2pdaf
  end if

  !  clean up the temp data from nc file
  ! ------------------------------------
  deallocate(local_dis)
  call clean_obs_nc()

END SUBROUTINE init_dim_obs_pdaf

