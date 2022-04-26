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
#ifdef CLMSA
       ONLY: obs_p, obs_index_p, dim_obs, obs_filename, dim_state_p, &
       pressure_obserr_p, clm_obserr_p, obs_nc2pdaf, screen, &
!hcp 
!CLMSA needs the physical  coordinates of the elements of state vector 
!and observation array.        
       longxy, latixy, longxy_obs, latixy_obs
!hcp end
#else
       ONLY: obs_p, obs_index_p, dim_obs, obs_filename, dim_state_p, &
       pressure_obserr_p, clm_obserr_p, obs_nc2pdaf, screen
#endif
  Use mod_read_obs, &
       only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, &
       read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, &
       z_idx_obs_nc, read_obs_nc_multi, read_obs_nc_multi_clm, clm_obs, &
       clmobs_lon, clmobs_lat,clmobs_layer, clmobs_dr, clm_obserr
  use mod_tsmp, &
#if defined CLMSA
       only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, point_obs
#else
       only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, point_obs
#endif


#if defined CLMSA
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

  IMPLICIT NONE
  ! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector
  integer :: ierror
  ! !CALLING SEQUENCE:
  ! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
  ! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
  ! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
  ! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
  ! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
  !EOP

  ! *** Local variables
  INTEGER :: i,j,k,count !,jj                  ! Counters
  logical :: flag, searching
  logical :: is_multi_observation_files
  character (len = 110) :: current_observation_filename
  integer,allocatable :: local_dis(:),local_dim(:)
#if defined CLMSA
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
  real    :: deltax, deltay !, deltaxy, y1 , x1, z1, x2, y2, z2, R, deltaxy_max
#endif

#if defined CLMSA
  lon   => clm3%g%londeg
  lat   => clm3%g%latdeg
  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)
#endif

  ! ****************************************
  ! *** Initialize observation dimension ***
  ! ****************************************

  !  if I'm root in filter, read the nc file
  is_multi_observation_files = .true.
  if (is_multi_observation_files) then
     write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step
#if defined CLMSA
     if (mype_filter .eq. 0) then
         if (screen > 2) then
             print *, "TSMP-PDAF mype(w)=", mype_world, ": read_obs_nc, CLMSA"
             print *, "TSMP-PDAF mype(w)=", mype_world, ": model=", model
         end if
        if(model == tag_model_parflow) then
           call read_obs_nc_multi(current_observation_filename)
        end if
        if(model == tag_model_clm)  then
           call read_obs_nc_multi_clm(current_observation_filename)
        end if
     end if
#else
     !hcp: This need to be changed later in LST DA with clm-pfl
     if (mype_filter==0 .and. screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": read_obs_nc, coupled"
         print *, "TSMP-PDAF mype(w)=", mype_world, ": model=", model
     end if
     if (mype_filter.eq.0) call read_obs_nc_multi(current_observation_filename)
#endif
  else
     if (mype_filter.eq.0) call read_obs_nc()
  end if

  if (mype_filter==0 .and. screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": broadcast obs vars"
      print *, "TSMP-PDAF mype(w)=", mype_world, ": dim_obs=", dim_obs
  end if
  ! broadcast dim_obs
  call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)
  ! broadcast multierr
  call mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)


  ! allocate for non-root procs
  if (mype_filter .ne. 0) then ! for all non-master proc
#ifndef CLMSA
     !if(model == tag_model_parflow) then
        if(allocated(idx_obs_nc))deallocate(idx_obs_nc)
        allocate(idx_obs_nc(dim_obs))
        if(allocated(pressure_obs))deallocate(pressure_obs)
        allocate(pressure_obs(dim_obs))
!        if((multierr.eq.1) .and. (.not.allocated(pressure_obserr))) allocate(pressure_obserr(dim_obs))
        if (multierr.eq.1) then
             if (allocated(pressure_obserr)) deallocate(pressure_obserr)
             allocate(pressure_obserr(dim_obs))
        endif
        if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
        allocate(x_idx_obs_nc(dim_obs))
        if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
        allocate(y_idx_obs_nc(dim_obs))
        if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
        allocate(z_idx_obs_nc(dim_obs))
    !end if
#endif

#if defined CLMSA
     if(model == tag_model_clm) then
        if(allocated(clm_obs))deallocate(clm_obs)
        allocate(clm_obs(dim_obs))
        if(allocated(clmobs_lon))deallocate(clmobs_lon)
        allocate(clmobs_lon(dim_obs))
        if(allocated(clmobs_lat))deallocate(clmobs_lat)
        allocate(clmobs_lat(dim_obs))
        if(allocated(clmobs_dr))deallocate(clmobs_dr)
        allocate(clmobs_dr(2))
        if(allocated(clmobs_layer))deallocate(clmobs_layer)
        allocate(clmobs_layer(dim_obs))
        if(multierr.eq.1) then 
           if(allocated(clm_obserr))deallocate(clm_obserr)               
           allocate(clm_obserr(dim_obs))
        end if 
     end if
#endif
  end if

#ifndef CLMSA
  ! boardcast the idx and pressure
  !if(model == tag_model_parflow) then ! for all non-master proc
     call mpi_bcast(pressure_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     if(multierr.eq.1) call mpi_bcast(pressure_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
     call mpi_bcast(idx_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
     ! broadcast xyz indices
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
     if(multierr.eq.1) call mpi_bcast(clm_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  end if
#endif

  ! select the obs in my domain
  dim_obs_p = 0
!hcp
!use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
!latixy, longxy_obs, latixy_obs 
#ifdef CLMSA
  call domain_def_clm(clmobs_lon, clmobs_lat, dim_obs, longxy, latixy, longxy_obs, latixy_obs)
#endif
!hcp end
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

#if defined CLMSA
  if(model .eq. tag_model_clm) then
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

  print *, "init_dim_obs_pdaf: dim_obs_p is", dim_obs_p

  !IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
  !IF (ALLOCATED(obs)) DEALLOCATE(obs)
  !ALLOCATE(obs_index(dim_obs_p))
  !ALLOCATE(obs(dim_obs_p))
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))

  ! allocate index for mapping between observations in nc input and sorted by pdaf
  if (allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
  allocate(obs_nc2pdaf(dim_obs))
  obs_nc2pdaf = 0
  allocate(local_dim(npes_filter))
  allocate(local_dis(npes_filter))
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dim, 1, MPI_INTEGER, comm_filter, ierror)
  local_dis(1) = 0
  do i = 2, npes_filter
     local_dis(i) = local_dis(i-1) + local_dim(i-1)
  end do


#ifndef CLMSA
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
            
#if defined CLMSA
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

  !  clean up the temp data from nc file
  call clean_obs_nc()
  deallocate(local_dim)
  deallocate(local_dis)

END SUBROUTINE init_dim_obs_pdaf

