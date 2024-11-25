!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz,Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!init_dim_obs_l_pdaf.F90: TSMP-PDAF implementation of routine
!                         'init_dim_obs_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_obs_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: LSEIK/LETKF/LESTKF
  !
  ! The routine is called during the loop over
  ! all local analysis domains. It has to set 
  ! the dimension of the local observation vector 
  ! for the current local analysis domain.
  !
  ! !REVISION HISTORY:
  ! 2013-02 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  !   USE mod_assimilation, &
  !        ONLY: nx, ny, local_dims, &
  !        cradius, coords_obs, coords_l, obs_index_p, obs_index_l
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter
  USE mod_assimilation, &
       ONLY: cradius, obs_index_l, dim_obs, obs_p, distance, obs_index_p, &
       dim_state, dim_obs_p, &
       longxy, latixy, longxy_obs, latixy_obs
  USE mod_assimilation, &
       ONLY: lon_var_id, ix_var_id, lat_var_id, iy_var_id
  USE mod_read_obs, &
       ONLY: x_idx_obs_nc, y_idx_obs_nc, z_idx_obs_nc, idx_obs_nc, clmobs_lon, &
       clmobs_lat, var_id_obs_nc, dim_nx, dim_ny 
  USE mod_tsmp, &
#if defined CLMSA
  ONLY: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, point_obs, model
#else
  ONLY: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, nx_glob, ny_glob, nz_glob, &
       xcoord, ycoord, zcoord, &
       xcoord_fortran, ycoord_fortran, zcoord_fortran, &
       point_obs, model
#endif

  USE, INTRINSIC :: iso_c_binding, ONLY: C_F_POINTER

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
  ! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
  ! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
  !EOP

  ! local variables
  INTEGER :: i, j, k, m, cnt ! Counters
  !  INTEGER :: idx, ix, iy, ix1, iy1
  REAL :: dist ! Distance between observation and analysis domain
  LOGICAL, ALLOCATABLE :: log_var_id(:) ! logical variable ID for setting location observation vector using remote sensing data
  INTEGER  :: domain_p_coord   ! Current local analysis domain for coord arrays

  !kuw
  integer :: dx,dy, max_var_id, ierror
  integer :: obsind(dim_obs)
  real    :: obsdist(dim_obs)
  ! kuw end

  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************

  ! Setting fortran pointer to ParFlow-coordinate arrays
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
!!#if (defined PARFLOW_STAND_ALONE || defined COUP_OAS_PFL)
  IF (model == tag_model_parflow) THEN
     !print *, "Parflow: converting xcoord to fortran"
     call C_F_POINTER(xcoord, xcoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(ycoord, ycoord_fortran, [enkf_subvecsize])
     call C_F_POINTER(zcoord, zcoord_fortran, [enkf_subvecsize])
  ENDIF

  ! Index for local analysis domain `domain_p` in coordinate array
  ! that only spans `enkf_subvecsize`.
  !
  ! Necessary condition: `domain_p` is initialized as an index of the
  ! process-local state dimension in `init_n_domains_pdaf`
  !
  ! Necessary condition II: the full state vector consists of sections
  ! of size `enkf_subvecsize`, where each section corresponds to a
  ! single coordinate array.
  domain_p_coord = MODULO(domain_p, enkf_subvecsize)
#endif
#endif

  ! Count observations within cradius

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
  obsind    = 0
  obsdist   = 0.0
  dim_obs_l = 0
  if(point_obs.eq.0) then
     if(model == tag_model_parflow) THEN
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     allocate(log_var_id(max_var_id))
     log_var_id(:) = .TRUE.

     do m = 1, dim_nx
        do k = 1, dim_ny   
           i = (m-1)* dim_ny + k
           do j = 1, max_var_id
              if(log_var_id(j) .and. var_id_obs_nc(k,m) == j) then
                 dx = abs(x_idx_obs_nc(i) - int(xcoord_fortran(domain_p_coord))-1)
                 dy = abs(y_idx_obs_nc(i) - int(ycoord_fortran(domain_p_coord))-1)
                 dist = sqrt(real(dx)**2 + real(dy)**2)
                 !obsdist(i) = dist
                 if (dist <= real(cradius) .AND. dist > 0) then
                    dim_obs_l = dim_obs_l + 1
                    obsind(i) = 1
                    log_var_id(j) = .FALSE.
                    dx = abs(ix_var_id(j) - int(xcoord_fortran(domain_p_coord))-1)
                    dy = abs(iy_var_id(j) - int(ycoord_fortran(domain_p_coord))-1)
                    obsdist(i) = sqrt(real(dx)**2 + real(dy)**2)
                 else if(dist == 0) then
                    dim_obs_l = dim_obs_l + 1
                    obsind(i) = 1
                    log_var_id(j) = .FALSE.
                    obsdist(i) = dist  
                 end if
              end if
           end do
        end do
     end do
     end if
  else   
     if(model == tag_model_parflow) THEN
        do i = 1,dim_obs
           dx = abs(x_idx_obs_nc(i) - int(xcoord_fortran(domain_p_coord))-1)
           dy = abs(y_idx_obs_nc(i) - int(ycoord_fortran(domain_p_coord))-1)
           dist = sqrt(real(dx)**2 + real(dy)**2)
           obsdist(i) = dist
           if (dist <= real(cradius)) then
              dim_obs_l = dim_obs_l + 1
              obsind(i) = 1
           end if
        end do
     endif
  endif
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
  obsind    = 0
  obsdist   = 0.0
  dim_obs_l = 0
  if(point_obs.eq.0) then
     max_var_id = MAXVAL(var_id_obs_nc(:,:))
     allocate(log_var_id(max_var_id))
     log_var_id(:) = .TRUE.

     if(model == tag_model_clm) THEN
     do m = 1, dim_nx
        do k = 1, dim_ny   
           i = (m-1)* dim_ny + k
           do j = 1, max_var_id
              if(log_var_id(j) .and. var_id_obs_nc(k,m) == j) then
                 dx = abs(longxy_obs(i) - longxy(domain_p))
                 dy = abs(latixy_obs(i) - latixy(domain_p))
                 dist = sqrt(real(dx)**2 + real(dy)**2)
                 !obsdist(i) = dist
                 if(dist == 0) then
                    dim_obs_l = dim_obs_l + 1
                    obsind(i) = 1
                    log_var_id(j) = .FALSE.
                    obsdist(i) = dist
                    EXIT
                 end if
              end if
           end do
        enddo
     enddo

     do m = 1, dim_nx
        do k = 1, dim_ny   
           i = (m-1)* dim_ny + k
           do j = 1, max_var_id
              if(log_var_id(j) .and. var_id_obs_nc(k,m) == j) then
                 dx = abs(longxy_obs(i) - longxy(domain_p))
                 dy = abs(latixy_obs(i) - latixy(domain_p))
                 dist = sqrt(real(dx)**2 + real(dy)**2)
                 !obsdist(i) = dist
                 if (dist <= real(cradius)) then
                    dim_obs_l = dim_obs_l + 1
                    obsind(i) = 1
                    log_var_id(j) = .FALSE.
                    dx = abs(lon_var_id(j) - longxy(domain_p))
                    dy = abs(lat_var_id(j) - latixy(domain_p))
                    obsdist(i) = sqrt(real(dx)**2 + real(dy)**2)
                 end if
              end if
           end do
        enddo
     enddo
     end if
  else 
     if(model == tag_model_clm) THEN
     do i = 1,dim_obs
        dx = abs(longxy_obs(i) - longxy(domain_p))
        dy = abs(latixy_obs(i) - latixy(domain_p))
        dist = sqrt(real(dx)**2 + real(dy)**2)
        obsdist(i) = dist
        if (dist <= real(cradius)) then
           dim_obs_l = dim_obs_l + 1
           obsind(i) = 1
        end if
     end do
     end if 
  end if
#endif
#endif  
!------------------------------------------------------------------------

  ! kuw: allocate and determine local observation index and distance
  !#ifndef CLMSA
  ! Initialize index array for local observations in full observed vector
  IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
  IF(dim_obs_l /= 0) ALLOCATE(obs_index_l(dim_obs_l))
  ! The array holds the distance of an observation from local analysis domain.
  IF (ALLOCATED(distance)) DEALLOCATE(distance)
  IF(dim_obs_l /= 0) ALLOCATE(distance(dim_obs_l))

  cnt = 1
  do i = 1,dim_obs
     if(obsind(i).eq.1) then
        obs_index_l(cnt) = i
        distance(cnt)    = obsdist(i)
        !print *,'mype_filter distance(cnt)  ', mype_filter, distance(cnt) 
        cnt = cnt + 1
     end if
  end do
  
  ! if allocated than deallocate logical variable ID log_var_id for setting location
  ! observation vector using remote sensing data
  IF (ALLOCATED(log_var_id)) DEALLOCATE(log_var_id)

END SUBROUTINE init_dim_obs_l_pdaf

