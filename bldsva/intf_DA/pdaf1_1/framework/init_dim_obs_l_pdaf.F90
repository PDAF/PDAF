!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz,Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!init_dim_obs_l_pdaf.F90: TerrSysMP-PDAF implementation of routine
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
!        local_range, coords_obs, coords_l, obs_index_p, obs_index_l
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter
  USE mod_parallel_model, &
       ONLY: model
  USE mod_assimilation, &
       ONLY: local_range, obs_index_l, dim_obs, obs_p, distance, obs_index_p, &
       dim_state, xcoord_fortran_g, ycoord_fortran_g, zcoord_fortran_g, dim_obs_p, &
       longxy, latixy, longxy_obs, latixy_obs  
  USE mod_read_obs, &
       ONLY: x_idx_obs_nc, y_idx_obs_nc, z_idx_obs_nc, idx_obs_nc, clmobs_lon, &
       clmobs_lat 
  USE mod_tsmp, &
#if defined CLMSA
       ONLY: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm
#else
       ONLY: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
       tag_model_clm, nx_glob, ny_glob, nz_glob, xcoord_fortran, ycoord_fortran, &
       zcoord_fortran
#endif

  USE, INTRINSIC :: iso_c_binding

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
  INTEGER :: i, j, cnt ! Counters
  REAL :: xcoord_obs, ycoord_obs, zcoord_obs ! store the coordinates for each observation point 
  REAL :: limits_x(2), limits_y(2), limits_z(2) ! Coordinate limits for observation domain
!  INTEGER :: idx, ix, iy, ix1, iy1
  REAL :: dist ! Distance between observation and analysis domain

  !kuw
  integer :: dx,dy
  integer :: obsind(dim_obs)
  real    :: obsdist(dim_obs)
  ! kuw end

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

! Count observations within local_range
#ifndef CLMSA
  obsind    = 0
  obsdist   = 0.0
  dim_obs_l = 0
  do i = 1,dim_obs
    dx = abs(x_idx_obs_nc(i) - int(xcoord_fortran(domain_p))-1)
    dy = abs(y_idx_obs_nc(i) - int(ycoord_fortran(domain_p))-1)
    dist = sqrt(real(dx)**2 + real(dy)**2)
    obsdist(i) = dist
    if (dist <= real(local_range)) then
      dim_obs_l = dim_obs_l + 1
      obsind(i) = 1
    end if 
  end do
#endif
  ! kuw end

! Count observations within local_range
! for clm stand alone model
#if defined CLMSA
  obsind    = 0
  obsdist   = 0.0
  dim_obs_l = 0
  do i = 1,dim_obs
    dx = abs(longxy_obs(i) - longxy(domain_p))
    dy = abs(latixy_obs(i) - latixy(domain_p))
    dist = sqrt(real(dx)**2 + real(dy)**2)
    obsdist(i) = dist
    if (dist <= real(local_range)) then
      dim_obs_l = dim_obs_l + 1
      obsind(i) = 1
    end if 
  end do
#endif

!print*,'obs_index_l ', obs_index_l

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
      cnt = cnt + 1
    end if
  end do

END SUBROUTINE init_dim_obs_l_pdaf

