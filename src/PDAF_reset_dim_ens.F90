! Copyright (c) 2004-2025 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!> Reset ensemble and re-allocate ensemble
!!
!! This routine resets the ensemble size and re-allocates
!! the ensemble array.
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
SUBROUTINE PDAF_reset_dim_ens(dim_ens_in, outflag)

  USE mpi
  USE PDAF_mod_filter, &
       ONLY: screen, dim_ens, dim_p, ens
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, filterpe, dim_ens_l, task_id, &
       COMM_couple

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_ens_in     ! Sub-type of filter
  INTEGER, INTENT(inout):: outflag      ! Status flag

! *** local variables ***
  INTEGER :: allocstat                  ! Status for allocate


! ***************************
! *** Reset ensemble size ***
! ***************************

  ! Initialize status flag
  outflag = 0

  dim_ens = dim_ens_in

  IF (dim_ens < 1) THEN
     WRITE (*,'(/5x,a/)') 'PDAF-ERROR(6): Invalid ensemble size!'
     outflag = 20
  END IF


! *********************************
! *** Re-allocate filter fields ***
! *********************************

  on_filterpe: IF (filterpe .AND. outflag==0) THEN

     ! Re-allocate full ensemble on filter-PEs
     IF (ALLOCATED(ens)) DEALLOCATE(ens)
     ALLOCATE(ens(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens'
        outflag = 20
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_ens - re-allocate ens of size ', &
          dim_ens, ' on pe(f) ', mype

  ELSEIF (outflag==0) THEN on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Re-allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        IF (ALLOCATED(ens)) DEALLOCATE(ens)
        ALLOCATE(ens(dim_p, dim_ens_l), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens on model-pe'
           outflag = 20
        END IF

        IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_ens - re-allocate ens of size ', &
             dim_ens_l, ' on pe(m) ', mype_model, ' of model task ',task_id
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_reset_dim_ens
