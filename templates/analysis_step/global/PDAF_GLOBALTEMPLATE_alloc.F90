! Copyright (c) 2004-2024 Lars Nerger
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
!
!>  Perform allocation of arrays for GLOBALTEMPLATE.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_GLOBALTEMPLATE_alloc(subtype, outflag)

  USE mpi
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: screen, incremental, dim_ens, dim_p, dim_lag, &
       state, state_inc, eofU, eofV, sens
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, filterpe, dim_ens_l, task_id, &
       COMM_couple

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
                                        !   (allocated arrays might depend on this)
  INTEGER, INTENT(out):: outflag        ! Status flag

! *** Local variables ***
  INTEGER :: allocstat                  ! Status for allocate


! ******************************
! *** Allocate filter fields ***
! ******************************

  on_filterpe: IF (filterpe) THEN
     ! On Filter-processes: allocate all arrays and full ensemble matrix

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                          +++
! +++ Arrays that need to be allocated for filterpe==.true. are:         +++
! +++  - state (size dim_p)                                              +++
! +++  - eofU (usually size (dim_ens, dim_ens) or (dim_ens-1, dim_ens-1) +++
! +++         if no transform matrix is used allocate with size (1,1))   +++
! +++  - eofV (size (dim_p, dim_ens))                                    +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_p)

     ALLOCATE(eofU(dim_ens, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofU'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_ens**2)

     ! Allocate full ensemble on filter-PEs
     ALLOCATE(eofV(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_ens)

     IF (screen > 2) WRITE (*,*) 'PDAF: GLOBALTEMPLATE_alloc - allocate eofV of size ', &
          dim_ens, ' on pe(f) ', mype

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                          +++
! +++ The following are optional                                         +++
! +++ - state_inc (state increment in case of incremental updating       +++
! +++ - sens (smoother ensemble array holding past ensembles; this is    +++
! +++         only need if a smoother is implemented and dim_lag>0)      +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     IF (incremental == 1) THEN
        ALLOCATE(state_inc(dim_p), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE_INC'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(1,'r',dim_p)

        state_inc = 0.0
     END IF

     ! Allocate array for past ensembles for smoothing on filter-PEs
     IF (dim_lag > 0) THEN
        ALLOCATE(sens(dim_p, dim_ens, dim_lag), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of sens'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_ens * dim_lag)
     ELSE
        ALLOCATE(sens(1, 1, 1), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,*) 'PDAF-ERROR(20): error in allocation of sens'
           outflag = 20
        END IF
     END IF

  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                          +++
! +++ The array that needs to be allocated for filterpe==.false. is:     +++
! +++  - eofV (size (dim_p, dim_ens_l); this array hold the sub-ensemble +++
! +++         forecasted by the model task to which the process belongs) +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        ALLOCATE(eofV(dim_p, dim_ens_l), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV on model-pe'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_ens_l)

        IF (screen > 2) WRITE (*,*) 'PDAF: GLOBALTEMPLATE_alloc - allocate eofV of size ', &
             dim_ens_l, ' on pe(m) ', mype_model, ' of model task ',task_id
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_GLOBALTEMPLATE_alloc
