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
!> General allocation routine of PDAF
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code from restructuring
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_alloc(dim_p, dim_ens, dim_ens_task, dim_es, dim_bias_p, &
     dim_lag, statetask, incremental, outflag)

  USE mpi
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: screen, state, &
       state_inc, Ainv, ens, sens, bias
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, mype_couple, filterpe, task_id, &
       COMM_couple

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p          !< Size of state vector
  INTEGER, INTENT(in) :: dim_ens        !< Ensemble size
  INTEGER, INTENT(in) :: dim_ens_task   !< Ensemble size handled by a model task
  INTEGER, INTENT(in) :: dim_es         !< Dimension of error space (size of Ainv)
  INTEGER, INTENT(in) :: dim_bias_p     !< Size of bias vector
  INTEGER, INTENT(in) :: dim_lag        !< Smoother lag
  INTEGER, INTENT(in) :: statetask      !< Task ID forecasting a single state
  INTEGER, INTENT(in) :: incremental    !< >0 to allocate state_inc_p
  INTEGER, INTENT(inout):: outflag      !< Status flag

! *** local variables ***
  INTEGER :: allocstat                  ! Status for allocate


! ****************************
! *** Allocate PDAF arrays ***
! ****************************

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_p)

     IF (incremental > 0) THEN
        ALLOCATE(state_inc(dim_p), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE_INC'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(1,'r',dim_p)

        state_inc = 0.0
     ELSE
        ALLOCATE(state_inc(1), stat = allocstat)
     END IF

     IF (dim_es > 1) THEN
        ALLOCATE(Ainv(dim_es, dim_es), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of Ainv'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(1, 'r', dim_es**2)
     ELSE
        ALLOCATE(Ainv(1, 1), stat = allocstat)
     END IF

     ! Allocate full ensemble on filter-PEs
     ALLOCATE(ens(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_ens)

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

     IF (dim_bias_p > 0) THEN
        ALLOCATE(bias(dim_bias_p), stat = allocstat)
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_bias_p)
        
        ! initialize bias field
        bias = 0.0
     ELSE
        ALLOCATE(bias(1), stat = allocstat)
     ENDIF
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of BIAS'
        outflag = 20
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate ens of size ', &
          dim_ens, ' on pe(f) ', mype
     
  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        ALLOCATE(ens(dim_p, dim_ens_task), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens on model-pe'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_ens_task)

        IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate ens of size ', &
             dim_ens_task, ' on pe(m) ', mype_model, ' of model task ',task_id

        ! Some of the model-PEs may integrate a central state
        IF (mype_couple+1 == statetask) THEN
           ALLOCATE(state(dim_p), stat = allocstat)
           IF (allocstat /= 0) THEN
              WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
              outflag = 20
           END IF
              
           IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate state of size ', &
                dim_p, ' on pe(m) ', mype_model, ' of model task ',task_id

           ! count allocated memory
           CALL PDAF_memcount(2, 'r', dim_p)
        END IF
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_alloc
