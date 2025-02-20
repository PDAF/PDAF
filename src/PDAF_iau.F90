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
!
!> Module for IAU
!!
!! This module provides variables and subroutines for incremental
!! updating in PDAF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_iau

  USE mpi
  USE PDAF_mod_filter, &
       ONLY: debug, screen

  IMPLICIT NONE
  SAVE

  ! Parameters set by user in PDAF_init_iau
  INTEGER :: type_iau=0           !< Control IAU: (0) run with direct updates; (1) use IAU
  INTEGER :: nsteps_iau=0         !< Number if time steps over which IAU is performed

  ! Internal parameters
  LOGICAL :: iau_now=.FALSE.      !< Whether to do IAU in the current forecast phase

  ! Ensemble array for incremental updating
  REAL, TARGET, ALLOCATABLE :: ens_iau(:,:) !< Matrix holding increment for IAU
  REAL, ALLOCATABLE :: iau_weight(:)        !< Vector holding the increment weights

CONTAINS

!-------------------------------------------------------------------------------
!> Initialize IAU
!!
!! Initialize parameters for IAU and allocate the IAU ensemble array.
!! The IAU ensemble array only exists on processes that are model tasks.
!!
  SUBROUTINE PDAF_iau_init(type_iau_in, nsteps_iau_in, flag)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &
         ONLY: mype_world, mype_model, filterpe, task_id, dim_ens_task
    USE PDAF_mod_filter, &
         ONLY: dim_p

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in) :: type_iau_in     !< Type of IAU, (0) no IAU
    INTEGER, INTENT(in) :: nsteps_iau_in   !< number of time steps in IAU
    INTEGER, INTENT(out) :: flag           !< Status flag

    ! *** Local variables ***
    INTEGER :: allocstat                   ! Status for allocate
    

    ! ********************************
    ! *** Initialize IAU variables ***
    ! ********************************

    ! Set status flag
    flag = 0

    ! Store IAU type
    type_iau = type_iau_in
    
    doiau: IF (type_iau > 0) THEN

       IF (nsteps_iau_in > 0) THEN
          nsteps_iau = nsteps_iau_in
       ELSE
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(30): Invalid value for IAU steps; required is >0!'
          flag = 30
       END IF

       ! Screen output
       IF (mype_world == 0 .AND. screen > 0) THEN
          WRITE (*, '(/a)') 'PDAF: Activate IAU'
          WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'IAU time steps', nsteps_iau
       END IF


    ! ******************************
    ! *** Initialize IAU weights ***
    ! ******************************

       CALL PDAF_iau_set_weights(type_iau, nsteps_iau)

       ! Screen output
       IF (mype_world == 0 .AND. screen > 0) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '! Note: IAU starts after first analysis. To start from'
          WRITE (*, '(a, 5x, a)') 'PDAF', '! initial time set initial increment with PDAF_iau_init_inc'
       END IF


    ! *****************************
    ! *** Allocate IAU ensemble ***
    ! *****************************

       on_filterpe: IF (filterpe) THEN

          ! Allocate task-local increment ensemble for IAU
          ! if the filter task is also model task
          IF (task_id > 0) THEN

             ALLOCATE(ens_iau(dim_p, dim_ens_task), stat = allocstat)
             IF (allocstat /= 0) THEN
                WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ENS_IAU'
                flag = 20
             END IF

             IF (screen > 2) WRITE (*,*) 'PDAF: PDAF_iau_init - allocate ENS_IAU of size ', &
                  dim_ens_task, ' on pe(m) ', mype_model, ' of model task ',task_id

             ! count allocated memory
             CALL PDAF_memcount(1,'r',dim_p*dim_ens_task)

          ELSE
             ALLOCATE(ens_iau(1,1), stat = allocstat)
          END IF

       ELSE

          ! Allocate task-local increment ensemble for IAU for all model tasks
          ALLOCATE(ens_iau(dim_p, dim_ens_task), stat = allocstat)
          IF (allocstat /= 0) THEN
             WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ENS_IAU'
             flag = 20
          END IF

          IF (screen > 2) WRITE (*,*) 'PDAF: PDAF_iau_init - allocate ENS_IAU of size ', &
               dim_ens_task, ' on pe(m) ', mype_model, ' of model task ',task_id
          
           ! count allocated memory
          CALL PDAF_memcount(1,'r',dim_p*dim_ens_task)

       END IF on_filterpe

    ELSE doiau

       ALLOCATE(ens_iau(1,1), stat = allocstat)

    END IF doiau
    
    ! Initialize ens_iau
    ens_iau = 0.0

  END SUBROUTINE PDAF_iau_init


!-------------------------------------------------------------------------------
!> Change IAU settings
!!
!! This routine allows to modify the IAU type and number of IAU
!! time steps during a run
!!
  SUBROUTINE PDAF_iau_reset(type_iau_in, nsteps_iau_in, flag)

    USE PDAF_mod_filtermpi, &
         ONLY: mype_world, mype_model, filterpe, task_id, dim_ens_task
    USE PDAF_mod_filter, &
         ONLY: dim_p

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in) :: type_iau_in     !< Type of IAU, (0) no IAU
    INTEGER, INTENT(in) :: nsteps_iau_in   !< number of time steps in IAU
    INTEGER, INTENT(out) :: flag           !< Status flag

    ! *** Local variables ***
    INTEGER :: allocstat                   ! Status for allocate
    

    ! ***************************
    ! *** Reset IAU variables ***
    ! ***************************

    ! Set status flag
    flag = 0

    ! Store IAU type
    type_iau = type_iau_in
    
    doiau: IF (type_iau > 0) THEN

       IF (nsteps_iau_in > 0) THEN
          nsteps_iau = nsteps_iau_in
       ELSE
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(30): Invalid value for IAU steps; required is >0!'
          flag = 30
       END IF

       ! Screen output
       IF (mype_world == 0 .AND. screen > 0) THEN
          WRITE (*, '(/a)') 'PDAF: Reset IAU'
          WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'IAU time steps', nsteps_iau
       END IF


    ! ******************************
    ! *** Initialize IAU weights ***
    ! ******************************

       CALL PDAF_iau_set_weights(type_iau, nsteps_iau)

    END IF doiau

  END SUBROUTINE PDAF_iau_reset

!-------------------------------------------------------------------------------
!> Set IAU weight vector
!!
!! This routine allocates the IAU weight vector and initializes
!! it according to the IAU type
!!
  SUBROUTINE PDAF_iau_set_weights(type_iau, nsteps_iau)

    USE PDAF_mod_filtermpi, &
         ONLY: mype_world, task_id, mype_model

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in) :: type_iau     !< Type of IAU, (0) no IAU
    INTEGER, INTENT(in) :: nsteps_iau   !< number of time steps in IAU

    ! *** Local variables ***
    INTEGER :: i          ! Counter
    INTEGER :: halfsteps  ! half of IAU period
    REAL :: norm          ! normalization of weight

    IF (nsteps_iau > 0) THEN

       ! Screen output
       IF (mype_world == 0 .AND. screen > 0) THEN
          WRITE (*, '(a, 5x, a, i3)') 'PDAF', 'IAU type', type_iau
          IF (type_iau == 1) THEN
             WRITE (*, '(a,8x, a, i3)') 'PDAF', '-- constant IAU weight' 
          END IF
       END IF

       ! *** Allocate IAU weight vector ***
       IF (ALLOCATED(iau_weight)) DEALLOCATE(iau_weight)
       ALLOCATE(iau_weight(nsteps_iau))

       ! *** Set weights ***
       IF (type_iau == 1) THEN
          ! Uniform weight
          iau_weight = 1.0 / REAL(nsteps_iau)

       ELSEIF (type_iau==2) THEN
          ! Linear with maximum in middle of IAU period (analogous to NEMO)

          ! Compute normalization
          norm = 0.0
          IF (MOD(nsteps_iau, 2)==0) THEN
             ! even number of time steps
             halfsteps = nsteps_iau/2
             DO i = 1, halfsteps
                norm = norm + REAL(i)
             END DO
             norm = 2.0 * norm
          ELSE
             ! odd number of time steps
             halfsteps = (nsteps_iau+1)/2
             DO i = 1, halfsteps-1
                norm = norm + REAL(i)
             END DO
             norm = 2.0 * norm + REAL (halfsteps)
          END IF
          norm = 1.0 / norm

          ! Initialize weights
          DO i = 1, halfsteps-1
             iau_weight(i) = REAL(i) * norm
          END DO
          DO i = halfsteps, nsteps_iau
             iau_weight(i) = REAL(nsteps_iau-i+1) * norm
          END DO
   
       END IF

       IF (debug > 0 .AND. task_id==1 .AND. mype_model==0) &
            WRITE (*,*) '++ PDAF-debug IAU:', debug, 'IAU weights', iau_weight
    END IF

  END SUBROUTINE PDAF_iau_set_weights


!-------------------------------------------------------------------------------
!> Deallocate IAU ensemble arrays
!!
  SUBROUTINE PDAF_iau_dealloc()

    IMPLICIT NONE

    IF (ALLOCATED(ens_iau)) DEALLOCATE(ens_iau)
    IF (ALLOCATED(iau_weight)) DEALLOCATE(iau_weight)

  END SUBROUTINE PDAF_iau_dealloc


!-------------------------------------------------------------------------------
!> Initialize the increment array from user side
!!
!! If the IAU should start from the initial time of a run,
!! the user needs to set the increment ensemble because PDAF
!! cannot compute it at this time.
!!
  SUBROUTINE PDAF_iau_init_inc(dim_p, dim_ens_l, ens_inc, flag)

    USE PDAF_mod_filtermpi, &
         ONLY: filterpe, task_id, dim_ens_task, mype_couple, &
         npes_couple, mype_world

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                    !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens_l                !< Task-local size of ensemble
    REAL, INTENT(in) :: ens_inc(dim_p, dim_ens_l)   !< PE-local increment ensemble
    INTEGER, INTENT(out) :: flag                    !< Status flag


    ! **********************
    ! *** Initialization ***
    ! **********************

    doiau: IF (type_iau > 0) THEN

       ! Activate IAU (it might have been deactivated at first forecast)
       iau_now = .TRUE.


    ! ********************************************************
    ! *** Initialize increment ensemble on each model task ***
    ! ********************************************************

       IF (dim_ens_task == dim_ens_l) THEN

          ! Screen output
          IF (mype_world == 0 .AND. screen > 0) THEN
             WRITE (*, '(a,5x,a)') 'PDAF', 'Setting initial increment for IAU - IAU activated'
          END IF

          on_filterpe: IF (filterpe) THEN

             ! Initialize increment if filter task is also model task
             IF (task_id > 0) THEN
                ens_iau(:, 1:dim_ens_task) = ens_inc(:, 1:dim_ens_task)
             END IF

          ELSE

             ! Initialize increment on all model tasks
             ens_iau(:, 1:dim_ens_task) = ens_inc(:, 1:dim_ens_task)

          END IF on_filterpe

       ELSE
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(30): Task-local increment ensemble size is inconsistent!'
          flag = 30
       END IF
    END IF doiau

  END SUBROUTINE PDAF_iau_init_inc


!-------------------------------------------------------------------------------
!> Store forecast ensemble in IAU ensemble array
!!
!! After the forecast phase we store the forecast ensemble into
!! the iau ensemble array. 
!!
  SUBROUTINE PDAF_iau_update_ens(ens)

    USE PDAF_mod_filtermpi, &
         ONLY: filterpe, task_id, dim_ens_task, mype_couple, &
         npes_couple

    IMPLICIT NONE

    ! *** Arguments ***
    REAL, INTENT(inout) :: ens(:, :)   !< PE-local state ensemble


    ! **********************
    ! *** Initialization ***
    ! **********************

    ! Activate IAU (it might have been deactivated at first forecast)
    iau_now = .TRUE.


    ! *******************************
    ! *** Store forecast ensemble ***
    ! *******************************

    doiau: IF (type_iau > 0) THEN

       on_filterpe: IF (filterpe) THEN

          ! Store ensemble if filter task is also model task
          IF (task_id > 0) THEN
             ens_iau(:, 1:dim_ens_task) = ens(:, 1:dim_ens_task)
          END IF

       ELSE
          ! Store ensemble on all model tasks
          ens_iau(:, 1:dim_ens_task) = ens(:, 1:dim_ens_task)

       END IF on_filterpe

    END IF doiau

  END SUBROUTINE PDAF_iau_update_ens



!-------------------------------------------------------------------------------
!> Compute DA increment and store it in IAU ensemble array
!!
!! After the analysis step we compute the DA increment and 
!! store it in the IAU ensemble array.
!!
  SUBROUTINE PDAF_iau_update_inc(ens_ana)

    USE PDAF_mod_filtermpi, &
         ONLY: filterpe, task_id, dim_ens_task, mype_couple, &
         npes_couple, mype_world

    IMPLICIT NONE

    ! *** Arguments ***
    REAL, INTENT(inout) :: ens_ana(:, :)   !< PE-local analysis ensemble


    ! ********************************************
    ! *** Compute and store analysis increment ***
    ! ********************************************

    doiau: IF (type_iau > 0 .AND. iau_now) THEN

       on_filterpe: IF (filterpe) THEN

          ! Store ensemble if filter task is also model task
          IF (task_id > 0) THEN
             ens_iau(:, 1:dim_ens_task) = ens_ana(:, 1:dim_ens_task) - ens_iau(:, 1:dim_ens_task)
          END IF

       ELSE
          ! Compute increment on all model tasks
          ens_iau(:, 1:dim_ens_task) = ens_ana(:, 1:dim_ens_task) - ens_iau(:, 1:dim_ens_task)

       END IF on_filterpe

    END IF doiau

  END SUBROUTINE PDAF_iau_update_inc


!-------------------------------------------------------------------------------
!> Apply IAU increment
!!
!! During the forecast phase add the IAU increments according to its IAU weight
!!
  SUBROUTINE PDAF_iau_apply_inc(step, dim_p, dim_ens_task, ens, U_collect_state, U_distribute_state)

    USE PDAF_mod_filtermpi, &
         ONLY: filterpe, task_id, mype_model

    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in) :: step         !< Time step
    INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
    INTEGER, INTENT(in) :: dim_ens_task !< Ensemble size of model task
    REAL, INTENT(inout) :: ens(:, :)    !< PE-local state ensemble

    ! *** External subroutines ***
    !  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_collect_state, &      !< Routine to collect a state vector
         U_distribute_state             !< Routine to distribute a state vector

    ! *** Local variables ***
    INTEGER :: member          ! Counter
    LOGICAL :: do_iau_inc      ! Flag whether the IAU update is done


    ! *******************************
    ! *** Store forecast ensemble ***
    ! *******************************

    doiau: IF (type_iau>0 .AND. iau_now) THEN

       ! Determine whether increment is added
       do_iau_inc = .FALSE.
       IF (step <= nsteps_iau) THEN
          IF (filterpe) THEN
             IF (task_id > 0) do_iau_inc = .TRUE.
          ELSE
             do_iau_inc = .TRUE.
          END IF
       END IF

       apply_iau: IF (do_iau_inc) THEN

          ! Screen output
          IF (task_id==1 .AND. mype_model==0 .AND. screen>0) &
               WRITE (*,'(a, 5x, a)') 'PDAF', 'Apply IAU'

          DO member = 1, dim_ens_task

             IF (debug>0 .AND. task_id>0 .AND. mype_model==0) THEN
                WRITE (*,*) '++ PDAF-debug IAU:', debug, &
                     'apply IAU, weight', iau_weight(step)
                WRITE (*,*) '++ PDAF-debug IAU:', debug, ' task: ', task_id, &
                     'call collect_state for IAU'
             END IF

             ! Store model fields in state vector
             CALL U_collect_state(dim_p, ens(:, member))

             ! Add increment
             ens(:, member) = ens(:, member) + iau_weight(step)*ens_iau(:, member)

             ! Write state vector back to model fields
             IF (debug>0 .AND. task_id>0 .AND. mype_model==0) &
                  WRITE (*,*) '++ PDAF-debug IAU:', debug, ' task: ', task_id, &
                  'call distribute_state for IAU'

             CALL U_distribute_state(dim_p, ens(:, member))

          END DO

       END IF apply_iau

    END IF doiau

  END SUBROUTINE PDAF_iau_apply_inc

END MODULE PDAF_iau
