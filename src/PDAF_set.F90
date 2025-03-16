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
!> Module providing PDAF_set routines
!!
MODULE PDAF_set

CONTAINS
!> Set root MPI communicator for PDAF
!!
!! Helper routine for PDAF.
!! This routine allows to set the overall (world)
!! MPI communicator for PDAF. By default this is 
!! MPI_COMM_WORLD. However, in the case that not
!! all processes call PDAF or participate in the DA
!! and forecasting - as for example in case of
!! an IO server that uses separate MPI tasks - a
!! separate communicator can be set.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_set_comm_pdaf(in_COMM_pdaf)

  USE PDAF_mod_parallel, &
       ONLY: isset_comm_pdaf, COMM_pdaf
  USE PDAF_mod_core, &
       ONLY: debug

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(in) :: in_COMM_pdaf    !< MPI communicator for PDAF

! *** Set ensemble member ***

  COMM_pdaf = in_COMM_pdaf

  isset_comm_pdaf = .true.

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'Set user-defined communicator for COMM_PDAF'

END SUBROUTINE PDAF_set_comm_pdaf

!--------------------------------------------------------------------------
!> Set debugging flag
!!
!! This routine set the debug flag for PDAF. 
!! One can set the flag dependent on the local analysis
!! domain, the MPI rank, or the OpenMP thread ID, or
!! and combination of them.
!!
!! For debugval>0 additional information is written by
!! the PDAF routine to stdout. One should activate the 
!! debugging before calling some selected routine(s) and
!! deactivate it with debugval=0 afterwards. This allows 
!! for a targeted checking of the functionality.
!!
!! Note: The debugging of PDAF is independent of that 
!! for PDAF-OMI.
!!
!!  This is a core routine of PDAF and
!!  should not be changed by the user   !
!!
!! __Revision history:__
!! * 2022-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_set_debug_flag(debugval)

  USE PDAF_mod_core, &
       ONLY: debug
  USE PDAF_mod_parallel, &
       ONLY: mype, filterpe, mype_world

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: debugval          !< Value for debugging flag

! *** Local variables ***
  INTEGER, SAVE :: debug_save=0            ! Store previous debug value to indicate deactivation


! *** Set debugging flag ***

  debug = debugval

  ! Print debug information
  IF (debug>0) THEN
     IF (filterpe) THEN
        WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_filter', mype, 'activate', debug
     ELSE
        WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_world', mype_world, 'activate', debug
     END IF
  ELSE 
     IF (debug_save>0 .AND. debug==0) THEN
        IF (filterpe) THEN
           WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_filter', mype, 'deactivate'
        ELSE
           WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_world', mype_world, 'deactivate'
        END IF
     END IF
  END IF

  ! Save current value of debug
  debug_save = debug

END SUBROUTINE PDAF_set_debug_flag


!--------------------------------------------------------------------------
!> Set pointer to ensemble array
!!
!! Routine to set the pointer to the PDAF-internal ensemble array.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2016-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_set_ens_pointer(ens_ptr, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_core, &
       ONLY: ens

  IMPLICIT NONE

! *** Arguments ***
  REAL, POINTER, INTENT(out) :: ens_ptr(:,:)  !< Pointer to ensemble array
  INTEGER, INTENT(out)       :: status        !< Status flag

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(ens)) THEN
     ens_ptr => ens

     status = 0
  ELSE
     status = 1
  END IF
  
END SUBROUTINE PDAF_set_ens_pointer


!--------------------------------------------------------------------------
!> PDAF_set_iparam --- Set integer parameter for PDAF
!!
!! This routine lets the user set the value of a
!! method-specific integer parameter. The routine
!! simply calls PDAF_set_iparam_filters, which 
!! includes the method-specific calls.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
SUBROUTINE PDAF_set_iparam(id, value, flag)

  USE PDAF_utils_filters, &
       ONLY: PDAF_set_iparam_filters
  USE PDAF_mod_core, &
       ONLY: subtype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: id       !< Index of parameter
  INTEGER, INTENT(in) :: value    !< Parameter value
  INTEGER, INTENT(inout) :: flag  !< Status flag: 0 for no error

! *** Local variable ***
  INTEGER :: inflag               ! Input flag value


! *******************************************
! *** Set filter-specific parameter value ***
! *******************************************

  IF (subtype_filter > -1) CALL PDAF_set_iparam_filters(id, value, inflag)

  ! Incremenat flag
  flag = flag + inflag

END SUBROUTINE PDAF_set_iparam


!--------------------------------------------------------------------------
!> Set ensemble index to e.g. force an analysis step
!!
!! Helper routine for PDAF.
!! The routine allows to overwrite the member index
!! of the ensemble state that is currently integrated.
!! The typical use is to set it to local_dim_ens to force
!! the analysis step at the next call to PDAF_put_state.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2021-02 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_set_memberid(memberid)

  USE PDAF_mod_core, &
       ONLY: member

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(inout) :: memberid    !< Index in the local ensemble


! *** Set ensemble member ***

  member = memberid

END SUBROUTINE PDAF_set_memberid


!--------------------------------------------------------------------------
!> Set offline mode of PDAF
!!
!! Helper routine for PDAF.
!!
!! This routine allows to activate the offline
!! mode of PDAF. Thus, the functionality of
!! PDAF to integrate an emsemble will be 
!! deactivated.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2023-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_set_offline_mode(screen)

  USE PDAF_mod_core, &
       ONLY: offline_mode, subtype_filter
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_utils_filters, &
       ONLY: PDAF_configinfo_filters

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(in) :: screen

! *** Set offline mode ***

  offline_mode = .true.

  IF (mype_world == 0 .AND. screen > 0) THEN
     WRITE (*,'(/a,4x,a)') 'PDAF','Activate PDAF offline mode'

     ! Print configuration info
     CALL PDAF_configinfo_filters(subtype_filter, 1)
  END IF

END SUBROUTINE PDAF_set_offline_mode


!--------------------------------------------------------------------------
!> PDAF_set_rparam --- Set real parameter for PDAF
!!
!! This routine lets the user set the value of a
!! method-specific real parameter. The routine
!! simply calls PDAF_set_rparam_filters, which 
!! includes the method-specific calls.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
SUBROUTINE PDAF_set_rparam(id, value, flag)

  USE PDAF_utils_filters, &
       ONLY: PDAF_set_rparam_filters
  USE PDAF_mod_core, &
       ONLY: subtype_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: id       !< Index of parameter
  REAL, INTENT(in) :: value       !< Parameter value
  INTEGER, INTENT(inout) :: flag  !< Status flag: 0 for no error

! *** Local variable ***
  INTEGER :: inflag               ! Input flag value


! *******************************************
! *** Set filter-specific parameter value ***
! *******************************************

  IF (subtype_filter > -1) CALL PDAF_set_rparam_filters(id, value, inflag)

  ! Incremenat flag
  flag = flag + inflag

END SUBROUTINE PDAF_set_rparam


!--------------------------------------------------------------------------
!> Set seedset for random number generation
!!
!! Helper routine for PDAF.
!! The routine allows to set the seedset index that
!! is used in PDAF_generate_rndmat. Values between
!! 1 and 20 are allowed.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2021-02 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_set_seedset(seedset_in)

  USE PDAF_mod_core, &
       ONLY: seedset, new_seedset

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(in) :: seedset_in    !< Seedset index (1-20)


! *** Set ensemble member ***

  IF (seedset_in>0 .AND. seedset_in<21) THEN
     seedset = seedset_in
     new_seedset = .TRUE.
  ELSE
     write (*,*) 'PDAF-ERROR: PDAF_set_seedset - Invalid value for seedset'
  END IF

END SUBROUTINE PDAF_set_seedset


!--------------------------------------------------------------------------
!> Set pointer to smoother ensemble
!!
!! Routine to set the pointer to the PDAF-internal smoother ensemble array
!! and to set the value of the maximum lah to be considered.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-05 - Lars Nerger - Initial code
!! *Other revisions - see repository log
!!
SUBROUTINE PDAF_set_smootherens(sens_point, maxlag, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_core, &
       ONLY: sens, cnt_maxlag, dim_lag

  IMPLICIT NONE

! *** Arguments ***
  REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  !< Pointer to smoother array
  INTEGER, INTENT(in)        :: maxlag  !< Number of past timesteps in sens
  INTEGER, INTENT(out)       :: status  !< Status flag, 
                                        !< 0: no error, 1: maxlag too large

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(sens)) THEN
     sens_point => sens

     status = 0
  END IF

  
! **************************************
! *** Set number of initialized lags ***
! **************************************

  IF (maxlag <= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     cnt_maxlag = maxlag
     
     status = 0
  ELSE
     ! Maxlag is larger than allocated smoother ensemble array
     cnt_maxlag = dim_lag

     status = 1
  END IF

END SUBROUTINE PDAF_set_smootherens

!--------------------------------------------------------------------------
!> Manually reset forgetting factor
!!
!! Helper routine for PDAF.
!! The routine allows to manually set the forgetting
!! factor to a new value. Usually this should be called
!! in assimilate_pdaf before calling the analysis step
!! routine.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-05 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_reset_forget(forget_in)

  IMPLICIT NONE
  
! *** Arguments ***
  REAL,INTENT(in) :: forget_in    !< New value of forgetting factor

! *** Local variables ***
  INTEGER :: flag      ! Status flag


! *** Set forgetting factor ***

  CALL PDAF_set_rparam(1, forget_in, flag)

END SUBROUTINE PDAF_reset_forget

END MODULE PDAF_set


!--------------------------------------------------------------------------
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
  USE PDAF_mod_core, &
       ONLY: screen, dim_ens, dim_p, ens
  USE PDAF_mod_parallel, &
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


!--------------------------------------------------------------------------
!> Reset state dimension and re-allocate state and ensemble
!!
!! Reset state dimension and re-allocate the state vector
!! and ensemble array.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_reset_dim_p(dim_p_in, outflag)

  USE mpi
  USE PDAF_mod_core, &
       ONLY: screen, dim_ens, dim_p, &
       state, ens !, type_iau, ens_iau
  USE PDAF_mod_parallel, &
       ONLY: mype, mype_model, filterpe, dim_ens_l, task_id, &
       COMM_couple, dim_ens_l, dim_ens_task
  USE PDAF_iau, &
       ONLY: type_iau, ens_iau, state_iau

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p_in       !< Sub-type of filter
  INTEGER, INTENT(inout):: outflag      !< Status flag

! *** local variables ***
  INTEGER :: allocstat                  ! Status for allocate


! ************************************
! *** Reset state vector dimension ***
! ************************************

  dim_p = dim_p_in


! *********************************
! *** Re-allocate filter fields ***
! *********************************

  ! Initialize status flag
  outflag = 0
  
  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     IF (ALLOCATED(state)) DEALLOCATE(state)
     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in re-allocation of STATE'
        outflag = 20
     END IF
 
     IF (task_id > 0 .AND. type_iau > 0) THEN
        IF (ALLOCATED(ens_iau)) DEALLOCATE(ens_iau)
        ALLOCATE(ens_iau(dim_p, dim_ens_task), stat = allocstat)
        IF (ALLOCATED(state_iau)) DEALLOCATE(state_iau)
        ALLOCATE(state_iau(dim_p), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ENS_IAU'
           outflag = 20
        END IF

        ens_iau = 0.0
        state_iau = 0.0
     ELSE
        IF (ALLOCATED(ens_iau)) DEALLOCATE(ens_iau)
        ALLOCATE(ens_iau(1, 1), stat = allocstat)
        IF (ALLOCATED(state_iau)) DEALLOCATE(state_iau)
        ALLOCATE(state_iau(1), stat = allocstat)
     END IF

     ! Allocate full ensemble on filter-PEs
     IF (ALLOCATED(ens)) DEALLOCATE(ens)
     ALLOCATE(ens(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens'
        outflag = 20
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_p - re-allocate ens of size ', &
          dim_ens, ' on pe(f) ', mype

  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        IF (ALLOCATED(ens)) DEALLOCATE(ens)
        ALLOCATE(ens(dim_p, dim_ens_l), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens on model-pe'
           outflag = 20
        END IF

        IF (type_iau > 0) THEN
           IF (ALLOCATED(ens_iau)) DEALLOCATE(ens_iau)
           ALLOCATE(ens_iau(dim_p, dim_ens_task), stat = allocstat)
           IF (ALLOCATED(state_iau)) DEALLOCATE(state_iau)
           ALLOCATE(state_iau(dim_p), stat = allocstat)
           IF (allocstat /= 0) THEN
              WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ENS_IAU'
              outflag = 20
           END IF

           ens_iau = 0.0
           state_iau = 0.0
        ELSE
           IF (ALLOCATED(ens_iau)) DEALLOCATE(ens_iau)
           ALLOCATE(ens_iau(1, 1), stat = allocstat)
           IF (ALLOCATED(state_iau)) DEALLOCATE(state_iau)
           ALLOCATE(state_iau(1), stat = allocstat)
        END IF

        IF (screen > 2) WRITE (*,*) 'PDAF: reset_dim_p - re-allocate ens of size ', &
             dim_ens_l, ' on pe(m) ', mype_model, ' of model task ',task_id
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_reset_dim_p
