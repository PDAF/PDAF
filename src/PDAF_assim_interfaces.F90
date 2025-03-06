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
!> Interface definitions for PDAF
!!
!! Module providing interface definition of the PDAF routines that
!! are called from the model code.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-05 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - Split for better coding overview
!! * Other revisions - see repository log
MODULE PDAF_assim_interfaces

  INTERFACE 
     SUBROUTINE PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
          param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, in_filterpe, U_init_ens, in_screen, &
          flag)
       INTEGER, INTENT(in) :: filtertype     ! Type of filter
       INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
       INTEGER, INTENT(in) :: stepnull       ! Initial time step of assimilation
       INTEGER, INTENT(in) :: dim_pint       ! Number of integer parameters
       INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
       INTEGER, INTENT(in) :: dim_preal      ! Number of real parameter 
       REAL, INTENT(inout) :: param_real(dim_preal) ! Real parameter array
       INTEGER, INTENT(in) :: COMM_model     ! Model communicator
       INTEGER, INTENT(in) :: COMM_couple    ! Coupling communicator
       INTEGER, INTENT(in) :: COMM_filter    ! Filter communicator
       INTEGER, INTENT(in) :: task_id        ! Id of my ensemble task
       INTEGER, INTENT(in) :: n_modeltasks   ! Number of parallel model tasks
       LOGICAL, INTENT(in) :: in_filterpe    ! Is my PE a filter-PE?
       INTEGER, INTENT(in) :: in_screen      ! Control screen output:
       INTEGER, INTENT(out):: flag           ! Status flag, 0: no error, error codes:
       EXTERNAL :: U_init_ens  ! User-supplied routine for ensemble initialization
     END SUBROUTINE PDAF_init
  END INTERFACE

END MODULE PDAF_assim_interfaces
