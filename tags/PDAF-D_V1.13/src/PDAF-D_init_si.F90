! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_init_si.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_init_si --- Interface routine to initialize PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_init_si(filtertype, subtype, stepnull, param_int, dim_pint, &
     param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
     task_id, n_modeltasks, in_filterpe, in_screen, outflag)

! !DESCRIPTION:
! Interface routine for initialization of PDAF.
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using pre-defined routine names.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  ! For valid and default values see PDAF-D_mod_filter.F90
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
                                        ! (0) none, (1) some, default, (2) extensive
  INTEGER, INTENT(out):: outflag        ! Status flag

! ! Names of external subroutines 
  EXTERNAL :: init_ens_pdaf        ! Routine to initialize the ensemble or modes
                                   ! for the assimilation algorithm.

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_init
!EOP


! *********************************************
! *** Call the full initialization routine  ***
! *********************************************

  CALL PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
       param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, in_filterpe, init_ens_pdaf, in_screen, &
       outflag)


END SUBROUTINE PDAF_init_si
