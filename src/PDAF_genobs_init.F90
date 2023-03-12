! Copyright (c) 2014-2021 Paul Kirchgessner
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
!BOP
!
! !ROUTINE: PDAF_genobs_init --- PDAF-internal initialization of observation generation
!
! !INTERFACE:
SUBROUTINE PDAF_genobs_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of NETF within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code based on code for ETKF
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype                ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis     ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP

! Local variables
  INTEGER :: subtype_dummy             ! Dummy variable to prevent compiler warning
  INTEGER :: param_int_dummy(dim_pint) ! Dummy variable to prevent compiler warning
  REAL :: param_real_dummy(dim_preal)  ! Dummy variable to prevent compiler warning


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Initialize dummy variables to prevent compiler warnings
  subtype_dummy = subtype
  param_int_dummy = param_int
  param_real_dummy = param_real

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++       PDAF Generator for synthetic observations       +++'
     WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
!     WRITE (*, '(/a, 4x, a)') 'PDAF', 'GENOBS configuration'

  END IF filter_pe2

  ! Set status flag
  outflag = 0

END SUBROUTINE PDAF_genobs_init
