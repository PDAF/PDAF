! Copyright (c) 2014-2016 Paul Kirchgessner
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
! !ROUTINE: PDAF_NETF_init --- PDAF-internal initialization of NETF
!
! !INTERFACE:
SUBROUTINE PDAF_NETF_init(subtype, param_int, dim_pint, param_real, dim_preal, &
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
  USE PDAF_mod_filter, &
       ONLY: incremental, forget, &
       type_forget, type_trans, dim_lag

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


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Size of lag considered for smoother
  IF (dim_pint>=3) THEN
     IF (param_int(3) > 0) THEN
        dim_lag = param_int(3)
     ELSE
        dim_lag = 0
     END IF
  END IF

  ! Store type for forgetting factor
  IF (dim_pint >= 5) THEN
     type_forget = param_int(5)
  END IF

  ! Type of ensemble transformation
  IF (dim_pint >= 6) THEN     
     type_trans = param_int(6)
  END IF


  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++      Nonlinear Ensemble Transform Filter (NETF)       +++'
     WRITE(*, '(a)')  'PDAF    +++                                                       +++'
     WRITE(*, '(a)')  'PDAF    +++                         by                            +++'
     WRITE(*, '(a)')  'PDAF    +++ J. Toedter, B. Ahrens, Mon. Wea. Rev. 143 (2015) 1347 +++'
     WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'NETF configuration'
     WRITE (*, '(a, 11x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') '--> NETF '
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 3
     END IF
     IF (type_trans == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
     ELSE IF (type_trans == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
     ELSE
        WRITE (*,'(/1x, a/)') &
             'PDAF-ERROR(9): Invalid setting for ensemble transformation!'
        outflag = 9
     END IF
     IF (incremental == 1) &
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
     IF (type_forget == 0) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
     ELSEIF (type_forget == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use global adaptive forgetting factor'
     ELSEIF (type_forget == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use local adaptive forgetting factors'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor!'
        outflag = 8
     ENDIF

  END IF filter_pe2

END SUBROUTINE PDAF_NETF_init
