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
!$Id: PDAF-D_seek_init.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_init --- PDAF-internal initialization of SEEK filter
!
! !INTERFACE:
SUBROUTINE PDAF_seek_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of SEEK within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - Print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_eof, forget, &
       int_rediag, epsilon, type_forget, dim_ens

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

  ! Whether incremental updating is performed
  if (dim_pint>=4) THEN
     incremental = param_int(4)
     IF (param_int(4) /= 0 .AND. param_int(4) /= 1) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(10): Invalid setting for incremental updating!'
        outflag = 10
     END IF
  END IF

  ! Interval to perform rediagonalization
  IF (dim_pint >= 3) THEN
     int_rediag = param_int(3)
     IF (int_rediag <= 0) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(101): Invalid setting for re-diag interval!'
        outflag = 101
     END IF
  END IF

  ! epsilon for approximated TLM evolution in SEEK
  IF (dim_pint >= 2 .AND. subtype /= 5) THEN
     epsilon = param_real(2)
     IF (epsilon <= 0.0) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(102): Invalid setting for epsilon in SEEK!'
        outflag = 102
     END IF
  END IF

  ! For fixed basis SEEK do not perform rediagonalization
  IF (subtype == 2 .OR. subtype == 3) THEN
     int_rediag = 0
  END IF

  ! Type of forgetting factor - not a choice for SEEK
  type_forget = 0
     
  ! Special for SEEK: Initialize number of modes
  dim_eof = dim_ens

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .FALSE.

  ! Initialize flag for fixed-basis filters
  IF (subtype == 2 .OR. subtype == 3) THEN
     fixedbasis = .TRUE.
  ELSE
     fixedbasis = .FALSE.
  END IF


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/a, 5x, a)') 'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++                  SEEK Filter                   +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++    Pham et al., J. Mar. Syst. 16 (1998) 323    +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++          This implementation follows           +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++      Nerger et al., Tellus 57A (2005) 715      +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'SEEK configuration'
     WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard SEEK with unit modes'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> SEEK with non-unit modes'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> fixed basis filter with update of matrix U'
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> no re-diagonalization of VUV^T'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> fixed basis filter & no update of matrix U'
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> no re-diagonalization of VUV^T'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 2
     END IF
     IF (incremental == 1) &
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
     IF (type_forget == 0) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor!'
        outflag = 8
     ENDIF
     WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> number of EOFs:', dim_eof
     IF (subtype /= 5) THEN
        IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) THEN
           IF (int_rediag == 1) THEN
              WRITE (*, '(a, 10x, a, i4, a)') 'PDAF', 'Re-diag at each analysis step'
           ELSE
              WRITE (*, '(a, 10x, a, i4, a)') 'PDAF', 'Re-diag at each ', int_rediag, &
                   '-th analysis step'
           END IF
        END IF
     ELSE
        IF (int_rediag == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform re-diagonalization'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', 'No re-diagonalization'
        END IF
     END IF
  END IF filter_pe2

END SUBROUTINE PDAF_seek_init
