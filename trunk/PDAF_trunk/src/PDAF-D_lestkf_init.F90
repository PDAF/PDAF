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
!$Id: PDAF-D_lestkf_init.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lestkf_init --- PDAF-internal initialization of LESTKF
!
! !INTERFACE:
SUBROUTINE PDAF_lestkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of LESTKF within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, rank, forget, &
       type_forget, type_trans, type_sqrt, dim_lag

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

  ! Whether incremental updating is performed
  if (dim_pint>=4) THEN
     incremental = param_int(4)
     IF (param_int(4) /= 0) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(10): ESTKF does not yet support incremental updating!'
        outflag = 10
     END IF
  END IF

  ! Rank of initial covariance matrix
  rank = dim_ens - 1

  ! Store type for forgetting factor
  IF (dim_pint >= 5) THEN
     type_forget = param_int(5)
  END IF

  ! Type of ensemble transformation
  IF (dim_pint >= 6) THEN     
     type_trans = param_int(6)
  END IF

  ! Store type of transform matrix square-root
  IF (dim_pint >= 7) THEN
     type_sqrt = param_int(7)
  END IF

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

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
  
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++  Local Error Subspace Transform Kalman Filter   +++'
     WRITE(*, '(a)')  'PDAF    +++                    (LESTKF)                     +++'
     WRITE(*, '(a)')  'PDAF    +++                                                 +++'
     WRITE(*, '(a)')  'PDAF    +++ Domain-localized implementation of the ESTKF by +++'
     WRITE(*, '(a)')  'PDAF    +++  Nerger et al., Mon. Wea. Rev. 140 (2012) 2335  +++'
     WRITE(*, '(a)')  'PDAF    +++           doi:10.1175/MWR-D-11-00102.1          +++'
     WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'LESTKF configuration'
     WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard LESTKF'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LESTKF with fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LESTKF with fixed state covariance matrix'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF ERROR(2): No valid sub type!'
        outflag = 2
     END IF
     IF (type_trans == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with deterministic Omega'
     ELSE IF (type_trans == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with random orthonormal Omega'
     ELSE IF (type_trans == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with product Omega'
     ELSE
        WRITE (*,'(/5x, a/)') &
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
     IF (dim_lag > 0) &
          WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
     WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens

  END IF filter_pe2

END SUBROUTINE PDAF_lestkf_init
