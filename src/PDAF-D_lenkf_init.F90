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
!BOP
!
! !ROUTINE: PDAF_lenkf_init --- PDAF-internal initialization of LEnKF
!
! !INTERFACE:
SUBROUTINE PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of EnKF within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2015-12 - Lars Nerger - Initial code by copying PDAF_enkf_init
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, forget, &
       rank_ana_enkf, type_forget, dim_lag

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(inout) :: subtype             ! Sub-type of filter
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
  ! We do not have incremental updating for EnKF!
  if (dim_pint>=4) THEN
     incremental = param_int(4)
     IF (param_int(4) /= 0 .AND. param_int(4) /= 1) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(10): EnKF does not yet support incremental updating!'
        outflag = 10
     END IF
  END IF

  ! Rank to be considered for inversion of HPH in analysis of EnKF
  IF (dim_pint >= 3) THEN
     IF (param_int(3)>=0 .AND. param_int(3) < dim_ens) THEN
        rank_ana_enkf = param_int(3)
     ELSE
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(103): Invalid setting of param_int(3)/rank_ana_enkf!'
        outflag = 103
        rank_ana_enkf = 0 ! Just for safety: Fall back to default
     END IF
  ELSE
     ! Default mode: Inversion by solving for representers
     rank_ana_enkf = 0
  END IF

  ! Size of lag considered for smoother
  IF (dim_pint>=5) THEN
     dim_lag = 0 !param_int(5)

     IF (param_int(5)>0) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(11): The local EnKF does not yet support smoothing!'
        outflag = 11
     END IF
  END IF

  ! Store type for forgetting factor (SEIK and LSEIK)
  ! We only have a fixed global forgetting factor for EnKF!
  type_forget = 0

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++    Localized Ensemble Kalman Filter (LEnKF)     +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++                                                 +++'     
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Evensen, J. Geophys. Res. 99C (1994) 10143    +++'     
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ using an ensemble of observations according to  +++'     
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ Burgers et al., Mon. Wea. Rev. 126 (1998) 1719  +++'     
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++          This implementation follows            +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++      Nerger et al., Tellus 57A (2005) 715       +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   The localization is covariance lozalization   +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++        of PH^T and HPH^T as described in        +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Houtekamer & Mitchell, MWR, 129 (2001) 123    +++'
     WRITE(*, '(a, 5x, a)') 'PDAF',  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 6x, a)') 'PDAF', 'local EnKF configuration'
     WRITE (*, '(a, 12x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 14x, a)') 'PDAF', '--> local EnKF (analysis for small observation dimension)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 14x, a)') 'PDAF', '--> offline mode'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF', 'PDAF-ERROR(2): No valid sub type!'
        outflag = 2
     END IF
     IF (dim_lag > 0) &
          WRITE (*, '(a, 14x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
     WRITE (*, '(a, 14x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens
     WRITE (*, '(a, 10x, a, f5.2)') 'PDAF', '--> forgetting factor:', forget
     IF (rank_ana_enkf > 0) THEN
        WRITE (*, '(a, 8x, a, i5)') &
             'PDAF', 'analysis with pseudo-inverse of HPH, rank:', rank_ana_enkf
     END IF

  END IF filter_pe2

END SUBROUTINE PDAF_lenkf_init
