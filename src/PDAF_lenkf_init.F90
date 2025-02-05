! Copyright (c) 2004-2024 Lars Nerger
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
!>  PDAF-internal initialization of ETKF
!!
!! Initialization of ETKF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! *  2015-12 - Lars Nerger - Initial code by copying PDAF_enkf_init
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

  USE PDAF_mod_filter, &
       ONLY: dim_ens, localfilter, dim_lag
  USE PDAF_lenkf, &
       ONLY: incremental, forget, rank_ana_enkf, &
       PDAF_lenkf_set_iparam, PDAF_lenkf_set_rparam
  USE PDAFobs, &
       ONLY: observe_ens

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: subtype               !< Sub-type of filter
  INTEGER, INTENT(in)    :: dim_pint              !< Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint)   !< Integer parameter array
  INTEGER, INTENT(in)    :: dim_preal             !< Number of real parameters 
  REAL, INTENT(inout)    :: param_real(dim_preal) !< Real parameter array
  LOGICAL, INTENT(out)   :: ensemblefilter        !< Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out)   :: fixedbasis            !< Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in)    :: verbose               !< Control screen output
  INTEGER, INTENT(inout) :: outflag               !< Status flag

! *** local variables ***
  INTEGER :: i                ! Counter
  INTEGER :: flagsum          ! Sum of status flags


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Set parameter default values
  ! (Other defaults are set in the module)
  incremental = 0
  observe_ens = .false.
  forget = 1.0
  dim_lag = 0
  rank_ana_enkf = 0

  ! Parse provided parameters
  flagsum = 0
  DO i=3, dim_pint
     CALL PDAF_lenkf_set_iparam(i, param_int(i), outflag)
     flagsum = flagsum+outflag
  END DO
  DO i=1, dim_preal
     CALL PDAF_lenkf_set_rparam(i, param_real(i), outflag)
     flagsum = flagsum+outflag
  END DO


  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Define whether filter is a domain-local filter
  localfilter = 0

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  writeout: IF (verbose == 1) THEN
  
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

     IF (flagsum== 0 ) THEN

        ! *** General output ***
        WRITE (*, '(/a, 6x, a)') 'PDAF', 'local EnKF configuration'
        WRITE (*, '(a, 12x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
        IF (subtype == 0) THEN
           WRITE (*, '(a, 14x, a)') 'PDAF', '--> local EnKF (analysis for small observation dimension)'
        END IF
        WRITE (*, '(a, 10x, a, f5.2)') 'PDAF', '--> forgetting factor:', forget
        IF (rank_ana_enkf > 0) THEN
           WRITE (*, '(a, 8x, a, i5)') &
                'PDAF', 'analysis with pseudo-inverse of HPH, rank:', rank_ana_enkf
        END IF
        WRITE (*, '(a, 14x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens
        IF (observe_ens) &
             WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
     END IF

  END IF writeout

END SUBROUTINE PDAF_lenkf_init
