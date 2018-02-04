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
!$Id: PDAF_lseik_options.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lseik_options --- Information output on options for LSEIK
!
! !INTERFACE:
SUBROUTINE PDAF_lseik_options()

! !DESCRIPTION:
! Subroutine to perform information output on options
! available for the LSEIK filter.

! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: PDAF_options_filters
!EOP
  
  WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                  LSEIK Filter                   +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                 +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++        Domain-localized SEIK filter by          +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++   Nerger et al., Ocean Dynamics 56 (2006) 634   +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++      based in the global SEIK filter by         +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++ Pham et al., C. R. Acad. Sci. II, 326(1998) 255 +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++    and Pham, Mon. Wea. Rev. 129 (2001) 1194     +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for LSEIK:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: full ensemble integration; left-sided application of T'
  WRITE(*, '(a, 7x, a)') 'PDAF', '2: Fixed error space basis'
  WRITE(*, '(a, 7x, a)') 'PDAF', '3: Fixed state covariance matrix'
  WRITE(*, '(a, 7x, a)') 'PDAF', '5: Offline mode'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(3): not used'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(4): 1 for incremental updating, 0 else; optional, default: 0'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(5): Type of forgetting factor; optional, default: 0'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: fixed forgetting factor'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: adaptive forgetting factor for full domain (experimental)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: locally adaptive forgetting factor (experimental)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(6): Type of ensemble transformation matrix; optional, default: 0'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: deterministic omega'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: random orthonormal omega orthogonal to (1,...,1)^T'
  WRITE(*, '(a, 11x, a)') &
       'PDAF', '2: use product of 0 with random orthonomal matrix with eigenvector (1,...,1)^T'
  WRITE(*, '(a, 14x, a)') &
       'PDAF', '(experimental; for random transformations, 1 is recommended)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(7): Type of transformation matrix square root; optional, default: 0'
  WRITE(*, '(a, 11x, a)') 'PDAF', '(Only relevant for subtype/=3)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: symmetric square root'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: Cholesky decomposition'


  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
  WRITE(*, '(a, 11x, a)') &
       'PDAF', '>=1 for subtypes 0 and 1; not larger than total number of processors'
  WRITE(*, '(a, 11x, a)') 'PDAF', '=1 required for subtypes 2 and 3'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
  WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Internal parameter (defined inside PDAF) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'Nm1vsN: Normalization of covariance matrix; default: 1'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: normalization with 1/(Ensemble size)'
  WRITE(*, '(a, 14x, a)') 'PDAF', '(original SEIK, mainly for compatibility with older studies)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: normalization with 1/(Ensemble size - 1)'
  WRITE(*, '(a, 14x, a)') 'PDAF', '(sample covariance matrix consistent with other EnKFs)'


  WRITE(*, '(a, 5x, a)') &
       'PDAF', '+++++++++ End of option overview for the LSEIK filter ++++++++++'
  
END SUBROUTINE PDAF_lseik_options
