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
!$Id: PDAF_enkf_options.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_options --- Information output on options for EnKF
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_options()

! !DESCRIPTION:
! Subroutine to perform information output on options
! available for the EnKF.

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
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++          Ensemble Kalman Filter (EnKF)          +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                 +++'     
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++   Evensen, J. Geophys. Res. 99C (1994) 10143    +++'     
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++ using an ensemble of observations according to  +++'     
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++ Burgers et al., Mon. Wea. Rev. 126 (1998) 1719  +++'     
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++          This implementation follows            +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++      Nerger et al., Tellus 57A (2005) 715       +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for EnKF:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: Full ensemble integration; analysis for 2*dim_obs>dim_ens'
  WRITE(*, '(a, 7x, a)') 'PDAF', '1: Full ensemble integration; analysis for 2*dim_obs<=dim_ens'
  WRITE(*, '(a, 7x, a)') 'PDAF', '5: Offline mode'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): maximum rank for inversion of HPH^T, optional, default=0'
  WRITE(*, '(a, 11x, a)') 'PDAF', '(for =0, HPH is inverted by solving the representer equation)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '(if set to >=ensemble size, it is reset to ensemble size - 1)'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(4): not used'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(5): Size of lag for smoothing'


  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
  WRITE(*, '(a, 11x, a)') &
       'PDAF', '(>=1; not larger than total number of processors)'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
  WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'


  WRITE(*, '(a, 5x, a)') &
       'PDAF', '+++++++++ End of option overview for the EnKF ++++++++++'
  
END SUBROUTINE PDAF_enkf_options
