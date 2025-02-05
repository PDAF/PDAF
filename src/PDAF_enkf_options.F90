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
!> Information output on options for EnKF
!!
!! Subroutine to perform information output on options
!! available for the EnKF filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2011-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_enkf_options()

  IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************
  
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
  WRITE(*, '(a, 7x, a)') 'PDAF', '5: Offline mode (deprecated, use PDAF_set_offline_mode)'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): rank_ana_enkf'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'maximum rank for inversion of HPH^T, optional, default=0'
  WRITE(*, '(a, 12x, a)') 'PDAF', 'for =0, HPH is inverted by solving the representer equation'
  WRITE(*, '(a, 12x, a)') 'PDAF', 'allowed range is 0 to ensemble size - 1'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(4): not used'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(5): dim_lag'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Size of smoothing lag (>=0), optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no smoothing (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '>0: apply smoother up to specified lag'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(6): not used'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(7): not used'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(8): observe_ens'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Application of observation operator H, optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: Apply H to ensemble mean to compute innovation'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: Apply H to ensemble states; then compute innovation from their mean (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '   param_int(8)=1 is the recomended choice for nonlinear H'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(9): type_obs_init'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Initialize observations before or after call to prepoststep_pdaf'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: Initialize observations before call to prepoststep_pdaf'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: Initialize observations after call to prepoststep_pdaf (default)'



  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): forget'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'

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
