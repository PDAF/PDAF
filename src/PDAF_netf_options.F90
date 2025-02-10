! Copyright (c) 2004-2025 Lars Nerger
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
!> Information output on options for NETF
!!
!! Subroutine to perform information output on options
!! available for the NETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2016-11 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_lnetf_options()

  IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************
  
  WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(a)')  'PDAF    +++      Nonlinear Ensemble Transform Filter (LNETF)      +++'
  WRITE(*, '(a)')  'PDAF    +++                                                       +++'
  WRITE(*, '(a)')  'PDAF    +++                         by                            +++'
  WRITE(*, '(a)')  'PDAF    +++ J. Toedter, B. Ahrens, Mon. Wea. Rev. 143 (2015) 1347 +++'
  WRITE(*, '(a)')  'PDAF    +++                                                       +++'
  WRITE(*, '(a)')  'PDAF    +++    with local analysis and observation localization   +++'
  WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for LNETF:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: Standard implementation with ensemble integration'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): dim_lag'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Size of smoothing lag (>=0), optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: no smoothing (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '>0: apply smoother up to specified lag'
  WRITE(*, '(a, 7x, a)') 'PDAF', &
       'param_int(4): type_noise'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble perturbations, optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: no perturbations (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: constant standard deviation'
  WRITE(*, '(a, 12x, a)') 'PDAF', '2: relative to ensemble standard deviation'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(5) type_forget'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of forgetting factor; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: forgetting factor on forecast ensemble (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: forgetting factor on forecast ensemble only observed domains'
  WRITE(*, '(a, 12x, a)') 'PDAF', '2: forgetting factor on analysis ensemble'
  WRITE(*, '(a, 12x, a)') 'PDAF', '3: forgetting factor on analysis ensemble only observed domains'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(6) type_trans'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble transformation matrix; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: random orthonormal matrix orthogonal to (1,...,1)^T (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: deterministic transformation'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(7): type_winf'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of weights inflation; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: no weights inflation (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: inflate so that N_eff/N > param_real(2)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(8): observe_ens'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Application of observation operator H, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Note: This parameter has not influence on the LNETF assimilation result'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: Apply H to ensemble mean to compute innovation'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: Apply H to ensemble states; then compute innovation from their mean (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '   param_int(8)=1 is the recomended choice for nonlinear H'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(9): type_obs_init'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Initialize observations before or after call to prepoststep_pdaf'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: Initialize observations before call to prepoststep_pdaf'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: Initialize observations after call to prepoststep_pdaf (default)'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): forget'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(2): limit_winf'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Limit for weigts inflation N_eff/N > param_real(2), optional, default=0.0'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(3): noise_amp'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Ensemble perturbation level (>0), required, only used if param_int(4)>0'

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

  WRITE(*, '(a, 5x, a)') &
       'PDAF', '+++++++++ End of option overview for the LNETF  ++++++++++'
  
END SUBROUTINE PDAF_lnetf_options
