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
!> Information output on options for LNETF
!!
!! Subroutine to perform information output on options
!! available for the LNETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2018-07 - Lars Nerger - Initial code based on code for LNETF
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_lknetf_options()

  IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************
  
  WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++  Local Hybrid Kalman-Nonlinear Ensemble Transform Filter  +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                           +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                Domain-localized LKNETF by                 +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++ L. Nerger, QJRMS, 148 (2022) 620-640, doi:10.1002/qj.4221 +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for LKNETF:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: HNK: 2-step LKNETF with NETF before LETKF'
  WRITE(*, '(a, 7x, a)') 'PDAF', '1: HKN: 2-step LKNETF with LETKF before NETF'
  WRITE(*, '(a, 7x, a)') 'PDAF', '4: HSync: LKNETF synchronous'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): not used'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(4): not used'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(5): type_forget'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of forgetting factor; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: inflate forecast ensemble by 1/forget (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: inflate forecast ensemble by 1/forget only observed domains'
  WRITE(*, '(a, 12x, a)') 'PDAF', '2: inflate analysis ensemble by 1/forget'
  WRITE(*, '(a, 12x, a)') 'PDAF', '3: inflate analysis ensemble by 1/forget only observed domains'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(6): type_trans'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble transformation matrix; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: random orthonormal matrix orthogonal to (1,...,1)^T (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: deterministic transformation'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(7): type_hyb'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of hybrid weight; optional'
  WRITE(*, '(a, 12x, a)') 'PDAF', '0: fixed value'
  WRITE(*, '(a, 12x, a)') 'PDAF', '1: gamma_lin: (1 - N_eff/N_e)*param_real(2) (default)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '2: gamma_alpha: hybrid weight from N_eff/N>=param_real(2)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '3: gamma_ska: 1 - min(s,k)/sqrt(param_real(3)) with N_eff/N>=param_real(2)'
  WRITE(*, '(a, 12x, a)') 'PDAF', '4: gamma_sklin: 1 - min(s,k)/sqrt(param_real(3)) >= 1-N_eff/N>=param_real(2)'


  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): forget'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(2): hyb_g'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'prescribed hybrid weight gamma (usually >0 and <=1), optional, default=0.95'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(3): hyb_k'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'hybrid norm kappa (>0), optional, default=dim_ens'

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
       'PDAF', '+++++++++ End of option overview for the LKNETF  ++++++++++'
  
END SUBROUTINE PDAF_lknetf_options
