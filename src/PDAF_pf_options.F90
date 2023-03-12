! Copyright (c) 2004-2021 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_pf_options --- Information output on options for PF
!
! !INTERFACE:
SUBROUTINE PDAF_pf_options()

! !DESCRIPTION:
! Subroutine to perform information output on options
! available for the PF.

! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: PDAF_options_filters
!EOP
  
  WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(a)')  'PDAF    +++           Particle Filter with resampling             +++'
  WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for PF:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: Standard implementation with ensemble integration'
  WRITE(*, '(a, 7x, a)') 'PDAF', '5: Offline mode'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): Resampling type, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: probabilistic resamping (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: stochastic universal resampling'
  WRITE(*, '(a, 11x, a)') 'PDAF', '3: residual resampling'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(4): Type of ensemble perturbations, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no perturbations (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: constant standard deviation'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: relative to ensemble standard deviation'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(5): Type of forgetting factor; optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: forgetting factor on forecast ensemble (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: forgetting factor on analysis ensemble'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(6): Type of weights inflation; optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no weights inflation (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: inflate so that N_eff/N > param_real(2)'


  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): Ensemble pert. level (>0), required, only used if param_int(4)>0'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(2): Forgetting factor (usually >0 and <=1), optional, default=1.0'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(3): Limit for weigts inflation N_eff/N > param_real(2), optional, default=0.0'

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
       'PDAF', '+++++++++ End of option overview for the PF  ++++++++++'
  
END SUBROUTINE PDAF_pf_options
