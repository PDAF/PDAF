! Copyright (c) 2004-2023 Lars Nerger
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
! !ROUTINE: PDAF_3dvar_options --- Information output on options for 3DVAR
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_options()

! !DESCRIPTION:
! Subroutine to perform information output on options
! available for the Ensemble Transform Kalman filter.

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
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                      3D-Var                     +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                 +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++      3D-Var variants implemented following      +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++      Bannister, Q. J. Royal Meteorol. Soc.,     +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++     143 (2017) 607-633, doi:10.1002/qj.2982     +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for 3D-Var:'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: incremental 3D-Var with parameterized covariance matrix'
  WRITE(*, '(a, 7x, a)') 'PDAF', '1: 3D ensemble Var using LESTKF for ensemble transformation'
  WRITE(*, '(a, 7x, a)') 'PDAF', '4: 3D ensemble Var using ESTKF for ensemble transformation'
  WRITE(*, '(a, 7x, a)') 'PDAF', '5: Offline mode  (deprecated, use PDAF_set_offline_mode)'
  WRITE(*, '(a, 7x, a)') 'PDAF', '6: hybrid 3D-Var using LESTKF for ensemble transformation'
  WRITE(*, '(a, 7x, a)') 'PDAF', '7: hybrid 3D-Var using ESTKF for ensemble transformation'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(3): Select optimization method (solver), required'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: LBFGS (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '1: CG+'
  WRITE(*, '(a, 11x, a)') 'PDAF', '2: direct implementation of CG'
  WRITE(*, '(a, 11x, a)') 'PDAF', '3: direct implementation of CG with decomposed control vector'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(4): size of parameterized control vector (for parameterized and hybrid 3D-Var), required'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(5): size of ensemble control vector (required for ensemble and hybrid 3D-Var), '
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(6): Solver-specific parameter, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: parameter m (default=5)'
  WRITE(*, '(a, 16x, a)') 'PDAF', 'Number of corrections used in limited memory matrix; 3<=m<=20'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: parameter method (default=2)'
  WRITE(*, '(a, 16x, a)') 'PDAF', '(1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'CG: maximum number of iterations (default=200)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_int(7): Solver-specific parameter, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: parameter irest (default=1)'
  WRITE(*, '(a, 16x, a)') 'PDAF', '(0) no restarts; (n>0) restart every n steps'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required;'
  WRITE(*, '(a, 11x, a)') 'PDAF', '(only used for ensemble and hybrid 3D-Var)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(2): hybrid weight beta, optional (only for hybrid 3D-Var)'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'range >=0.0 and <=1.0, =1.0 for pure ensemble 3D-var  (default=0.5)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(3): Solver-specific parameter, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: Limit for stopping iterations (pgtol, default=1.0e-5)'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: convergence parameter eps (default=1.0e-5)'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'CG: convergence parameter eps (default=1.0e-6)'
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(4): Solver-specific parameter, optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: Tolerance in termination test (factr, default=1.0e+7)'

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
       'PDAF', '+++++++++ End of option overview for 3DVAR ++++++++++'
  
END SUBROUTINE PDAF_3dvar_options
