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
!
!> Information output on options for GLOBALTEMPLATE
!!
!! Subroutine to perform information output on options
!! available for your GLOBALTEMPLATE DA method.
!!
!! This routine is called by PDAF_options_filters.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_GLOBALTEMPLATE_options()

  IMPLICIT NONE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                          +++
! +++ Upon implementation of a DA method one should adapt this routine   +++
! +++ to provide a complete overview of the available options            +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(a, 4x, a)')  'PDAF', '+++      TEMPLATE FOR GLOBALTEMPLATE DA METHOD      +++'
  WRITE(*, '(a, 4x, a)')  'PDAF', '+++ replace this string upon implemeting the method +++'
  WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for GLOBALTEMPLATE:'

! +++ Describe available subtypes (at least subtype=0 is required)
  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
  WRITE(*, '(a, 7x, a)') 'PDAF', '0: describe subtype 0'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
! +++ param_int(1) and param_int(2) are mandatory and should not be changed
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
! +++ additional method-specific integer parameters
  WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): Size of smoothing lag (>=0), optional'
  WRITE(*, '(a, 11x, a)') 'PDAF', '0: no smoothing (default)'
  WRITE(*, '(a, 11x, a)') 'PDAF', '>0: apply smoother up to specified lag'

  WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
! +++ param_real(1) is mandatory and should not be changed
  WRITE(*, '(a, 7x, a)') &
       'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required'
! +++ here one can add the method-specific real parameters

! +++ Further generic output that usually does not need changes
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
       'PDAF', '+++++++++ End of option overview for the GLOBALTEMPLATE ++++++++++'
  
END SUBROUTINE PDAF_GLOBALTEMPLATE_options
