! Copyright (c) 2014-2018 Paul Kirchgessner
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
! !ROUTINE: PDAF_etkf_options --- Information output on options for ETKF
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_options()

! !DESCRIPTION:
! Subroutine to perform information output on options
! available for the Ensemble Transform Kalman filter.

! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: PDAF_options_filters
!EOP
  
  WRITE(*, '(/8x, a)') '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  WRITE(*, '(8x, a)')  '+++       EQUVIVALENT WEIGHTS PARTICLE FILTER       +++'
  WRITE(*, '(8x, a)')  '+++                                                 +++'
  WRITE(*, '(8x, a)')  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  WRITE(*, '(/8x, a)') 'Available options:'

  WRITE(*, '(/8x, a)') 'Sub-types (Parameter subtype)'
  WRITE(*, '(10x, a)') '0: Equal weights particle filter with nudging'
  
  WRITE(*, '(/8x, a)') 'Integer parameters (Array param_int)'
  WRITE(*, '(10x, a)') 'param_int(1): Dimension of state vector (>0), required'
  WRITE(*, '(10x, a)') 'param_int(2): Ensemble size (>0), required'
  WRITE(*, '(10x, a)') 'param_int(3): Type of resampling'
  WRITE(*, '(/8x, a)') 'Further parameters'
  WRITE(*, '(10x, a)') 'n_modeltasks: Number of parallel model integration tasks'
  WRITE(*, '(14x, a)') &
       '>=1; not larger than total number of processors'
  WRITE(*, '(10x, a)') 'screen: Control verbosity of PDAF'
  WRITE(*, '(14x, a)') '0: no outputs'
  WRITE(*, '(14x, a)') '1: basic output (default)'
  WRITE(*, '(14x, a)') '2: 1 plus timing output'
  WRITE(*, '(14x, a)') '3: 2 plus debug output'

  WRITE(*, '(/8x, a)') &
       '+++++++++ End of option overview for the EWPF ++++++++++'
  
END SUBROUTINE PDAF_ewpf_options
