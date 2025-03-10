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
!> Abstract interfaces for call-back routines
!!
!! This module provides the abstract interfaces for the 
!! user-supplied call-back routines. Together with using
!! procedure declarations they allow to obtain an
!! explicit interface for the call-back routines
MODULE PDAF_cb_procedures

  IMPLICIT NONE

  ABSTRACT INTERFACE 
     SUBROUTINE init_ens_cb(filtertype, dim_p, dim_ens, &
          state_p, Uinv, ens_p, flag)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: filtertype      !< Type of filter to initialize
       INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens         !< Size of ensemble
       REAL, INTENT(inout) :: state_p(dim_p)              !< PE-local model state
       REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1)   !< Array not referenced for ensemble filters
       REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)       !< PE-local state ensemble
       INTEGER, INTENT(inout) :: flag         !< PDAF status flag
     END SUBROUTINE init_ens_cb
  END INTERFACE


END MODULE PDAF_cb_procedures
