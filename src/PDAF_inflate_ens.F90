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
!$Id: PDAF_inflate_ens.F90 1629 2016-08-10 12:48:23Z lnerger $
!BOP
!
! !ROUTINE: PDAF_inflate_ens --- Inflate an ensemble 
!
! !INTERFACE:
SUBROUTINE PDAF_inflate_ens(dim, dim_ens, meanstate, ens, forget)

! !DESCRIPTION:
! This routine modifies an input ensemble such that its covariance 
! is inflated by the factor 1/forget.  The ensemble perturbations
! are inflated by 1/sqrt(forget). The ensemble mean is unchanged:
!      Xnew = Xmean + 1/sqrt(forget) * (X - Xmean) 
! The routine returns the inflated ensemble, replacing the input ensemble
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-11 - Julian Toedter - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! dimension of states
  INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
  REAL, INTENT(inout) :: meanstate(dim)    ! state vector to hold ensemble mean
  REAL, INTENT(inout) :: ens(dim, dim_ens) ! Input/output ensemble matrix
  REAL, INTENT(in)    :: forget            ! Forgetting factor

! !CALLING SEQUENCE:
! Called by: PDAF_lnetf_update
! Called by: PDAF_netf_analysis
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  REAL :: invdimens    ! Inverse of ensemble size
  REAL :: infl         ! inflation factor of perturbations = 1/sqrt(forget)


! **********************
! *** INITIALIZATION ***
! **********************

  invdimens = 1.0 / REAL(dim_ens)

  ! Compute ensemble mean state 
  meanstate   = 0.0
  DO col = 1, dim_ens
     DO row = 1, dim
        meanstate(row) = meanstate(row) + ens(row, col)
     END DO
  END DO
  meanstate = invdimens * meanstate


! **********************************************
! ***  Subtract mean and apply inflation     ***
! **********************************************

  ! Get perturbation matrix X'=X-xmean
  DO col = 1, dim_ens
     DO row = 1, dim
        ens(row, col) = ens(row, col) - meanstate(row)
     END DO
  END DO
  
  infl = 1.0/SQRT(forget)

  ! Inflation is done by X'new = infl * X 
  ! Add mean again to get Xnew = X'new + ensmean
  DO col = 1, dim_ens
     DO row = 1, dim
        ens(row, col) = infl * ens(row, col) + meanstate(row)
     END DO
  END DO

END SUBROUTINE PDAF_inflate_ens
