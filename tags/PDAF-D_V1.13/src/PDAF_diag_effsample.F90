! Copyright (c) 2014-2016 Paul Kirchgessner, paul.kirchgessner@awi.de
!
! This routine is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! This code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this software.  If not, see <http://www.gnu.org/licenses/>.
!BOP
!
!
! !ROUTINE: PDAF_diag_effsample --- Compute effective sample size
!
! !INTERFACE:
  SUBROUTINE PDAF_diag_effsample(dim_sample, weights, n_eff)

! !DESCRIPTION:
! This routine computes the effective sample size of a particle
! filter as defined in Doucet et al. 2001 p. 333 


! !REVISION HISTORY:
! 2014-06 - Paul Kirchgessner - Initial code for SANGOMA 
! 2016-08 - L. Nerger - adaption for PDAF
!
! !USES:

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: dim_sample         ! Sample size
  REAL, INTENT(in)    :: weights(dim_sample) ! Weights of the samples
  REAL, INTENT(out)   :: n_eff               ! Effecfive sample size
!EOP

! *** local variables ***
  INTEGER :: i     ! Counters

  n_eff = 0
  DO i = 1, dim_sample
     n_eff = n_eff + weights(i)*weights(i)
  ENDDO
  
  IF (n_eff/=0) n_eff = 1.0/n_eff

END SUBROUTINE PDAF_diag_effsample
