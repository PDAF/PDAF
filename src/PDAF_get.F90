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
!
!> Module providing PDAF_get routines
MODULE PDAF_get

CONTAINS

!> Query whether assimilation was performed at current time step
!!
!! Helper routine for PDAF.
!! The routine allows to query whether observations were assimilated 
!! at the most recent to to PDAF_assimilate_X.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2018-03 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_get_assim_flag(did_assim)

  USE PDAF_mod_filter, &
       ONLY: assim_flag

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(out) :: did_assim    !< Flag: (1) for assimilation; (0) else


! *** Set ensemble member ***

  did_assim = assim_flag

END SUBROUTINE PDAF_get_assim_flag

!--------------------------------------------------------------------------
!> Set pointer to ensemble statistics
!!
!! Routine to set the pointer to the PDAF-internal array of skewness and kurtosis.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2020-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_ensstats(skew_ptr, kurt_ptr, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: skewness, kurtosis

  IMPLICIT NONE

! *** Arguments ***
  REAL, POINTER, INTENT(out) :: skew_ptr(:)  !< Pointer to skewness array
  REAL, POINTER, INTENT(out) :: kurt_ptr(:)  !< Pointer to kurtosis array
  INTEGER, INTENT(out)       :: status  !< Status flag 

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(skewness)) THEN
     skew_ptr => skewness
     kurt_ptr => kurtosis

     status = 0
  END IF

END SUBROUTINE PDAF_get_ensstats


!--------------------------------------------------------------------------
!> Query whether chosen filter is domain-localized
!!
!! Routine to return the information whether the current filter
!! is domain-localized. The valu eof localfilter is set in
!! the initialization routine of a filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_localfilter(localfilter_out)

  USE PDAF_mod_filter, &
       ONLY: localfilter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: localfilter_out  !< Whether the filter is domain-localized

  
! ***********************
! *** Set localfilter ***
! ***********************

  localfilter_out = localfilter
  
END SUBROUTINE PDAF_get_localfilter


!--------------------------------------------------------------------------
!> Query whether chosen filter uses localization
!!
!! Routine to return the information whether the current filter
!! uses domain or covariance localization. The values of
!! localfilter and covarloc are set the initialization routine
!! of a filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_local_type(localtype)

  USE PDAF_mod_filter, &
       ONLY: localfilter, covarloc

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: localtype !< Localization type of the filter
                                    !< * (0) no localization; global filter
                                    !< * (1) domain localization (LESTKF, LETKF, LNETF, LSEIK)
                                    !< * (2) covariance localization (LEnKF)
                                    !< * (3) covariance loc. but observation handling like domain localization (ENSRF)

  
! *********************
! *** Set localtype ***
! *********************

  localtype = 0
  IF (localfilter == 1 .AND. covarloc == 0) THEN
     localtype = 1
  ELSEIF (localfilter == 0 .AND. covarloc == 1) THEN
     localtype = 2
  ELSEIF (localfilter == 1 .AND. covarloc == 1) THEN
     localtype = 3
  END IF

END SUBROUTINE PDAF_get_local_type


!--------------------------------------------------------------------------
!> Query ensemble index of the current member
!!
!! Helper routine for PDAF.
!! The routine allows to query the member index of the ensemble
!! state that is currently integrated.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_memberid(memberid)

  USE PDAF_mod_filter, &
       ONLY: member_save

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(inout) :: memberid    !< Index in the local ensemble

! *** Set ensemble member ***

  memberid = member_save

END SUBROUTINE PDAF_get_memberid


!--------------------------------------------------------------------------
!> Query ensemble index of the member calling U_obs_op
!!
!! Helper routine for PDAF.
!! The routine allows to query the member index of the ensemble
!! state that is currently integrated.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_obsmemberid(memberid)

  USE PDAF_mod_filter, &
       ONLY: obs_member

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(inout) :: memberid    !< Index in the local ensemble


! *** Set ensemble member ***

  memberid = obs_member

END SUBROUTINE PDAF_get_obsmemberid


!--------------------------------------------------------------------------
!> Set pointer to smoother ensemble
!!
!! Routine to set the pointer to the PDAF-internal smoother ensemble array.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-05 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_get_smootherens(sens_point, maxlag, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: sens, cnt_maxlag, dim_lag

  IMPLICIT NONE

! *** Arguments ***
  REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  !< Pointer to smoother array
  INTEGER, INTENT(out)       :: maxlag  !< Number of past timesteps processed in sens
  INTEGER, INTENT(out)       :: status  !< Status flag 

  
! *******************
! *** Set pointer ***
! *******************

  status = 1

  IF (allocated(sens)) THEN
     sens_point => sens

     status = 0
  END IF

  ! Set number of initialized lags
  IF (cnt_maxlag-1 >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     maxlag = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     maxlag = cnt_maxlag-1
  END IF

END SUBROUTINE PDAF_get_smootherens

END MODULE PDAF_get
