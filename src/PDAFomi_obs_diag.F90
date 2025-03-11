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
!> PDAF-OMI observation diagnostics
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFomi_obs_diag

  USE PDAFomi_obs_f, &
       ONLY: obs_diag, obs_f_all, n_obstypes, have_obsmean_diag, &
       have_obsens_diag

CONTAINS

!-------------------------------------------------------------------------------
!!> Set observation diagnostics flag
!!
!! This routine sets the flag for observation diagnostics 
!! in PDAF-OMI. One should set this flag before calling
!! PDAFomi_gather_obs.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_obs_diag(diag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: diag          !< Value for observation diagnostics mode

    ! Set value for observation diagnostics type
    obs_diag = diag

  END SUBROUTINE PDAFomi_set_obs_diag



!-------------------------------------------------------------------------------
!!> Compute RMS deviation beetween observation and observed ensemble mean
!!
!! This routine computes the RMSD between the observation and the
!! observed ensemble mean.
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_obs_rmsd(nobs, rmsd)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: nobs          !< Number of observation types
    REAL, INTENT(inout) :: rmsd(:)          !< Vector of RMSD values

! *** Local variables ***
    INTEGER :: i         ! Counter


    ! Set number of obstypes
    nobs = n_obstypes

    rmsd = 0.0

    IF (have_obsmean_diag>0 .OR. have_obsens_diag>0) THEN

       DO i = 1, n_obstypes

!          IF (have_obsmean_diag>0) THEN
             write (*,*) 'OBStype', i, obs_f_all(i)%ptr%obs_diag_p
!          END IF

       END DO

    END IF

  END SUBROUTINE PDAFomi_get_obs_rmsd

END MODULE PDAFomi_obs_diag
