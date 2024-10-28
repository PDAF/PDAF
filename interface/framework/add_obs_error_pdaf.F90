!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!add_obs_error_pdaf.F90: TSMP-PDAF implementation of routine
!                        'add_obs_error_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: add_obs_error_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: add_obs_error_pdaf --- Add observation error covariance matrix
!
! !INTERFACE:
SUBROUTINE add_obs_error_pdaf(step, dim_obs_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: EnKF
!
! The routine is called during the analysis step
! by PDAF\_enkf\_analysis_X (X=rlm or rsm).  It 
! has to add the observation error covariance 
! matrix to the provided matrix C_p for the 
! PE-local domain .
! 
! Implementation for TSMP-PDAF.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs, obs_nc2pdaf

  USE mod_read_obs, ONLY: multierr,clm_obserr, pressure_obserr
  USE mod_parallel_pdaf, ONLY: mype_world
  USE mod_parallel_pdaf, ONLY: abort_parallel
  USE mod_tsmp, ONLY: point_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p  ! Dimension of observation vector
  REAL, INTENT(inout) :: C_p(dim_obs_p,dim_obs_p) ! Matrix to that
                                    ! observation covariance R is added

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm   (as U_add_obs_err)
! Called by: PDAF_enkf_analysis_rsm   (as U_add_obs_err)
!EOP


! *** local variables ***
  INTEGER :: i          ! index of observation component
  REAL :: variance_obs  ! variance of observations


! **********************
! *** INITIALIZATION ***
! **********************

    variance_obs = rms_obs ** 2


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

  if(multierr.ne.1) then
    DO i = 1, dim_obs_p
       C_p(i, i) = C_p(i, i) + variance_obs
    ENDDO
  endif

 
  if(multierr.eq.1) then

    ! Check that point observations are used
    if (.not. point_obs .eq. 1) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR(3) `point_obs.eq.1` needed for using obs_nc2pdaf."
      call abort_parallel()
    end if

    do i=1,dim_obs_p
#if defined CLMSA
      C_p(i,i) = C_p(i,i) + clm_obserr(obs_nc2pdaf(i))*clm_obserr(obs_nc2pdaf(i))
#else
      C_p(i,i) = C_p(i,i) + pressure_obserr(obs_nc2pdaf(i))*pressure_obserr(obs_nc2pdaf(i))
#endif
    enddo
  endif

END SUBROUTINE add_obs_error_pdaf
