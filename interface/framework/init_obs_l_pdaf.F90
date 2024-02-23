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
!init_obs_l_pdaf.F90: TSMP-PDAF implementation of routine
!                     'init_obs_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_obs_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_obs_l_pdaf --- Initialize local observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain in 
! PDAF\_lseik\_analysis.  It has to initialize 
! the local vector of observations for the 
! current local analysis domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code 
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
        ONLY: obs, obs_index_l 

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain index
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_init_obs_l)
! Called by: PDAF_lestkf_analysis  (as U_init_obs_l)
! Called by: PDAF_letkf_analysis   (as U_init_obs_l)
!EOP

! *** local variables ***
  INTEGER :: i       ! counter

! *******************************************
! *** Initialize local observation vector ***
! *******************************************
  
  ! Intialize local obseravtion vector
  observation_l(:) = 0.0

!#ifndef CLMSA
  do i=1,dim_obs_l
    observation_l(i) = obs(obs_index_l(i)) 
  end do

END SUBROUTINE init_obs_l_pdaf

