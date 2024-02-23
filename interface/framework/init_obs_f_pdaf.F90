!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and M.Pondkule (Forschungszentrum Juelich GmbH)
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
!init_obs_f_pdaf.F90: TSMP-PDAF implementation of routine
!                     'init_obs_f_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_obs_f_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_obs_f_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_lseik\_update
! before the loop over all local analysis domains
! is entered. It has to provide the full observation 
! vector according to current time step (where 'full' 
! means 'all observations required for the localized 
! analysis on the PE-local domain).  This routine 
! is only used for LSEIK if a globally adaptive 
! forgetting factor is requested, rather than an 
! individual forgetting factor for each analysis 
! domain. This routine has to be implemented 
! consistently with the routines for the full 
! observation dimension and the full observation 
! operator. The forgetting factor will only be 
! globally adaptive, if the full observation vector 
! is the global observation vector.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code 
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_obs)
! Called by: PDAF_lestkf_update  (as U_init_obs)
! Called by: PDAF_letkf_update   (as U_init_obs)
!EOP


! ******************************************
! *** Initialize full observation vector ***
! ******************************************
  
  observation_f = obs

END SUBROUTINE init_obs_f_pdaf

