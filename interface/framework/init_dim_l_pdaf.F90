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
!init_dim_l_pdaf.F90: TSMP-PDAF implementation of routine
!                     'init_dim_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the loop over all local analysis domain.
! It has to set the dimension of local model 
! state on the current analysis domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_tsmp, ONLY: tag_model_parflow, &
       tag_model_clm, model
  USE mod_tsmp, &
       ONLY: init_parf_l_size
#ifdef CLMSA
  USE enkf_clm_mod, &
       ONLY: init_clm_l_size
#endif
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
!EOP

! ****************************************
! *** Initialize local state dimension ***
! ****************************************
#if (defined PARFLOW_STAND_ALONE || defined COUP_OAS_PFL)
  if (model.eq.tag_model_parflow) then
     ! Set the size of the local analysis domain 
     call init_parf_l_size(dim_l)
  end if
#endif  

#ifndef CLMSA
  if (model.eq.tag_model_clm) then
     ! Set the size of the local analysis domain   
     dim_l = 1     
  end if   
#else
  ! Set the size of the local analysis domain  
  ! for clm stand alone mode only
  call init_clm_l_size(dim_l)
#endif
 
END SUBROUTINE init_dim_l_pdaf
