!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!init_n_domains_pdaf.F90: TSMP-PDAF implementation of routine
!                        'init_n_domains_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_n_domains_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
!
! !INTERFACE:
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to set the number of local analysis 
! domains for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_tsmp, ONLY: tag_model_parflow, &
      tag_model_clm, model
  USE mod_assimilation, &
      ONLY: dim_state_p
  USE mod_tsmp, &
      ONLY: init_n_domains_size
#if defined CLMSA
#if defined CLMFIVE
  USE decompMod, ONLY: get_proc_bounds
#else
  USE decompMod, ONLY: get_proc_bounds_atm
#endif
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains
#if defined CLMSA
  INTEGER :: begg, endg   ! per-proc gridcell ending gridcell indices
#endif

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_n_domains)
! Called by: PDAF_lestkf_update  (as U_init_n_domains)
! Called by: PDAF_letkf_update   (as U_init_n_domains)
!EOP

! ************************************
! *** Initialize number of domains ***
! ************************************
#if (defined PARFLOW_STAND_ALONE || defined COUP_OAS_PFL)
  if (model.eq.tag_model_parflow) then
     ! Here simply the process-local state dimension
     call init_n_domains_size(n_domains_p)
  end if
#endif   

#ifndef CLMSA
  if (model == tag_model_clm) then
     ! Here simply the process-local state dimension  
     n_domains_p = dim_state_p
  end if   
#else
  ! beg and end gridcell for atm
#if defined CLMFIVE
  call get_proc_bounds(begg, endg)
#else  
  call get_proc_bounds_atm(begg, endg)
#endif
  ! Here simply the process-local state dimension  
  n_domains_p = endg - begg + 1
#endif

END SUBROUTINE init_n_domains_pdaf
