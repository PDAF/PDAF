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
!l2g_state_pdaf.F90: TSMP-PDAF implementation of routine
!                     'l2g_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: l2g_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: l2g_state_pdaf --- Initialize full state from local analysis
!
! !INTERFACE:
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_X\_update 
! after the analysis and ensemble transformation 
! on a single local analysis domain. It has to 
! initialize elements of the PE-local full state 
! vector from the provided analysis state vector 
! on the local analysis domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_tsmp, ONLY: tag_model_parflow, &
       tag_model_clm, model
  USE mod_tsmp, &
       ONLY: nx_local, ny_local
  USE iso_c_binding, ONLY: c_loc

#if defined CLMSA
  use enkf_clm_mod, ONLY: l2g_state_clm
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  REAL, TARGET, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL, TARGET, INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

  INTEGER :: i, n_domain, nshift_p
  INTEGER :: begg, endg   ! per-proc gridcell ending gridcell indices
! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
! Called by: PDAF_lestkf_update   (as U_l2g_state)
! Called by: PDAF_letkf_update    (as U_l2g_state)
!EOP

! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************
#ifndef CLMSA
  if (model == tag_model_parflow) then
     n_domain = nx_local * ny_local
     DO i = 0, dim_l-1
        nshift_p = domain_p + i * n_domain
        state_p(nshift_p) = state_l(i+1)
     ENDDO
  else  if (model == tag_model_clm) then
     state_p(domain_p) = state_l(dim_l) 
  end if   
  !call l2g_state(domain_p, c_loc(state_p), dim_l, c_loc(state_l))
#else
  call l2g_state_clm(domain_p, dim_l, state_l, dim_p, state_p)
#endif
  
END SUBROUTINE l2g_state_pdaf
