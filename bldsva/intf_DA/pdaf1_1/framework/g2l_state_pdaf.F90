!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!g2l_state_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                    'g2l_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: g2l_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_model, ONLY: model
  USE mod_tsmp, ONLY: tag_model_parflow, &
       tag_model_clm
  USE mod_tsmp, &
       ONLY: nx_local, ny_local
#ifndef PARFLOW_STAND_ALONE 
  USE decompMod, ONLY: get_proc_bounds_atm
#endif
  USE iso_c_binding, ONLY: c_loc

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, TARGET, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, TARGET, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

  INTEGER :: i, n_domain, nshift_p
  INTEGER :: begg, endg   ! per-proc gridcell ending gridcell indices
! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
!EOP

! *************************************
! *** Initialize local state vector ***
! *************************************
#ifdef CLMSA 
  ! beg and end gridcell for atm
  call get_proc_bounds_atm(begg, endg)
  n_domain = endg - begg + 1
  DO i = 0, dim_l-1
     nshift_p = domain_p + i * n_domain
     state_l(i+1) = state_p(nshift_p)
  ENDDO
#else
#ifdef PARFLOW_STAND_ALONE
     n_domain = nx_local * ny_local
     DO i = 0, dim_l-1
        nshift_p = domain_p + i * n_domain
        state_l(i+1) = state_p(nshift_p)
     ENDDO
#else
  if (model == tag_model_parflow) then
     n_domain = nx_local * ny_local
     DO i = 0, dim_l-1
        nshift_p = domain_p + i * n_domain
        state_l(i+1) = state_p(nshift_p)
     ENDDO
  end if
   
  if (model == tag_model_clm) then
  ! beg and end gridcell for atm
    call get_proc_bounds_atm(begg, endg)
    n_domain = endg - begg + 1
    DO i = 0, dim_l-1
       nshift_p = domain_p + i * n_domain
       state_l(i+1) = state_p(nshift_p)
    ENDDO
  end if
#endif
#endif 


END SUBROUTINE g2l_state_pdaf
