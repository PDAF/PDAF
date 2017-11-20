!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
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
!collect_state_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                        'collect_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: collect_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X after the 
! propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    use mod_tsmp, &
        only: pf_statevec_fortran, tag_model_parflow, tag_model_clm
    use mod_parallel_model, &
        only: model
#if defined CLMSA
    !kuw: get access to clm variables
    USE clmtype      , only : clm3
    USE clm_varpar   , only : nlevsoi
    use shr_kind_mod, only: r8 => shr_kind_r8
    use enkf_clm_mod, only: clm_statevec
    !kuw end
#endif

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

  ! local variables
  integer i
#if defined CLMSA
  real(r8), pointer :: sw_c(:,:)
#endif
! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

 if (model == tag_model_parflow) then
     !print *, "Parflow: collect_state_pdaf, from subvecs to state_p"
     do i=1,dim_p
       state_p(i) = pf_statevec_fortran(i)
     enddo
 end if

#if defined CLMSA
 !kuw: define state vector for clm
 if (model == tag_model_clm) then
     !sw_c  => clm3%g%l%c%cws%h2osoi_vol
     !state_p = reshape(sw_c,shape(state_p))
     state_p = clm_statevec
 end if
 !kuw end
#endif
  
END SUBROUTINE collect_state_pdaf
