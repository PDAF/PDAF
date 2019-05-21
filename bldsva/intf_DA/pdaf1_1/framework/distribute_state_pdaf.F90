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
!distribute_state_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                           'distribute_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------


!$Id: distribute_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! During the forecast phase of the filter this
    ! subroutine is called from PDAF\_get\_state
    ! supplying a model state which has to be evolved.
    ! The routine has to initialize the fields of the
    ! model (typically available through a module) from
    ! the state vector of PDAF. With parallelization,
    ! MPI communication might be required to
    ! initialize all subdomains on the model PEs.
    !
    ! The routine is executed by each process that is
    ! participating in the model integrations.
    !
    ! !REVISION HISTORY:
    ! 2004-10 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    use mod_tsmp, &
        only: pf_statevec_fortran, tag_model_parflow, tag_model_clm
    use mod_parallel_model, &
        only: model, mype_world
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
    REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

    ! local variables
    integer :: i
#if defined CLMSA
    real(r8), pointer :: sw_c(:,:)
#endif
    ! !CALLING SEQUENCE:
    ! Called by: PDAF_get_state   (as U_dist_state)
    !EOP


    ! *******************************************
    ! *** Initialize model fields from state  ***
    ! *** Each model PE knows his sub-state   ***
    !********************************************

    !print *, "Distributing state"
    if ((model == tag_model_parflow)) then
        !print *, "Parflow: distrubute_state_pdaf, from state_p to subvec"
        do i=1,dim_p 
          pf_statevec_fortran(i) = state_p(i)
        enddo
    end if

#if defined CLMSA
    !kuw: distribute state vector to clm
    if (model == tag_model_clm) then
        !sw_c  => clm3%g%l%c%cws%h2osoi_vol
        !sw_c = reshape(state_p,shape(sw_c))
        clm_statevec = state_p
    end if
    !kuw end
#endif

END SUBROUTINE distribute_state_pdaf
