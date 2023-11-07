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
!init_ens.F90: TerrSysMP-PDAF implementation of routine
!              'init_ens' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_ens.F90 1444 2013-10-04 10:54:08Z lnerger $
!BOP
!
! !ROUTINE: init_ens --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_ens(filtertype, dim_p, dim_ens, state_p, Uinv, &
    ens_p, flag)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! The routine is called when the filter is
    ! initialized in PDAF\_filter\_init.  It has
    ! to initialize an ensemble of dim\_ens states.
    ! Typically, the ensemble will be directly read from files.
    !
    ! The routine is called by all filter processes and
    ! initializes the ensemble for the PE-local domain.
    !
    ! !REVISION HISTORY:
    ! 2013-02 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    !   USE mod_model, &
    !        ONLY: nx, ny, nx_p
    USE mod_parallel_model, &
        ONLY: mype_model, model, &
        mype_world
    USE mod_parallel_pdaf, &
        ONLY: mype_filter, task_id
    use mod_tsmp, &
        only: tag_model_parflow, pf_statevecsize, pf_statevec, pf_statevec_fortran
    use mod_assimilation, &
        only: screen

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
    ! It is not necessary to initialize the array 'state_p' for SEIK.
    ! It is available here only for convenience and can be used freely.
    REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
    REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
    INTEGER, INTENT(inout) :: flag                 ! PDAF status flag
    ! !CALLING SEQUENCE:
    ! Called by: PDAF_filter_init    (as U_ens_init)
    !EOP

    ! *** local variables ***
    INTEGER :: i, j, member  ! Counters


    ! **********************
    ! *** INITIALIZATION ***
    ! **********************

    ! *** Generate full ensemble on filter-PE 0 ***
    IF (mype_filter==0 .and. screen > 0) THEN
        WRITE (*, '(/9x, a)') 'Initialize state ensemble'
        WRITE (*, '(9x, a)') '--- read ensemble from files'
        WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
    END IF



    ! ********************************
    ! *** Read ensemble from files ***
    ! ********************************

    !    WRITE (*,*) 'TEMPLATE init_ens.F90: Initialize ensemble array ens_p!'

    !    convert pf_statevec to fortran pointer

    if (model == tag_model_parflow) then


        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": Parflow component: initialize ensemble array ens_p"
            print *, "TSMP-PDAF mype(w)=", mype_world, ": my dim_p is", dim_p
        end if

        do i = 1, dim_ens
            ens_p(:, i) = 10 + mype_model + i
        end do
    else
        do i = 1, dim_ens
            ens_p(:, i) = 1.1
        end do

        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, " CLM component: initialize ensemble array ens_p"
            print *, "TSMP-PDAF mype(w)=", mype_world, " my dim_p is", dim_p
        end if

    end if
! ****************
! *** clean up ***
! ****************


END SUBROUTINE init_ens
