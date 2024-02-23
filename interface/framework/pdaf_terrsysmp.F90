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
!pdaf_terrsysmp.F90: Driver program for TSMP-PDAF
!-------------------------------------------------------------------------------------------

!> @author Wolfgang Kurtz, Guowei He
!> @date 28.02.2022
!> @brief Main program for TSMP-PDAF
!> @details
!> Main TSMP-PDAF program.
PROGRAM pdaf_terrsysmp

! !USES:
    USE mpi, &
      ONLY: MPI_INIT, MPI_FINALIZE
    !, MPI_BARRIER, MPI_COMM_WORLD, MPI_SUCCESS
  
    USE mod_parallel_pdaf, &
        ONLY : mype_world, MPIerr

    USE mod_tsmp, &
      ONLY: initialize_tsmp, integrate_tsmp, update_tsmp, finalize_tsmp, &
      tcycle, total_steps

    USE mod_assimilation, &
      ONLY: screen

    IMPLICIT NONE

! !CALLING SEQUENCE:
! Calls: MPI_INIT
! Calls: init_parallel_pdaf
! Calls: initialize_tsmp
! Calls: init_pdaf
! Calls: integrate_tsmp
! Calls: assimilate_pdaf
! Calls: update_tsmp
! Calls: finalize_tsmp
! Calls: MPI_FINALIZE

    ! initialize mpi
    CALL MPI_INIT(MPIerr)

    ! intitialize parallel pdaf (communicators et al.)
    CALL init_parallel_pdaf(0, 1)

    ! Read TSMP-PDAF input from "enkfpf.par"
    ! initialize TSMP instances
    CALL initialize_tsmp()

    ! initialize pdaf variables
    CALL init_pdaf()

    ! barrier before model integration starts
    ! call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    ! if (MPIerr .ne. MPI_SUCCESS) then
    !     print *, "barrier failed"
    ! end if

    ! time loop
    DO tcycle = 1, total_steps
        IF (mype_world > -1 .AND. screen > 2) THEN
            PRINT *, "TSMP-PDAF mype(w)=", mype_world, ": time loop", tcycle
        ENDIF

        ! forward simulation of component models
        CALL integrate_tsmp()

        ! assimilation step
        CALL assimilate_pdaf()

        !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
        !print *,"Finished assimilation", tcycle

        !call print_update_pfb()
        CALL update_tsmp()

        !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
        !print *,"Finished complete assimilation cycle", tcycle
    ENDDO

    ! barrier after model integrations
    !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    !print *, "pdaf: advancing finished, rank ", mype_world

    ! close pdaf
    CALL finalize_pdaf()
    !print *, "pdaf: finalized, rank ", mype_world

    ! close model instances
    CALL finalize_tsmp()
    !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    !print *, "model: finalized, rank ", mype_world

    ! close mpi
    CALL MPI_FINALIZE(MPIerr)

END PROGRAM pdaf_terrsysmp
