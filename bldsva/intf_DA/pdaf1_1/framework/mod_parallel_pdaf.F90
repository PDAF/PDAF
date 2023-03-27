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
!mod_parallel_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                       'mod_parallel_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: mod_parallel_pdaf.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_parallel_pdaf

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! to be shared between model-related routines. The are variables
! that are used in the model, even without PDAF and additional
! variables that are only used, if data assimialtion with PDAF
! is performed.
! In addition methods to initialize and finalize MPI are provided.
! The initialization routine is only for the model itself, the 
! more complex initialization of communicators for xecution with
! PDAF is peformed in init\_parallel\_pdaf.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use iso_c_binding, only: c_int

  IMPLICIT NONE
  SAVE 

  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  ! Additional variables for use with PDAF
  INTEGER(c_int), bind(c) :: n_modeltasks         ! Number of parallel model tasks
  INTEGER :: n_filterpes  = 1         ! Number of PEs for filter analysis
  INTEGER :: COMM_filter ! MPI communicator for filter PEs 
  INTEGER(c_int), bind(c) :: mype_filter ! PE rank in COMM_filter
  INTEGER(c_int), bind(c) :: npes_filter ! # PEs in COMM_filter
  INTEGER :: COMM_couple ! MPI communicator for coupling filter and model
  LOGICAL :: modelpe     ! Whether we are on a PE in a COMM_model
  LOGICAL :: filterpe    ! Whether we are on a PE in a COMM_filter
  INTEGER :: task_id     ! Index of my model task (1,...,n_modeltasks)
  INTEGER :: MPIerr      ! Error flag for MPI
  INTEGER :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  INTEGER, ALLOCATABLE :: local_npes_model(:) ! # PEs per ensemble
!  TODO: for statistics
  bind(c) :: COMM_couple
  bind(c) :: task_id
!EOP

END MODULE mod_parallel_pdaf
