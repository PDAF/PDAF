! Copyright (c) 2004-2016 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF-D_enkf_gather_resid.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_gather_resid --- perform gathering of local residuals
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)

! !DESCRIPTION:
! This routine performs an allgather operation during 
! the analysis step of the domain-decomposed EnKF. 
! This operation is separated into a subroutine 
! for compactness of the analysis routine.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: npes_filter, MPIerr, COMM_filter, MPI_INTEGER, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_obs   ! Global observation dimension
  INTEGER, INTENT(in) :: dim_obs_p ! PE-local observation dimension
  INTEGER, INTENT(in) :: dim_ens   ! Ensemble size
  REAL, INTENT(out) :: resid(dim_obs, dim_ens)    ! Global residual matrix
  REAL, INTENT(in) :: resid_p(dim_obs_p, dim_ens) ! PE-local residual matrix

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm
! Called by: PDAF_enkf_analysis_rsm
! Calls: PDAF_memcount
! Calls: MPI_allgather (MPI)
! Calls: MPI_AllgatherV (MPI)
!EOP

! *** local variables ***
  INTEGER :: i                     ! Counter   
  INTEGER, SAVE :: allocflag = 0   ! Flag for first-time allocation
  INTEGER, ALLOCATABLE :: local_dim_obs(:) ! Array of PE-local observation dimensions
  INTEGER, ALLOCATABLE :: local_dis(:)     ! Array of PE-local displacements


! *******************************************************
! *** Allgather field of local observation dimensions ***
! *******************************************************
  ALLOCATE(local_dim_obs(npes_filter))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', npes_filter)

  IF (npes_filter>1) THEN
     CALL MPI_allgather(dim_obs_p, 1, MPI_INTEGER, local_dim_obs, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)
  ELSE
     ! This is a work around for working with nullmpi.F90
     local_dim_obs = dim_obs_p
  END IF


! *****************************************************
! *** Allgather residual matrix                     ***
! *** We use simple ordering here!!                 ***
! *** The chosed ordering is unimportans as long as ***
! *** the ordering for computing HPH and HP is the  ***
! *** same as for the residuals                     ***
! *****************************************************
  ALLOCATE(local_dis(npes_filter))
  IF (allocflag == 0) THEN
     CALL PDAF_memcount(3, 'i', npes_filter)
     allocflag = 1
  END IF

  ! Init array of displacements
  local_dis(1) = 0
  DO i = 2, npes_filter
     local_dis(i) = local_dis(i - 1) + local_dim_obs(i - 1)
  END DO

  DO i = 1, dim_ens
     CALL MPI_AllGatherV(resid_p(1 : dim_obs_p, i), dim_obs_p, MPI_REALTYPE, &
          resid(1 : dim_obs, i), local_dim_obs, local_dis, MPI_REALTYPE, &
          COMM_filter, MPIerr)
  END DO

! *** Clean up ***
  DEALLOCATE(local_dis, local_dim_obs)

END SUBROUTINE PDAF_enkf_gather_resid
