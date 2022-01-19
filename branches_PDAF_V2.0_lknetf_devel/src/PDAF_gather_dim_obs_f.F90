! Copyright (c) 2004-2021 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_gather_dim_obs_f --- get full observation dimension
!
! !INTERFACE:
SUBROUTINE PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

! !DESCRIPTION:
! If the local filter is used with a domain-decomposed model,
! the observational information from different sub-domains
! has to be combined into the full observation vector. 
! This routine is called as a first step to compute the
! full observation dimension from the process-local
! observation dimensions.
! Practically, the operation is a simple MPI_Allreduce. This
! is encapsulated here to simplify the operation for the users. 
! In addition an array storing the pe-local observation dimensions
! and an array of displacements for gathering are initialized.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: COMM_filter, MPIerr, npes_filter, &
       all_dim_obs_p, all_dis_obs_p, dimobs_p, dimobs_f

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in)  :: dim_obs_p    ! PE-local observation dimension
  INTEGER, INTENT(out) :: dim_obs_f    ! Full observation dimension
  
! !CALLING SEQUENCE:
! Called by: user code, usually init_dim_obs_f
! Calls: MPI_Allreduce
! Calls: MPI_Allgather
!EOP

! local variables
  INTEGER :: i  ! Counter


! **********************************************************
! *** Compute global sum of local observation dimensions ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_Allreduce(dim_obs_p, dim_obs_f, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     dim_obs_f = dim_obs_p
  END IF

  ! Store dimensions inside PDAF for use in PDAF_gather_obs_f
  dimobs_p = dim_obs_p
  dimobs_f = dim_obs_f


! ****************************************************************************
! *** Gather and store array of process-local dimensions and displacements ***
! ****************************************************************************

  IF (ALLOCATED(all_dim_obs_p)) DEALLOCATE(all_dim_obs_p)
  ALLOCATE(all_dim_obs_p(npes_filter))
  IF (ALLOCATED(all_dis_obs_p)) DEALLOCATE(all_dis_obs_p)
  ALLOCATE(all_dis_obs_p(npes_filter))

  IF (npes_filter>1) THEN
     CALL MPI_Allgather(dim_obs_p, 1, MPI_INTEGER, all_dim_obs_p, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)

     ! Init array of displacements for observation vector
     all_dis_obs_p(1) = 0
     DO i = 2, npes_filter
        all_dis_obs_p(i) = all_dis_obs_p(i-1) + all_dim_obs_p(i-1)
     END DO
  ELSE
     all_dim_obs_p = dim_obs_p
     all_dis_obs_p = 0
  END IF

END SUBROUTINE PDAF_gather_dim_obs_f
