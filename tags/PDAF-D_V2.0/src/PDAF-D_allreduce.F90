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
! !ROUTINE: PDAF_allreduce --- Perform an allreduce in COMM_filter
!
! !INTERFACE:
SUBROUTINE PDAF_allreduce(val_p, val_g, mpitype, mpiop, status)

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
! 2019-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: COMM_filter, MPIerr, npes_filter

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in)  :: val_p    ! PE-local value
  INTEGER, INTENT(out) :: val_g    ! reduced global value
  INTEGER, INTENT(in)  :: mpitype  ! MPI data type
  INTEGER, INTENT(in)  :: mpiop    ! MPI operator
  INTEGER, INTENT(out) :: status   ! Status flag: (0) no error
  
! !CALLING SEQUENCE:
! Called by: user code, usually init_dim_obs_f
! Calls: MPI_Allreduce
!EOP

! **********************************************************
! *** Compute global sum of local observation dimensions ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_Allreduce(val_p, val_g, 1, mpitype, mpiop, &
          COMM_filter, MPIerr)

     status = MPIerr
  ELSE
     val_g = val_p
     status = 0
  END IF



END SUBROUTINE PDAF_allreduce
