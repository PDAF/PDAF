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
! !ROUTINE: PDAF_gather_obs_f2_flex --- Gather a full observation array
!
! !INTERFACE:
SUBROUTINE PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p, coords_f, &
     nrows, status)

! !DESCRIPTION:
! If the local filter is used with a domain-decomposed model,
! the observational information from different sub-domains
! has to be combined into the full observation vector. 
! In this routine the process-local parts of a coordinate array
! accompanying the observation vector are gathered into a full
! array of coordinates. 
! The routine is for the case that the observation coordinates
! are stored column-wise, i.e. each column is the set of coordinates
! for one observation. This should be the usual case, as in this
! case the set of coordinates of one observations are stored
! next to each other in memory. If the coordinates are stored row-
! wise, the routine PDAF_gather_obs_f can be used, but has to be
! called separately for each column. 
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
       ONLY: COMM_filter, MPIerr, mype_filter, npes_filter

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local observation dimension
  INTEGER, INTENT(in) :: dim_obs_f    ! Full observation dimension
  INTEGER, INTENT(in) :: nrows        ! Number of rows in array
  REAL, INTENT(in)  :: coords_p(nrows, dim_obs_p)  ! PE-local array
  REAL, INTENT(out) :: coords_f(nrows, dim_obs_f)  ! Full gathered array
  INTEGER, INTENT(out) :: status   ! Status flag: (0) no error

! !CALLING SEQUENCE:
! Called by: user code
! Calls: MPI_Allreduce
! Calls: MPI_Allgather
! Calls: MPI_AllgatherV
!EOP

! local variables
  INTEGER :: i                              ! Counter
  INTEGER :: dimobs_f                       ! full dimension of observation vector obtained from allreduce
  INTEGER, ALLOCATABLE :: all_dim_obs_p(:)  ! PE-Local observation dimensions
  INTEGER, ALLOCATABLE :: all_dim_obs_p2(:) ! local-dims for multi-row array
  INTEGER, ALLOCATABLE :: all_dis_obs_p2(:) ! displacements to gather multi-row array


! **********************************************************
! *** Compute global sum of local observation dimensions ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_Allreduce(dim_obs_p, dimobs_f, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     dimobs_f = dim_obs_p
  END IF


! ****************************************************************************
! *** Gather and store array of process-local dimensions and displacements ***
! ****************************************************************************

  ALLOCATE(all_dim_obs_p(npes_filter))

  IF (npes_filter>1) THEN
     CALL MPI_Allgather(dim_obs_p, 1, MPI_INTEGER, all_dim_obs_p, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)
  ELSE
     all_dim_obs_p = dim_obs_p
  END IF


! **********************************************************
! *** Gather full observation coordinates array          ***
! **********************************************************

  IF (npes_filter>1) THEN
     ALLOCATE(all_dis_obs_p2(npes_filter))
     ALLOCATE(all_dim_obs_p2(npes_filter))

     ! Init array of local dimensions
     do i = 1, npes_filter
        all_dim_obs_p2(i) = nrows * all_dim_obs_p(i)
     end do

     ! Init array of displacements for observation vector
     all_dis_obs_p2(1) = 0
     DO i = 2, npes_filter
        all_dis_obs_p2(i) = all_dis_obs_p2(i-1) + all_dim_obs_p2(i-1)
     END DO

     CALL MPI_AllGatherV(coords_p, all_dim_obs_p2(mype_filter+1), MPI_REALTYPE, &
          coords_f, all_dim_obs_p2, all_dis_obs_p2, MPI_REALTYPE, &
          COMM_filter, MPIerr)

     DEALLOCATE(all_dim_obs_p2, all_dis_obs_p2)

     status = MPIerr
  ELSE
     coords_f = coords_p
     
     status = 0
  END IF

END SUBROUTINE PDAF_gather_obs_f2_flex
