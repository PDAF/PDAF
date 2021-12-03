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
! !ROUTINE: PDAF_gather_obs_f2 --- Gather a full observation array
!
! !INTERFACE:
SUBROUTINE PDAF_gather_obs_f2(coords_p, coords_f, nrows, status)

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
! 2017-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: COMM_filter, MPIerr, mype_filter, npes_filter, &
       all_dim_obs_p, dimobs_p, dimobs_f

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: nrows     ! Number of rows in array
  REAL, INTENT(in)  :: coords_p(nrows, dimobs_p)  ! PE-local array
  REAL, INTENT(out) :: coords_f(nrows, dimobs_f)  ! Full gathered array
  INTEGER, INTENT(out) :: status   ! Status flag: 
                                   ! (0) no error
                                   ! (1) when PDAF_gather dim_obs_f not executed before

! !CALLING SEQUENCE:
! Called by: user code, usually init_dim_obs_f and obs_op_f
! Calls: MPI_AllgatherV
!EOP

! local variables
  INTEGER :: i                              ! Counter
  INTEGER, ALLOCATABLE :: all_dim_obs_p2(:) ! local-dims for multi-row array
  INTEGER, ALLOCATABLE :: all_dis_obs_p2(:) ! displacements to gather multi-row array


! **********************************************************
! *** Gather full observation coordinates array          ***
! **********************************************************

  IF (ALLOCATED(all_dim_obs_p)) THEN

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
     ELSE
        coords_f = coords_p
     END IF

     status = 0

  ELSE
     ! ERROR: all_dim_obs_p not allocated 
     ! probably PDAF_gather_dim_obs_f was not run before
     status = 1
     
  END IF

END SUBROUTINE PDAF_gather_obs_f2
