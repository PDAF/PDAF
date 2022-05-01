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
! !ROUTINE: PDAF_gather_obs_f --- Gather a full observation vector
!
! !INTERFACE:
SUBROUTINE PDAF_gather_obs_f(obs_p, obs_f, status)

! !DESCRIPTION:
! If the local filter is used with a domain-decomposed model,
! the observational information from different sub-domains
! has to be combined into the full observation vector. 
! In this routine the process-local parts of the observation
! vector are gathered into a full observation vector. 
! The routine requires that PDAF_gather_dim_obs_f was executed
! before, because this routine initializes dimensions that are 
! used here. 
! The routine can also be used to gather full arrays of coordinates.
! It is however, only usable if the coordinates are stored row-
! wise, i.e. each row represents the set of coordinates for one
! observation point. It has to be called separately for each column. 
! A  better alternative is the row-wise storage of coordinates. In this
! case the routine PDAF_gather_dim_obs_f allows the gather the full
! coordinate array in one step.
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
       all_dim_obs_p, all_dis_obs_p, dimobs_p, dimobs_f

  IMPLICIT NONE
  
! !ARGUMENTS:
  REAL, INTENT(in)  :: obs_p(dimobs_p)  ! PE-local vector
  REAL, INTENT(out) :: obs_f(dimobs_f)  ! Full gathered vector
  INTEGER, INTENT(out) :: status   ! Status flag: 
                                   ! (0) no error
                                   ! (1) when PDAF_gather_dim_obs_f not executed before

! !CALLING SEQUENCE:
! Called by: user code, usually init_dim_obs_f and obs_op_f
! Calls: MPI_AllgatherV
!EOP


! **********************************************************
! *** Gather full observation vector                     ***
! **********************************************************

  IF (ALLOCATED(all_dim_obs_p)) THEN

     IF (npes_filter>1) THEN
        CALL MPI_AllGatherV(obs_p, all_dim_obs_p(mype_filter+1), MPI_REALTYPE, &
             obs_f, all_dim_obs_p, all_dis_obs_p, MPI_REALTYPE, &
             COMM_filter, MPIerr)
     ELSE
        obs_f = obs_p
     END IF

     status = 0

  ELSE
     ! ERROR: all_dim_obs_p not allocated 
     ! probably PDAF_gather_dim_obs_f was not run before
     status = 1
     
  END IF

END SUBROUTINE PDAF_gather_obs_f
