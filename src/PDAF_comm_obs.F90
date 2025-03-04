! Copyright (c) 2004-2025 Lars Nerger
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
!> Routines for communicating observations when not using OMI
!!
!! When not using PDAF-OMI one needs to gather observation information
!! on the user side. This module provides functionality for
!! these parallel operations.
!!
!! !  This are core routines of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2025-03 - Lars Nerger - Initial code by joining several files
!! Other revisions - see repository log
!!
MODULE PDAF_comm_obs

CONTAINS

!--------------------------------------------------------------------------
!> get full observation dimension
!!
!! If the local filter is used with a domain-decomposed model,
!! the observational information from different sub-domains
!! has to be combined into the full observation vector. 
!! This routine is called as a first step to compute the
!! full observation dimension from the process-local
!! observation dimensions.
!! Practically, the operation is a simple MPI_Allreduce. This
!! is encapsulated here to simplify the operation for the users. 
!! In addition an array storing the pe-local observation dimensions
!! and an array of displacements for gathering are initialized.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2017-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, npes_filter, &
         all_dim_obs_p, all_dis_obs_p, dimobs_p, dimobs_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: dim_obs_p    !< PE-local observation dimension
    INTEGER, INTENT(out) :: dim_obs_f    !< Full observation dimension

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


!--------------------------------------------------------------------------
!> Gather a full observation vector
!!
!! If the local filter is used with a domain-decomposed model,
!! the observational information from different sub-domains
!! has to be combined into the full observation vector. 
!! In this routine the process-local parts of the observation
!! vector are gathered into a full observation vector. 
!! The routine requires that PDAF_gather_dim_obs_f was executed
!! before, because this routine initializes dimensions that are 
!! used here. 
!!
!! The routine can also be used to gather full arrays of coordinates.
!! It is however, only usable if the coordinates are stored row-
!! wise, i.e. each row represents the set of coordinates for one
!! observation point. It has to be called separately for each column. 
!! A  better alternative is the row-wise storage of coordinates. In this
!! case the routine PDAF_gather_dim_obs_f allows the gather the full
!! coordinate array in one step.
!!
!! __Revision history:__
!! * 2017-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_obs_f(obs_p, obs_f, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, mype_filter, npes_filter, &
         all_dim_obs_p, all_dis_obs_p, dimobs_p, dimobs_f

    IMPLICIT NONE

! *** Arguments
    REAL, INTENT(in)  :: obs_p(dimobs_p)  !< PE-local vector
    REAL, INTENT(out) :: obs_f(dimobs_f)  !< Full gathered vector
    INTEGER, INTENT(out) :: status        !< Status flag: 
                                          !< * (0) no error
                                          !< * (1) when PDAF_gather_dim_obs_f not executed before


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


!--------------------------------------------------------------------------
!> Gather a full observation array
!!
!! If the local filter is used with a domain-decomposed model,
!! the observational information from different sub-domains
!! has to be combined into the full observation vector. 
!! In this routine the process-local parts of a coordinate array
!! accompanying the observation vector are gathered into a full
!! array of coordinates. 
!! The routine is for the case that the observation coordinates
!! are stored column-wise, i.e. each column is the set of coordinates
!! for one observation. This should be the usual case, as in this
!! case the set of coordinates of one observations are stored
!! next to each other in memory. If the coordinates are stored row-
!! wise, the routine PDAF_gather_obs_f can be used, but has to be
!! called separately for each column. 
!!
!! __Revision history:__
!! * 2017-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_obs_f2(coords_p, coords_f, nrows, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, mype_filter, npes_filter, &
         all_dim_obs_p, dimobs_p, dimobs_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nrows     !< Number of rows in array
    REAL, INTENT(in)  :: coords_p(nrows, dimobs_p)  !< PE-local array
    REAL, INTENT(out) :: coords_f(nrows, dimobs_f)  !< Full gathered array
    INTEGER, INTENT(out) :: status   !< Status flag: 
                                     !< (0) no error
                                     !< (1) when PDAF_gather dim_obs_f not executed before

! *** local variables ***
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

!--------------------------------------------------------------------------
!> Gather a full observation vector
!!
!! If the local filter is used with a domain-decomposed model,
!! the observational information from different sub-domains
!! has to be combined into the full observation vector. 
!! In this routine the process-local parts of the observation
!! vector are gathered into a full observation vector. 
!! The routine requires that PDAF_gather_dim_obs_f was executed
!! before, because this routine initializes dimensions that are 
!! used here. 
!! The routine can also be used to gather full arrays of coordinates.
!! It is however, only usable if the coordinates are stored row-
!! wise, i.e. each row represents the set of coordinates for one
!! observation point. It has to be called separately for each column. 
!! A  better alternative is the row-wise storage of coordinates. In this
!! case the routine PDAF_gather_dim_obs_f allows the gather the full
!! coordinate array in one step.
!!
!! __Revision history:__
!! * 2019-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, obs_f, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, mype_filter, npes_filter

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_p       !< PE-local observation dimension
    INTEGER, INTENT(in) :: dim_obs_f       !< Full observation dimension
    REAL, INTENT(in)  :: obs_p(dim_obs_p)  !< PE-local vector
    REAL, INTENT(out) :: obs_f(dim_obs_f)  !< Full gathered vector
    INTEGER, INTENT(out) :: status         !< Status flag: (0) no error

! *** Local variables ***
    INTEGER :: i                              ! Counter
    INTEGER :: dimobs_f                       ! full dimension of observation vector obtained from allreduce
    INTEGER, ALLOCATABLE :: all_dim_obs_p(:)  ! PE-Local observation dimensions
    INTEGER, ALLOCATABLE :: all_dis_obs_p(:)  ! PE-Local observation displacements


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


! **********************************************************
! *** Gather full observation vector                     ***
! **********************************************************

    IF (npes_filter>1) THEN
       CALL MPI_AllGatherV(obs_p, all_dim_obs_p(mype_filter+1), MPI_REALTYPE, &
            obs_f, all_dim_obs_p, all_dis_obs_p, MPI_REALTYPE, &
            COMM_filter, MPIerr)

       status = MPIerr
    ELSE
       obs_f = obs_p

       status = 0
    END IF


! ****************
! *** Clean up ***
! ****************

    DEALLOCATE(all_dim_obs_p, all_dis_obs_p)

  END SUBROUTINE PDAF_gather_obs_f_flex

!--------------------------------------------------------------------------
!> Gather a full observation array
!!
!! If the local filter is used with a domain-decomposed model,
!! the observational information from different sub-domains
!! has to be combined into the full observation vector. 
!! In this routine the process-local parts of a coordinate array
!! accompanying the observation vector are gathered into a full
!! array of coordinates. 
!! The routine is for the case that the observation coordinates
!! are stored column-wise, i.e. each column is the set of coordinates
!! for one observation. This should be the usual case, as in this
!! case the set of coordinates of one observations are stored
!! next to each other in memory. If the coordinates are stored row-
!! wise, the routine PDAF_gather_obs_f can be used, but has to be
!! called separately for each column. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, coords_p, coords_f, &
       nrows, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, mype_filter, npes_filter

    IMPLICIT NONE

! *** Arguments
    INTEGER, INTENT(in) :: dim_obs_p    !< PE-local observation dimension
    INTEGER, INTENT(in) :: dim_obs_f    !< Full observation dimension
    INTEGER, INTENT(in) :: nrows        !< Number of rows in array
    REAL, INTENT(in)  :: coords_p(nrows, dim_obs_p)  !< PE-local array
    REAL, INTENT(out) :: coords_f(nrows, dim_obs_f)  !< Full gathered array
    INTEGER, INTENT(out) :: status      !< Status flag: (0) no error

! *** local variables ***
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

!--------------------------------------------------------------------------
!> Perform MPI_allreduce of a single value in COMM_filter
!!
!! This routine simlpy performs an MPI_Allreduce for a single
!! value of the specified MPI type. It can be used to 
!! obtain the global number of observations.
!!
!! __Revision history:__
!! 2019-03 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF_allreduce(val_p, val_g, mpitype, mpiop, status)

!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, npes_filter

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: val_p    !< PE-local value
    INTEGER, INTENT(out) :: val_g    !< reduced global value
    INTEGER, INTENT(in)  :: mpitype  !< MPI data type
    INTEGER, INTENT(in)  :: mpiop    !< MPI operator
    INTEGER, INTENT(out) :: status   !< Status flag: (0) no error


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

END MODULE PDAF_comm_obs
