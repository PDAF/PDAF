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
!$Id: PDAF-D_enkf_obs_ensemble.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_obs_ensemble --- Generate ensemble of observations for EnKF
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens, m_ens_p, &
     U_init_obs, U_init_obs_covar, screen, flag)

! !DESCRIPTION:
! This routine generates an ensemble of observations 
! from a mean observation with prescribed error
! covariance matrix for the EnKF94/98.
! For a diagonal matrix the column vectors are 
! directly used. For a non-diagonal matrix the 
! eigen-decomposition is computed.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2001-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter, MPI_INTEGER

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p  ! Local dimension of current observation
  INTEGER, INTENT(in) :: dim_obs    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens    ! Size of ensemble
  REAL, INTENT(out)   :: m_ens_p(dim_obs_p,dim_ens) ! PE-local obs. ensemble 
  INTEGER, INTENT(in) :: screen     ! Verbosity flag
  INTEGER, INTENT(inout) :: flag    ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_obs, & ! Initialize observation vector
       U_init_obs_covar     ! Initialize observation error covariance matrix

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm
! Called by: PDAF_enkf_analysis_rsm
! Calls: U_init_obs
! Calls: U_init_obs_covar
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: larnvTYPE (BLAS; dlarnv or slarnv dependent on precision)
! Calls: MPI_allgather (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, j, member         ! Counters
  REAL :: randval                 ! Value of random number
  REAL :: m_state_p(dim_obs_p)    ! Observation vector
  REAL :: covar(dim_obs, dim_obs) ! Observation covariance matrix
  INTEGER :: syev_info            ! Output flag of eigenproblem routine
  INTEGER, SAVE :: allocflag = 0  ! Flag for first-time allocation
  INTEGER, SAVE :: iseed(4)       ! Seed for random number generator LARNV
  INTEGER, SAVE :: first = 1      ! Flag for setting of random-number seed
  LOGICAL :: isdiag               ! Is the observation error cov. matrix diagonal?
  REAL, ALLOCATABLE :: eigenv(:)  ! Vector of eigenvalues
  REAL, ALLOCATABLE :: workarray(:) ! Workarray for eigenproblem routine
  INTEGER, ALLOCATABLE :: local_dim_obs(:) ! Array of local dimensions
  INTEGER, ALLOCATABLE :: local_dis(:)     ! Array of local displacements


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- Generate ensemble of observations'

  IF (first == 1) THEN
     ! Initialize seed
     iseed(1) = 1
     iseed(2) = 5
     iseed(3) = 7
     iseed(4) = 9
     first = 2
  END IF

  ! allocate memory for temporary fields
  ALLOCATE(eigenv(dim_obs))
  ALLOCATE(workarray(3 * dim_obs))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', 4 * dim_obs)
     allocflag = 1
  END IF


! *************************************
! *** generate observation ensemble ***
! *************************************

  ! *** get current observation vector ***
  CALL U_init_obs(step, dim_obs_p, m_state_p)

  ! *** Get current observation covariance matrix ***
  ! *** We initialize the global observation error covariance matrix
  ! *** to avoid a parallelization of the possible eigendecomposition.
  CALL U_init_obs_covar(step, dim_obs, dim_obs_p, covar, m_state_p, &
       isdiag)

  diagA: IF (.NOT. isdiag) THEN
     ! *** compute Eigendecomposition of covariance matrix ***
     ! *** We use the LAPACK routine SYEV                  ***
     CALL syevTYPE('v', 'l', dim_obs, covar, dim_obs, &
          eigenv, workarray, 3 * dim_obs, syev_info)
  ELSE diagA
     ! *** Do not perform eigendecomposition here, since COVAR is diagonal
     IF (mype == 0 .AND. screen > 0) &
         WRITE (*, '(a, 5x, a)') 'PDAF', '--- use diagonal observation eror cov. matrix'
     syev_info = 0
  END IF diagA

  ! check if eigendecomposition was successful
  ensemble: IF (syev_info /= 0) THEN
     ! Eigendecomposition failed

     WRITE (*, '(/5x, a/)') &
          'PDAF-ERROR(3): Problem in eigendecomposition of observation covariance !!!'
     flag = 3

  ELSE
     ! Eigendecomposition OK, continue with ensemble generation

     diagB: IF (.NOT. isdiag) THEN
        ! rescale eigenvectors (if EVP was computed)
        DO j = 1, dim_obs
           DO i = 1, dim_obs
              covar(i, j) = covar(i, j) * SQRT(eigenv(j))
           END DO
        END DO
     ELSE diagB
        ! Only compute square-roots of variances here
        DO i = 1, dim_obs
           covar(i, i) = SQRT(covar(i, i))
        END DO
     END IF diagB

     ! gather array of observation dimensions
     ALLOCATE(local_dim_obs(npes_filter))
     ALLOCATE(local_dis(npes_filter))

     CALL MPI_allgather(dim_obs_p, 1, MPI_INTEGER, local_dim_obs, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)

     ! Init array of displacements
     local_dis(1) = 0
     DO i = 2, npes_filter
        local_dis(i) = local_dis(i - 1) + local_dim_obs(i - 1)
     END DO

     ! generate random states for local domain
     members: DO member = 1, dim_ens
        m_ens_p(:, member) = m_state_p(:)
        eigenvectors: DO j = 1, dim_obs
           CALL larnvTYPE(3, iseed, 1, randval)
           components: DO i = 1, dim_obs_p
              m_ens_p(i, member) = m_ens_p(i, member) &
                   + covar(i + local_dis(mype + 1), j) * randval
           END DO components
        END DO eigenvectors
     END DO members

     DEALLOCATE(local_dim_obs, local_dis)
  END IF ensemble


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(eigenv, workarray)

END SUBROUTINE PDAF_enkf_obs_ensemble
