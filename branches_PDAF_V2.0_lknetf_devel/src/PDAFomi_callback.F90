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

!> Generic PDAFomi callback routines
!!
!! This file provides generic call-back routines for OMI. The OMI structure
!! allow us to provide several of the observation-related call-back 
!! routines for PDAF in a generic form. These routines are collected in
!! this file.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!! 
!! __Revision history:__
!! * 2020-06 - Lars Nerger - Initial code splitting callback_obs_pdafomi.F90
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_obs_f
!!
!! This routine calls the routine PDAFomi_init_obs_f
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obs_f_cb(step, dim_obs_f, observation_f)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_f

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector

! *** local variables ***
  INTEGER :: i                ! Loop counter
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  ! Initialize offset (it will be incremented in PDAFomi_init_obs_f)
  offset_obs_f = 0

  ! The order of the calls has to be consistent with those in obs_op_f_pdafomi
  DO i=1, n_obstypes
     CALL PDAFomi_init_obs_f(obs_f_all(i)%ptr, dim_obs_f, observation_f, offset_obs_f)
  END DO

END SUBROUTINE PDAFomi_init_obs_f_cb



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar
!!
!! This routine calls the routine PDAFomi_init_obsvar_f
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obsvar_cb(step, dim_obs_p, obs_p, meanvar)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_f

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
  REAL, INTENT(out)   :: meanvar       !< Mean observation error variance

! *** Local variables ***
  INTEGER :: i                ! Loop counter
  INTEGER :: cnt_obs_f        ! Count observations for offset


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_f = 0

  DO i=1, n_obstypes
     CALL PDAFomi_init_obsvar_f(obs_f_all(i)%ptr, meanvar, cnt_obs_f)
  END DO

END SUBROUTINE PDAFomi_init_obsvar_cb



!-------------------------------------------------------------------------------
!> Call-back routine for g2l_obs
!!
!! This routine calls the routine PDAFomi_g2l_obs
!! for each observation type
!!
SUBROUTINE PDAFomi_g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_g2l_obs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local observation vector
  INTEGER, INTENT(in) :: dim_obs_l  !< Dimension of local observation vector
  REAL, INTENT(in)    :: ostate_f(dim_obs_f)   !< Full PE-local obs.ervation vector
  REAL, INTENT(out)   :: ostate_l(dim_obs_l)   !< Observation vector on local domain

! *** local variables ***
  INTEGER :: i                      ! Loop counter


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  DO i=1, n_obstypes
     CALL PDAFomi_g2l_obs(obs_l_all(i)%ptr, obs_f_all(i)%ptr, ostate_f, ostate_l)
  END DO

END SUBROUTINE PDAFomi_g2l_obs_cb



!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_l
!!
!! This routine calls the routine PDAFomi_init_obs_l
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obs_l_cb(domain_p, step, dim_obs_l, observation_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector

! *** local variables ***
  INTEGER :: i                      ! Loop counter


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  DO i=1, n_obstypes
     CALL PDAFomi_init_obs_l(obs_l_all(i)%ptr, obs_f_all(i)%ptr,  observation_l)
  END DO

END SUBROUTINE PDAFomi_init_obs_l_cb



!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar_l
!!
!! This routine calls the routine PDAFomi_init_obsvar_l
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l, meanvar_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance

! *** Local variables ***
  INTEGER :: i                         ! Loop counter
  INTEGER :: cnt_obs_l                 ! Count local observations for offset


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_l = 0

  ! The order of the calls is not relevant
  DO i=1, n_obstypes
     CALL PDAFomi_init_obsvar_l(obs_l_all(i)%ptr, obs_f_all(i)%ptr, meanvar_l, cnt_obs_l)
  END DO

END SUBROUTINE PDAFomi_init_obsvar_l_cb



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA_l
!!
!! This routine calls the routine PDAFomi_prodRinvA_l
!! for each observation type
!!
SUBROUTINE PDAFomi_prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA_l

  ! Include filter process rank
  USE PDAF_mod_filterMPI, ONLY: mype_filter
  ! Include verbosity information
  USE PDAF_mod_filter, ONLY: screen
#if defined (_OPENMP)
  ! Include OpenMP function to determine verbosity for OpenMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix

! *** local variables ***
  INTEGER :: i                       ! Loop counter
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  INTEGER, SAVE :: mythread          ! Thread variable for OpenMP

!$OMP THREADPRIVATE(mythread, domain_save)


! **********************
! *** INITIALIZATION ***
! **********************

  ! For OpenMP parallelization, determine the thread index
#if defined (_OPENMP)
  mythread = omp_get_thread_num()
#else
  mythread = 0
#endif

  ! Set verbosity flag (Screen output for first analysis domain)
  IF (screen > 0) THEN
     IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
        verbose = 1

        ! In case of OpenMP, let only thread 0 write output to the screen
        IF (mythread>0) verbose = 0
     ELSE
        verbose = 0
     END IF
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p


! *************************************
! *** Compute                       ***
! ***                  -1           ***
! ***           C = W R   A         ***
! ***                               ***
! *** where W are the localization  ***
! *** weights.                      ***
! *************************************

  DO i=1, n_obstypes
     CALL PDAFomi_prodRinvA_l(obs_l_all(i)%ptr, obs_f_all(i)%ptr, dim_obs_l, rank, &
          A_l, C_l, verbose)
  END DO
  
END SUBROUTINE PDAFomi_prodRinvA_l_cb


!-------------------------------------------------------------------------------
!> Call-back routine for likelihood_l
!!
!! This routine calls the routine PDAFomi_likelihood_l
!! for each observation type
!!
SUBROUTINE PDAFomi_likelihood_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all, obs_l_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_likelihood_l

  ! Include filter process rank
  USE PDAF_mod_filterMPI, ONLY: mype_filter
  ! Include verbosity information
  USE PDAF_mod_filter, ONLY: screen
#if defined (_OPENMP)
  ! Include OpenMP function to determine verbosity for OpenMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p           ! Current local analysis domain
  INTEGER, INTENT(in) :: step               !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l          !< PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs_l(dim_obs_l)   !< PE-local vector of observations
  REAL, INTENT(inout) :: resid_l(dim_obs_l) !< Input vector of residuum
  REAL, INTENT(out)   :: lhood_l            !< Output vector - log likelihood

! *** local variables ***
  INTEGER :: i                       ! Loop counter
  INTEGER :: verbose                 ! Verbosity flag
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  INTEGER, SAVE :: mythread          ! Thread variable for OpenMP

!$OMP THREADPRIVATE(mythread, domain_save)


! **********************
! *** INITIALIZATION ***
! **********************

  ! For OpenMP parallelization, determine the thread index
#if defined (_OPENMP)
  mythread = omp_get_thread_num()
#else
  mythread = 0
#endif

  ! Set verbosity flag (Screen output for first analysis domain)
  IF (screen > 0) THEN
     IF ((domain_p < domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
        verbose = 1

        ! In case of OpenMP, let only thread 0 write output to the screen
        IF (mythread>0) verbose = 0
     ELSE
        verbose = 0
     END IF
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p


! ********************************
! *** Compute local likelihood ***
! ********************************

  ! Initialize likelihood value before starting computation
  lhood_l = 0.0

  ! Increment likelihood
  DO i=1, n_obstypes
     CALL PDAFomi_likelihood_l(obs_l_all(i)%ptr, obs_f_all(i)%ptr, resid_l, &
          lhood_l, verbose)
  END DO

END SUBROUTINE PDAFomi_likelihood_l_cb



!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA
!!
!! This routine calls the routine PDAFomi_prodRinvA
!! for each observation type
!!
SUBROUTINE PDAFomi_prodRinvA_cb(step, dim_obs_p, ncol, obs_p, A_p, C_p)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  INTEGER, INTENT(in) :: ncol              !< Number of columns in A_p and C_p
  REAL, INTENT(in)    :: obs_p(dim_obs_p)  !< PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p, ncol) !< Input matrix
  REAL, INTENT(out)   :: C_p(dim_obs_p, ncol) !< Output matrix

! *** local variables ***
  INTEGER :: i                ! Loop counter


! *************************************
! *** Compute                       ***
! ***                -1             ***
! ***           C = R   A           ***
! *************************************

  DO i=1, n_obstypes
     CALL PDAFomi_prodRinvA(obs_f_all(i)%ptr, ncol, A_p, C_p)
  END DO
  
END SUBROUTINE PDAFomi_prodRinvA_cb



!-------------------------------------------------------------------------------
!> Call-back routine for likelihood
!!
!! This routine calls the routine PDAFomi_likelihood
!! for each observation type
!!
SUBROUTINE PDAFomi_likelihood_cb(step, dim_obs, obs, resid, lhood)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_likelihood

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: dim_obs          !< PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs(dim_obs)     !< PE-local vector of observations
  REAL, INTENT(in)    :: resid(dim_obs)   !< Input vector of residuum
  REAL, INTENT(out)   :: lhood            !< Output vector - log likelihood

! *** local variables ***
  INTEGER :: i                ! Loop counter


! **************************
! *** Compute likelihood ***
! **************************

  ! Initialize likelihood value before starting computation
  lhood = 0.0

  ! Increment likelihood
  DO i=1, n_obstypes
     CALL PDAFomi_likelihood(obs_f_all(i)%ptr, dim_obs, obs, resid, lhood)
  END DO

END SUBROUTINE PDAFomi_likelihood_cb



!-------------------------------------------------------------------------------
!> Call-back routine for add_obs_error
!!
!! This routine calls the routine PDAFomi_add_obs_error
!! for each observation type
!!
SUBROUTINE PDAFomi_add_obs_error_cb(step, dim_obs_p, C_p)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_add_obs_error

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  REAL, INTENT(inout) :: C_p(dim_obs_p,dim_obs_p) ! Matrix to which R is added

! *** local variables ***
  INTEGER :: i                ! Loop counter


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

  DO i=1, n_obstypes
     CALL PDAFomi_add_obs_error(obs_f_all(i)%ptr, dim_obs_p, C_p)
  END DO
  
END SUBROUTINE PDAFomi_add_obs_error_cb



!-------------------------------------------------------------------------------
!> Call-back routine for init_obscovar
!!
!! This routine calls the routine PDAFomi_init_obscovar
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obscovar_cb(step, dim_obs, dim_obs_p, covar, m_state_p, &
     isdiag)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obscovar

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of observation vector
  INTEGER, INTENT(in) :: dim_obs_p            !< PE-local dimension of obs. vector
  REAL, INTENT(out) :: covar(dim_obs,dim_obs) !< Observation error covar. matrix
  REAL, INTENT(in) :: m_state_p(dim_obs_p)    !< Observation vector
  LOGICAL, INTENT(out) :: isdiag              !< Whether matrix R is diagonal

! *** local variables ***
  INTEGER :: i                ! Loop counter


! *************************************
! ***   Initialize covariances      ***
! *************************************

  DO i=1, n_obstypes
     CALL PDAFomi_init_obscovar(obs_f_all(i)%ptr, dim_obs_p, covar, isdiag)
  END DO
  
END SUBROUTINE PDAFomi_init_obscovar_cb



!-------------------------------------------------------------------------------
!> Call-back routine for init_obserr_f_pdaf
!!
!! This routine calls the routine PDAFomi_init_obserr_f
!! for each observation type
!!
SUBROUTINE PDAFomi_init_obserr_f_cb(step, dim_obs_f, obs_f, obserr_f)

  ! Include overall pointer to observation variables
  use PDAFomi, only: n_obstypes, obs_f_all
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obserr_f

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f           !< Full dimension of observation vector
  REAL, INTENT(in)    :: obs_f(dim_obs_f)    !< Full observation vector
  REAL, INTENT(out)   :: obserr_f(dim_obs_f) !< Full observation error stddev

! *** local variables ***
  INTEGER :: i                ! Loop counter


! *****************************************************************************
! *** Initialize vector of observation errors for generating synthetic obs. ***
! *****************************************************************************

  DO i=1, n_obstypes
     CALL PDAFomi_init_obserr_f(obs_f_all(i)%ptr, obserr_f)
  END DO
  
END SUBROUTINE PDAFomi_init_obserr_f_cb
