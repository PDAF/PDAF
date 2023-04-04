!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!mod_assimilation.F90: TerrSysMP-PDAF implementation of routine
!                     'mod_assimilation' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: mod_assimilation.F90 1444 2013-10-04 10:54:08Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assimilation

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! Implementation for the 2D online example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP

! *** Variables specific for online tutorial example ***

  INTEGER :: dim_state           ! Global model state dimension
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
  ! gw
  INTEGER, ALLOCATABLE :: dim_state_p_count(:)
  INTEGER, ALLOCATABLE :: dim_state_p_stride(:)
  ! gw end
  REAL, ALLOCATABLE    :: obs(:)          ! Vector holding all observations for Global domain 
  REAL, ALLOCATABLE    :: obs_p(:)        ! Vector holding observations for PE-local domain
  INTEGER, ALLOCATABLE :: obs_index_p(:)  ! Vector holding state-vector indices of observations for PE-local domain
  INTEGER, ALLOCATABLE :: obs_index_l(:)  ! Vector holding local state-vector indices of observations
  INTEGER, ALLOCATABLE :: obs_interp_indices_p(:,:)  ! Vector holding state-vector indices of grid cells surrounding interpolation for PE-local domain
  INTEGER, ALLOCATABLE :: obs_interp_weights_p(:,:)  ! Vector holding weights of grid cells surrounding observation for PE-local domain
  INTEGER, ALLOCATABLE :: local_dims_obs(:) ! Array for process-local observation dimensions
  INTEGER, ALLOCATABLE :: obs_nc2pdaf(:)  ! mapping ordering of obs between netcdf input and internal ordering in pdaf
  REAL, ALLOCATABLE :: pressure_obserr_p(:) ! Vector holding observation errors for paraflow run at each PE-local domain 
  REAL, ALLOCATABLE :: clm_obserr_p(:)    ! Vector holding  observation errors for CLM run at each PE-local domain  
  REAL, ALLOCATABLE :: distance(:)        ! Localization distance
  REAL, ALLOCATABLE :: xcoord_fortran_g(:), & ! Global coordinates for the domain are stored,
                       ycoord_fortran_g(:), & ! been gathered from local domains. used when 
                       zcoord_fortran_g(:)    ! local filter analysis is selected. 
  INTEGER, ALLOCATABLE :: global_to_local(:)  ! Vector to map global index to local domain index
  INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:) ! longitude and latitude of grid cells and observation cells
  INTEGER, ALLOCATABLE :: longxy_obs_floor(:), latixy_obs_floor(:) ! indices of grid cells with smaller lon/lat than observation location
  INTEGER, ALLOCATABLE :: var_id_obs(:)   ! for remote sensing data the variable identifier to group  
                                          ! variables distributed over a grid surface area 
  !kuw
  INTEGER, ALLOCATABLE :: obs_id_p(:) ! ID of observation point in PE-local domain
  INTEGER, ALLOCATABLE :: m_id_f(:)   ! index for mapping mstate to local domain
  !kuw end

  INTEGER :: dim_obs_p     ! Process-local number of observations
  INTEGER :: dim_obs_f     ! Full number of observations
  INTEGER, ALLOCATABLE :: maxlon(:), minlon(:), maxlat(:), minlat(:), & ! store the maximum and minimum coordinates limits 
                          maxix(:), minix(:), maxiy(:), miniy(:)        ! for remote sensing data with the same variable identity 
  INTEGER :: dim_nx, dim_ny ! the dimension along the x and y direction
                            ! for remote sensing data
  REAL, ALLOCATABLE :: lon_var_id(:), ix_var_id(:)
  REAL, ALLOCATABLE :: lat_var_id(:), iy_var_id(:)

  ! *** User defined observation filename ***
  character (len = 110) :: obs_filename

! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs      ! time step interval between assimilation steps
  REAL    :: rms_obs       ! RMS error size for observation generation
  INTEGER :: dim_obs       ! Number of observations
  INTEGER :: toffset      ! offset time step to shift all the assimilation steps

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
  INTEGER :: subtype      ! Subtype of filter algorithm
                          !   SEEK: 
                          !     (0) evolve normalized modes
                          !     (1) evolve scaled modes with unit U
                          !     (2) fixed basis (V); variable U matrix
                          !     (3) fixed covar matrix (V,U kept static)
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   EnKF:
                          !     (0) analysis for large observation dimension
                          !     (1) analysis for small observation dimension
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
!    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/LSEIK/ETKFS
  INTEGER :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
!    ! LSEIK/LETKF
  REAL    :: local_range   ! Range for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  REAL    :: srange        ! Support range for 5th order polynomial
                           !   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     ! For SEIK: Definition of ensemble covar matrix
                           ! (0): Factor (r+1)^-1 (or N^-1)
                           ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           ! This setting is only for the model part; The definition
                           ! of P has also to be specified in PDAF_filter_init.
                           ! Only for upward-compatibility of PDAF!
  REAL    :: time          ! model time

END MODULE mod_assimilation
