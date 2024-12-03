!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!mod_assimilation.F90: TSMP-PDAF implementation of routine
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
! Implementation for TSMP-PDAF.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP

! *** Model- and data specific variables ***

  INTEGER :: dim_state           ! Global model state dimension
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain

  INTEGER :: dim_obs_p                    ! Process-local number of observations
  REAL, ALLOCATABLE    :: obs_p(:)        ! Vector holding process-local observations
  INTEGER, ALLOCATABLE :: obs_index_p(:)  ! Vector holding state-vector indices of observations
  ! REAL, ALLOCATABLE    :: obs_f(:)        ! Vector holding full vector of observations
  ! REAL, ALLOCATABLE :: coords_obs_f(:,:)  ! Array for full observation coordinates

! *** Variables specific for TSMP-PDAF ***

  ! gw
  INTEGER, ALLOCATABLE :: dim_state_p_count(:) !Vector holding local state vector dimensions for processors of a single model communicator
  ! gw end
  REAL, ALLOCATABLE    :: obs(:)          ! Vector holding all observations for Global domain 
  INTEGER, ALLOCATABLE :: obs_index_l(:)  ! Vector holding local state-vector indices of observations
  INTEGER, ALLOCATABLE :: obs_interp_indices_p(:,:)  ! Vector holding state-vector indices of grid cells surrounding interpolation for PE-local domain
  INTEGER, ALLOCATABLE :: obs_interp_weights_p(:,:)  ! Vector holding weights of grid cells surrounding observation for PE-local domain
  INTEGER, ALLOCATABLE :: local_dims_obs(:) ! Array for process-local observation dimensions
  INTEGER, ALLOCATABLE :: local_disp_obs(:) ! Observation displacement array for gathering. Displacement: #obs before current PE
  ! pdaf-ordered index: determined by domain-decomposition
  ! nc-ordered index:   pure ordering of observation in NetCDF observation file
  INTEGER, ALLOCATABLE :: obs_pdaf2nc(:)  ! index mapping from a pdaf-ordered index to a nc-ordered index
  REAL, ALLOCATABLE :: pressure_obserr_p(:) ! Vector holding observation errors for paraflow run at each PE-local domain 
  !hcp
  !type :: scoltype
  !     integer, dimension(:), allocatable :: scol_obs_in
  !endtype
  integer, dimension(:,:), allocatable :: sc_p  !soil moisture of a soil column distributed over the procs of PE-local (i_z, i_obs)
  real, allocatable    :: idx_obs_nc_p(:)        
  INTEGER :: toffset      ! offset time step to shift all the assimilation steps
  !end hcp  
  REAL, ALLOCATABLE :: clm_obserr_p(:)    ! Vector holding  observation errors for CLM run at each PE-local domain  
  REAL, ALLOCATABLE :: distance(:)        ! Localization distance
  INTEGER, ALLOCATABLE :: global_to_local(:)  ! Vector to map global index to local domain index
  INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:) ! longitude and latitude of grid cells and observation cells
  INTEGER, ALLOCATABLE :: longxy_obs_floor(:), latixy_obs_floor(:) ! indices of grid cells with smaller lon/lat than observation location
  INTEGER, ALLOCATABLE :: var_id_obs(:)   ! for remote sensing data the variable identifier to group  
                                          ! variables distributed over a grid surface area 
  !kuw
  INTEGER, ALLOCATABLE :: obs_id_p(:) ! ID of observation point in PE-local domain
  INTEGER, ALLOCATABLE :: obs_nc2pdaf_deprecated(:)   ! index for mapping mstate to local domain
  INTEGER, ALLOCATABLE :: obs_nc2pdaf(:)   ! index for mapping mstate to local domain
  !kuw end

  ! Multi-scale DA

  ! store the maximum and minimum limits for remote sensing data with
  ! the same variable identity (CLM: coordinate limits, ParFlow: index
  ! limits)
  INTEGER, ALLOCATABLE :: maxlon(:) !CLM
  INTEGER, ALLOCATABLE :: minlon(:)
  INTEGER, ALLOCATABLE :: maxlat(:)
  INTEGER, ALLOCATABLE :: minlat(:)
  INTEGER, ALLOCATABLE :: maxix(:) !ParFlow
  INTEGER, ALLOCATABLE :: minix(:)
  INTEGER, ALLOCATABLE :: maxiy(:)
  INTEGER, ALLOCATABLE :: miniy(:)

  ! Dimension along the x and y direction for remote sensing data
  INTEGER :: dim_nx
  INTEGER :: dim_ny

  REAL, ALLOCATABLE :: lon_var_id(:)
  REAL, ALLOCATABLE :: lat_var_id(:)
  REAL, ALLOCATABLE :: ix_var_id(:)
  REAL, ALLOCATABLE :: iy_var_id(:)

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

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
  INTEGER :: filtertype   ! Select filter algorithm:
                          !   SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                          !   ESTKF (6), LESTKF (7), LEnKF (8)
                          !   NETF (9), LNETF (10), LKNETF(11), PF (12), 3DVAR (200)
  INTEGER :: subtype      ! Subtype of filter algorithm
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
                          !     (0) LETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   ESTKF:
                          !     (0) standard ESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LESTKF:
                          !     (0) standard LESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   LEnKF:
                          !     (0) Standard form of EnKF with covariance localization
                          !   NETF:
                          !     (0) standard NETF 
                          !   LNETF:
                          !     (0) standard LNETF 
                          !   LKNETF:
                          !     (0) HNK: 2-step LKNETF with NETF before LETKF
                          !     (1) HKN: 2-step LKNETF with LETKF before NETF
                          !     (4) HSync: LKNETF synchronous
                          !     (5) Offline mode - HNK: 2-step LKNETF with NETF before LETKF
                          !   PF:
                          !     (0) standard PF 
                          !   3D-Var:
                          !     (0) parameterized 3D-Var
                          !     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !     (4) 3D Ensemble Var using ESTKF for ensemble update
                          !     (6) hybrid 3D-Var using LESTKF for ensemble update
                          !     (7) hybrid 3D-Var using ESTKF for ensemble update
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
                          !  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                          !   (0) fixed
                          !   (1) global adaptive
                          !   (2) local adaptive for LSEIK/LETKF/LESTKF
                          !  NETF/LNETF/PF
                          !   (0) apply inflation on forecast ensemble
                          !   (2) apply inflation on analysis ensemble
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
!    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
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
                           ! ESTKF/LESTKF:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! NETF/LNETF:
                           ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           ! (1) use identity transformation
                           ! LKNETF:
                           ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           ! (1) use identity transformation
!    ! LSEIK/LETKF/LESTKF
  REAL    :: cradius       ! Cut-off radius for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                           ! For LESTKF, LETKF, and LSEIK
                           !   (0) constant weight of 1
                           !   (1) exponentially decreasing with SRADIUS
                           !   (2) use 5th-order polynomial
                           !   (3) regulated localization of R with mean error variance
                           !   (4) regulated localization of R with single-point error variance
                           ! For LEnKF
                           !   (0) constant weight of 1
                           !   (1) exponentially decreasing with SRADIUS
                           !   (2) 5th-order polynomial weight function
  REAL    :: sradius       ! Support radius for 5th order polynomial
                           !   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                           !   (0) symmetric square root
                           !   (1) Cholesky decomposition
!    ! NETF/LNETF
  INTEGER :: type_winf     ! Set weights inflation: (1) activate
                           !   (0) no weights inflation
                           !   (1) use N_eff/N>limit_winf
  REAL    :: limit_winf    ! Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  INTEGER :: type_hyb      ! Type of hybrid weight:
                    !   (0) use fixed hybrid weight hyb_gamma
                    !   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                    !   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                    !   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                    !   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  REAL    :: hyb_gamma     ! Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  REAL    :: hyb_kappa     ! Hybrid norm for using skewness and kurtosis
!    ! Particle filter
  INTEGER :: pf_res_type   ! Resampling type for PF
                           ! (1) probabilistic resampling
                           ! (2) stochastic universal resampling
                           ! (3) residual resampling        
  INTEGER :: pf_noise_type    ! Resampling type for PF
                           ! (0) no perturbations, (1) constant stddev, 
                           ! (2) amplitude of stddev relative of ensemble variance
  REAL :: pf_noise_amp     ! Noise amplitude (>=0.0, only used if pf_noise_type>0)
!    ! 3D-Var
  INTEGER :: type_opt      !< * Type of minimizer for 3DVar
                           !<   (1) LBFGS (default)
                           !<   (2) CG+
                           !<   (3) plain CG
                           !<   (12) CG+ parallelized
                           !<   (13) plain CG parallelized
  INTEGER :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  INTEGER :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  INTEGER :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  REAL :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
  INTEGER :: solver_iparam1 = 2 ! Solver specific parameter
                                !  LBFGS: parameter m (default=5)
                                !       Number of corrections used in limited memory matrix; 3<=m<=20
                                !  CG+: parameter method (default=2)
                                !       (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
                                !  CG: maximum number of iterations (default=200)
  INTEGER :: solver_iparam2 = 1 ! Solver specific parameter
                                !  LBFGS: - not used - 
                                !  CG+: parameter irest (default=1)
                                !       (0) no restarts; (n>0) restart every n steps
                                !  CG: - not used -
  REAL :: solver_rparam1 = 1.0e-6 ! Solver specific parameter
                                !  LBFGS: limit for stopping iterations 'pgtol' (default=1.0e-5)
                                !  CG+: convergence parameter 'eps' (default=1.0e-5)
                                !  CG: conpergence parameter 'eps' (default=1.0e-6)
  REAL :: solver_rparam2 = 1.0e+7 ! Solver specific parameter
                                !  LBFGS: tolerance in termination test 'factr' (default=1.0e+7) 
                                !  CG+: - not used -
                                !  CG: - not used -

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time          ! model time
!   REAL :: coords_l(2)      ! Coordinates of local analysis domain
!   INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
!   INTEGER, ALLOCATABLE :: id_lobs_in_fobs(:)  ! Indices of local observations in full obs. vector
!   REAL, ALLOCATABLE    :: distance_l(:)   ! Vector holding distances of local observations

! !$OMP THREADPRIVATE(coords_l, id_lstate_in_pstate, id_lobs_in_fobs, distance_l)

END MODULE mod_assimilation
