!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
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
!init_pdaf.F90: TSMP-PDAF implementation of routine
!               'init_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_pdaf.F90 1444 2013-10-04 10:54:08Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the online mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !TSMP-PDAF-DESCRIPTION:
! This routine initializes a pointer to the state vector that is set
! by component-model-specific routines in `initialize_tsmp`. 
!
! This routine sets the local and global state vector dimension.
!
! Debug output for this routine is turned on by preprocessor flag
! `PDAF_DEBUG`.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_model, &             ! Model variables
!        ONLY: nx, ny, nx_p
  
  USE mod_parallel_pdaf, &     ! Parallelization variables for
    ! assimilation
        ONLY: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe, &
        abort_parallel, &
        mype_world, COMM_model, npes_model, &
        mpi_success, mpi_comm_world, mpi_integer, mype_model
  USE mod_assimilation, &      ! Variables for assimilation
        ONLY: dim_state_p, dim_state, screen, filtertype, subtype, toffset,&
        dim_ens, rms_obs, model_error, model_err_amp, incremental, &
        type_forget, forget, dim_bias, rank_analysis_enkf, &
        locweight, cradius, sradius, filename, &
        type_trans, type_sqrt, delt_obs, toffset, dim_state_p_count, &
        dim_lag, &
        type_winf, limit_winf, &
        type_hyb, hyb_gamma, hyb_kappa, &
        pf_res_type, pf_noise_type, pf_noise_amp
  USE mod_tsmp, &
        ONLY: pf_statevecsize, nprocpf, tag_model_parflow, tag_model_clm, nprocclm, pf_statevec, pf_statevec_fortran, &
        idx_map_subvec2state, idx_map_subvec2state_fortran, model
#if defined CLMSA
  ! kuw: get access to clm variables
#ifndef CLMFIVE    
  USE shr_kind_mod , only : r8 => shr_kind_r8
  USE clm_atmlnd   , only : clm_l2a, atm_l2a, clm_mapl2a
  USE clmtype      , only : clm3, nameg
  USE subgridAveMod, only : p2g, c2g
  USE domainMod    , only : latlon_type
  USE clm_varpar   , only : nlevsoi
  USE decompMod    , only : get_proc_global, get_proc_bounds, adecomp
  USE spmdGathScatMod , only : gather_data_to_master
  USE spmdMod      , only : masterproc
#endif
  use enkf_clm_mod, only: clm_statevecsize
#endif
  ! kuw end

  use, intrinsic :: iso_c_binding

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
! Calls: PDAF_get_state
! Calls: PDAF_set_debug_flag
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation
  integer :: ierror

  ! External subroutines
  EXTERNAL :: init_ens         ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time,
                                ! and dimension of next observation
    distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
    prepoststep_ens_pdaf            ! User supplied pre/poststep routine
  integer :: i, j
  !kuw: dimensions for clm variables
  integer :: numg           ! total number of gridcells across all processors
  integer :: numl           ! total number of landunits across all processors
  integer :: numc           ! total number of columns across all processors
  integer :: nump           ! total number of pfts across all processors
  integer :: begg,endg      ! local beg/end gridcells gdc
  integer :: begl,endl      ! local beg/end landunits
  integer :: begc,endc      ! local beg/end columns
  integer :: begp,endp      ! local beg/end pfts
  !kuw end
  CHARACTER (LEN = 30) :: fn    !TSMP-PDAF: function name for output of index array

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

! *** Pointer initialization for ParFlow-type state vector ***
  if (model == tag_model_parflow) then
    ! Parflow: Initialize Fortran-pointer on pf_statevec
    call C_F_POINTER(pf_statevec, pf_statevec_fortran, [pf_statevecsize])

    ! Parflow: Initialize Fortran-pointer on idx_mapping_subvec2state
    call C_F_POINTER(idx_map_subvec2state, idx_map_subvec2state_fortran, [pf_statevecsize])
  end if

#ifdef PDAF_DEBUG
  ! Index array only needed for a single modeltask
  IF(filterpe) THEN

    ! idx_map_subvec2state is only set for ParFlow processes
    IF(model == tag_model_parflow) THEN
      ! Debug output: index manipulation array subvec->state
      WRITE(fn, "(a,i5.5,a)") "idx_map_subvec2state_", mype_world, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, pf_statevecsize
        WRITE (71,"(i10)") idx_map_subvec2state_fortran(i)
      END DO
      CLOSE(71)
    END IF

  END IF
#endif

! *** Define state dimension ***

! *** Setting local state vector dimension ***
  if (model == tag_model_parflow) then
    ! Parflow component, setting local state dimension `dim_state_p`
    ! and later dim_state from `pf_statevecsize` from `initialize_tsmp
    ! -> parflow_oasis_init`.
    dim_state_p = pf_statevecsize  ! Local state dimension
  else
    ! CLM/COSMO component, setting dummy dim_state_p and dim_state
    dim_state_p = 1  ! Local state dimension
  end if

#if defined CLMSA
  if (model == tag_model_clm) then

    ! CLM component: setting local state dimension from
    ! `clm_statevecsize` from `initialize_tsmp -> clm(5)_init ->
    ! define_clm_statevec`
    dim_state_p = clm_statevecsize

  end if
#endif

! *** Setting global state vector dimension ***
  IF (allocated(dim_state_p_count)) deallocate(dim_state_p_count)
  allocate(dim_state_p_count(npes_model))
  call MPI_Gather(dim_state_p, 1, MPI_INTEGER, dim_state_p_count, 1, MPI_INTEGER, 0, COMM_model, ierror)

#ifdef PDAF_DEBUG
  ! Debug output: local state dimension array
  if (mype_model == 0) WRITE(*, '(a,x,a,i5,x,a,x)', advance="no") "TSMP-PDAF-debug", "mype(w)=", mype_world, "init_pdaf: dim_state_p_count in modified:"
  if (mype_model == 0) WRITE(*, *) dim_state_p_count
#endif

  if (mype_model == 0) then
    dim_state = sum(dim_state_p_count)
  end if
  call MPI_BCAST(dim_state, 1, MPI_INTEGER, 0, COMM_model, IERROR)

#ifdef PDAF_DEBUG
  ! Debug output: global state dimension
  WRITE(*, '(a,x,a,i5,x,a,x,i9)') "TSMP-PDAF-debug", "mype(w)=", mype_world, "init_pdaf: my local state vector dimension dim_state_p:", dim_state_p
  WRITE(*, '(a,x,a,i5,x,a,2x,i9)') "TSMP-PDAF-debug", "mype(w)=", mype_world, "init_pdaf: my global state vector dimension dim_state:", dim_state
#endif

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 2    ! Write screen output (1) for output, (2) add timings

! *** Ensemble size ***
  dim_ens = n_modeltasks        ! Size of ensemble for all ensemble filters

! *** Options for filter method

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see MOD_ASSIMILATION +++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++

  filtertype = 2     ! Type of filter
  subtype = 0        ! Subtype of filter

  forget  = 1.0      ! Forgetting factor value for inflation
  type_forget = 0    ! Type of forgetting factor

  type_trans = 0     ! Type of ensemble transformation (deterministic or random)
  type_sqrt = 0      ! SEIK/LSEIK/ESTKF/LESTKF: Type of transform matrix square-root
  incremental = 0    ! SEIK/LSEIK: (1) to perform incremental updating

  !EnKF
  rank_analysis_enkf = 0  ! EnKF: rank to be considered for inversion of HPH in analysis step

  ! NETF/LNETF/PF
  type_winf = 0      ! NETF/LNETF/PF: Type of weights inflation
  limit_winf = 0.0   ! NETF/LNETF/PF: Limit for weights inflation

  ! LKNETF
  type_hyb = 0       ! LKNETF: Type of hybrid weight
  hyb_gamma =  1.0   ! LKNETF: Hybrid filter weight for state (1.0: LETKF, 0.0: LNETF)
  hyb_kappa = 30.0   ! LKNETF: Hybrid norm for using skewness and kurtosis (type_hyb 3 or 4)

  ! PF
  pf_res_type = 1    ! PF: Resampling type for particle filter
  pf_noise_type = 0  ! PF: Type of pertubing noise
  pf_noise_amp = 0.0 ! PF: Noise amplitude for particle filter


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 2     ! Number of time steps between analysis/assimilation steps

! *** specifications for observations ***
  rms_obs = 0.5    ! Observation error standard deviation

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  cradius = 0.0     ! Cut-off radius for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting

! *** File names
  filename = 'output.dat'

! *** TSMP-PDAF-specific inputs
  !kuw: add smoother support
  dim_lag = 0
  !kuw end

  ! hcp
  toffset = 0      ! offset of time steps shifting all analysis/assimilation steps
  ! hcp end
  
! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()

! *** Switch on debug output ***
! *** for main process        ***
#ifdef PDAF_DEBUG
  IF (mype_world == 0) CALL PDAF_set_debug_flag(1)
#endif

! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  whichinit: IF (filtertype == 2 .or. filtertype==8) THEN
    ! *** EnKF with Monte Carlo init ***
    filter_param_i(1) = dim_state_p ! State dimension
    filter_param_i(2) = dim_ens     ! Size of ensemble
    filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
    filter_param_i(4) = incremental ! Whether to perform incremental analysis
    filter_param_i(5) = 0           ! Smoother lag (not implemented here)
    !kuw: add smoother support
    filter_param_i(5) = dim_lag     ! Smoother lag (not implemented here)
    !kuw end
    filter_param_r(1) = forget      ! Forgetting factor

    !hcp 0-> toffset
    CALL PDAF_init(filtertype, subtype, toffset, &
      filter_param_i, 6,&
      filter_param_r, 2, &
      COMM_model, COMM_filter, COMM_couple, &
      task_id, n_modeltasks, filterpe, init_ens, &
      screen, status_pdaf)
  ELSE
    ! *** All other filters                       ***
    ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
    filter_param_i(1) = dim_state_p ! State dimension
    filter_param_i(2) = dim_ens     ! Size of ensemble
    filter_param_i(3) = 0           ! Smoother lag (not implemented here)
    filter_param_i(4) = incremental ! Whether to perform incremental analysis
    filter_param_i(5) = type_forget ! Type of forgetting factor
    filter_param_i(6) = type_trans  ! Type of ensemble transformation
    filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
    filter_param_r(1) = forget      ! Forgetting factor

    !hcp 0-> toffset
    CALL PDAF_init(filtertype, subtype, toffset, &
      filter_param_i, 7,&
      filter_param_r, 2, &
      COMM_model, COMM_filter, COMM_couple, &
      task_id, n_modeltasks, filterpe, init_ens, &
      screen, status_pdaf)
  END IF whichinit


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

! *** Switch off debug output ***
#ifdef PDAF_DEBUG
  CALL PDAF_set_debug_flag(0)
#endif

  if (mype_world == 0 .and. screen > 2) then
    print *, "TSMP-PDAF INITIALIZE PDAF FINISHED"
  end if

END SUBROUTINE init_pdaf
