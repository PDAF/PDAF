!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
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
!init_pdaf.F90: TerrSysMP-PDAF implementation of routine
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
    ! !REVISION HISTORY:
    ! 2008-10 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    !   USE mod_model, &             ! Model variables
    !        ONLY: nx, ny, nx_p

    USE mod_parallel_model, &    ! Parallelization variables for model
        ONLY: mype_world, COMM_model, abort_parallel,  model, npes_model, &
        mpi_success, mpi_comm_world, mpi_integer, mype_model
    use mod_tsmp, &
        only: pf_statevecsize, nprocpf, tag_model_parflow, tag_model_clm, nprocclm, pf_statevec, pf_statevec_fortran, &
        idx_map_subvec2state, idx_map_subvec2state_fortran
    USE mod_parallel_pdaf, &     ! Parallelization variables fro assimilation
        ONLY: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe
    USE mod_assimilation, &      ! Variables for assimilation
        ONLY: dim_state_p, dim_state, screen, filtertype, subtype, toffset,&
        dim_ens, rms_obs, model_error, model_err_amp, incremental, &
        covartype, type_forget, forget, dim_bias, rank_analysis_enkf, &
        locweight, local_range, srange, int_rediag, filename, &
        type_trans, type_sqrt, delt_obs, toffset, dim_state_p_count, dim_state_p_stride,&
        dim_lag
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
    !include "mpif.h"
    ! !CALLING SEQUENCE:
    ! Called by: main
    ! Calls: init_pdaf_parse
    ! Calls: init_pdaf_info
    ! Calls: PDAF_init
    ! Calls: PDAF_get_state
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

    ! ***************************
    ! ***   Initialize PDAF   ***
    ! ***************************

    IF (mype_world == 0) THEN
        WRITE (*,'(a)') 'TSMP-PDAF INITIALIZE PDAF - ONLINE MODE'
    END IF

    ! **********************************************************
    ! ***               CONTROL OF PDAF                      ***
    ! **********************************************************

    ! *** Forecast length (time interval between analysis steps) ***
    delt_obs = 2     ! Number of time steps between analysis/assimilation steps
    toffset = 0     ! offset of time steps shifting all analysis/assimilation steps

    ! *** IO options ***
    screen      = 2  ! Write screen output (1) for output, (2) add timings

    ! *** specifications for observations ***
    ! avg. observation error (used for assimilation)
    rms_obs = 0.5    ! This error is the standard deviation
    ! for the Gaussian distribution

    ! *** Filter specific variables
    filtertype = 2    ! Type of filter
    !   (1) SEIK
    !   (2) EnKF
    !   (3) LSEIK
    !   (4) ETKF
    !   (5) LETKF
    !   (6) ESTKF
    !   (7) LESTKF
    dim_ens = n_modeltasks       ! Size of ensemble for all ensemble filters
    ! Number of EOFs to be used for SEEK
    subtype = 0       ! subtype of filter:
    !   ESTKF:
    !     (0) Standard form of ESTKF
    !   LESTKF:
    !     (0) Standard form of LESTKF
    type_trans = 0    ! Type of ensemble transformation
    !   SEIK/LSEIK and ESTKF/LESTKF:
    !     (0) use deterministic omega
    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
    !     (2) use product of (0) with random orthonormal matrix with
    !         eigenvector (1,...,1)^T
    !   ETKF/LETKF:
    !     (0) use deterministic symmetric transformation
    !     (2) use product of (0) with random orthonormal matrix with
    !         eigenvector (1,...,1)^T
    type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !   (0) fixed
    !   (1) global adaptive
    !   (2) local adaptive for LSEIK/LETKF/LESTKF
    forget  = 1.0     ! Forgetting factor
    type_sqrt = 0     ! Type of transform matrix square-root
    !   (0) symmetric square root, (1) Cholesky decomposition
    incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
    covartype = 1     ! Definition of factor in covar. matrix used in SEIK
    !   (0) for dim_ens^-1 (old SEIK)
    !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
    !   This parameter has also to be set internally in PDAF_init.
    rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition
    int_rediag = 1    ! Interval of analysis steps to perform
    !    re-diagonalization in SEEK
    locweight = 0     ! Type of localizating weighting
    !   (0) constant weight of 1
    !   (1) exponentially decreasing with SRANGE
    !   (2) use 5th-order polynomial
    !   (3) regulated localization of R with mean error variance
    !   (4) regulated localization of R with single-point error variance
    local_range = 0  ! Range in grid points for observation domain in local filters
    srange = local_range  ! Support range for 5th-order polynomial
    ! or range for 1/e for exponential weighting

    ! *** File names
    filename = 'output.dat'

    !kuw: add smoother support
    dim_lag = 0
    !kuw end

    ! ***********************************
    ! *** Some optional functionality ***
    ! ***********************************

    ! *** Parse command line options   ***
    ! *** This is optional, but useful ***

    call init_pdaf_parse()


    ! *** Initial Screen output ***
    ! *** This is optional      ***

    IF (mype_world == 0) call init_pdaf_info()


    ! *** Define state dimension ***

    if (model == tag_model_parflow) then
        if (screen > 2) then
            print *, "Parflow: converting pf_statevec to fortran"
        end if

        call C_F_POINTER(pf_statevec, pf_statevec_fortran, [pf_statevecsize])

        if (screen > 2) then
            print *, "Parflow: converting idx_mapping_subvec2state to fortran"
        end if

        call C_F_POINTER(idx_map_subvec2state, idx_map_subvec2state_fortran, [pf_statevecsize])

        ! if (screen > 2) then
        !     print *, "Parflow: first several elements of the idx:", idx_map_subvec2state_fortran
        ! end if
    end if

    if (model == tag_model_parflow) then
        dim_state_p = pf_statevecsize  ! Local state dimension

        if (screen > 2) then
            print *,""
            print *, "TSMP-PDAF mype(w)=", mype_world, ": Parflow component, setting correct dim_state_p and dim_state"
        end if
    else
        if (screen > 2) then
            print *,""
            print *, "TSMP-PDAF mype(w)=", mype_world, ": CLM component, setting dummy dim_state_p and dim_state"
        end if

        dim_state_p = 1  ! Local state dimension
    end if

#if defined CLMSA
    if (model == tag_model_clm) then
       ! comment only CLMSA
       !call get_proc_global(numg,numl,numc,nump)
       !call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
       !dim_state_p =  (endg-begg+1) * nlevsoi

        dim_state_p = clm_statevecsize

        if (screen > 2) then
            print *,"TSMP-PDAF mype(w)=", mype_world, ": CLM: dim_state_p is ",dim_state_p
        end if
    end if
#endif

    IF (allocated(dim_state_p_count)) deallocate(dim_state_p_count)
    allocate(dim_state_p_count(npes_model))
    call MPI_Gather(dim_state_p, 1, MPI_INTEGER, dim_state_p_count, 1, MPI_INTEGER, 0, comm_model, ierror)

    if (mype_model == 0 .and. screen > 2) print *, "TSMP-PDAF mype(w)=", mype_world, ": init_pdaf: dim_state_p_count in modified: ", dim_state_p_count
    IF (allocated(dim_state_p_stride)) deallocate(dim_state_p_stride)
    allocate(dim_state_p_stride(npes_model))
    do i = 1, npes_model
        dim_state_p_stride(i) = 0
        do j = 1, i - 1
            dim_state_p_stride(i) = dim_state_p_count(j) + dim_state_p_stride(i)
        end do
    end do
    if (mype_model == 0 .and. screen > 2) print *, "TSMP-PDAF mype(w)=", mype_world, ": init_pdaf: dim_state_p_stride in modified: ", dim_state_p_stride

    if (mype_model == 0) then
        dim_state = sum(dim_state_p_count)
    end if
    call MPI_BCAST(dim_state, 1, MPI_INTEGER, 0, comm_model, IERROR)
    !print  *, "my local state vector dimension is :" , dim_state_p
    !print  *, "my global state vector dimension is :" , dim_state
    !print *,""
    call MPI_Barrier(MPI_COMM_WORLD, ierror)

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

    if (mype_world == 0 .and. screen > 2) then
        print *, "TSMP-PDAF INITIALIZE PDAF FINISHED"
    end if

END SUBROUTINE init_pdaf
