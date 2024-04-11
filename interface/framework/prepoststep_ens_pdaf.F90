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
!prepoststep_ens_pdaf.F90: TSMP-PDAF implementation of routine
!                          'prepoststep_ens_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: prepoststep_ens_pdaf.F90 1444 2013-10-04 10:54:08Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
    state_p, Uinv, ens_p, flag)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! The routine is called for global filters (e.g. SEIK)
    ! before the analysis and after the ensemble transformation.
    ! For local filters (e.g. LSEIK) the routine is called
    ! before and after the loop over all local analysis
    ! domains.
    ! The routine provides full access to the state
    ! estimate and the state ensemble to the user.
    ! Thus, user-controlled pre- and poststep
    ! operations can be performed here. For example
    ! the forecast and the analysis states and ensemble
    ! covariance matrix can be analized, e.g. by
    ! computing the estimated variances.
    ! For the offline mode, this routine is the place
    ! in which the writing of the analysis ensemble
    ! can be performed.
    !
    ! If a user considers to perform adjustments to the
    ! estimates (e.g. for balances), this routine is
    ! the right place for it.
    !
    ! !REVISION HISTORY:
    ! 2013-09 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_assimilation, &
        ONLY: dim_state, dim_state_p_count
    USE mod_parallel_pdaf, &
        ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
        MPIerr, MPIstatus, filterpe, mype_model, npes_model, mype_world, &
        MPI_COMM_WORLD, MPI_SUCCESS
    use mod_tsmp, &
        only: tag_model_parflow, pf_statevecsize, nprocclm, model


    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: step        ! Current time step (negative for call after forecast)
    INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
    INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
    INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
    REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
    ! The array 'state_p' is not generally not initialized in the case of SEIK.
    ! It can be used freely here.
    REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
    INTEGER, INTENT(in) :: flag        ! PDAF status flag

    ! !CALLING SEQUENCE:
    ! Called by: PDAF_get_state      (as U_prepoststep)
    ! Called by: PDAF_X_update       (as U_prepoststep)
    ! Calls: MPI_send
    ! Calls: MPI_recv
    !EOP

    ! *** local variables ***
    INTEGER :: i, j, member, domain     ! Counters
    INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
    LOGICAL, SAVE :: firstio = .TRUE.   ! File output is peformed for first time?
    LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
    REAL :: invdim_ens                  ! Inverse ensemble size
    REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
    REAL :: rmserror_est                ! estimated RMS error
    REAL, ALLOCATABLE :: variance_p(:)  ! model state variances
    !REAL, ALLOCATABLE :: field(:,:)     ! global model field
    CHARACTER(len=2) :: ensstr          ! String for ensemble member
    CHARACTER(len=2) :: stepstr         ! String for time step
    CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
    ! Variables for parallelization - global fields
    INTEGER :: offset   ! Row-offset according to domain decomposition
    INTEGER, ALLOCATABLE :: dim_state_p_stride(:) ! local state vector sizes summation array
    REAL, ALLOCATABLE :: variance(:)    ! local variance
    !REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
    !REAL, ALLOCATABLE :: state(:)       ! global state vector
    !REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
    !REAL,ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain
    integer :: ierror

    ! **********************
    ! *** INITIALIZATION ***
    ! **********************
    if (2 .eq. 1) then
    IF (mype_filter == 0) THEN
        IF (firsttime) THEN
            WRITE (*, '(8x, a)') 'Analize initial state ensemble'
            anastr = 'ini'
        ELSE
            IF (step<0) THEN
                WRITE (*, '(8x, a)') 'Analize and write forecasted state ensemble'
                anastr = 'for'
            ELSE
                WRITE (*, '(8x, a)') 'Analize and write assimilated state ensemble'
                anastr = 'ana'
            END IF
        END IF
    END IF
    ! Allocate fields
    ALLOCATE(variance_p(dim_p))
    ALLOCATE(variance(dim_state))

    ! Initialize numbers
    rmserror_est  = 0.0
    invdim_ens    = 1.0 / REAL(dim_ens)
    invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


    ! **************************************************************
    ! *** Perform prepoststep for SEIK with re-inititialization. ***
    ! *** The state and error information is completely in the   ***
    ! *** ensemble.                                              ***
    ! *** Also performed for SEIK without re-init at the initial ***
    ! *** time.                                                  ***
    ! **************************************************************

    ! *** Compute mean state
    IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'

    !    state_p = 0.0
    state_p = 0.0
    DO member = 1, dim_ens
        DO i = 1, dim_p
            state_p(i) = state_p(i) + ens_p(i, member)
        END DO
    END DO
    state_p(:) = invdim_ens * state_p(:)
    print *, "mean state_p is:", state_p(1)

    ! *** Compute sampled variances ***
    variance_p(:) = 0.0
    DO member = 1, dim_ens
        DO j = 1, dim_p
            variance_p(j) = variance_p(j) &
                + (ens_p(j, member) - state_p(j)) &
                * (ens_p(j, member) - state_p(j))
        END DO
    END DO
    variance_p(:) = invdim_ensm1 * variance_p(:)

    ! ******************************************************
    ! *** Assemble global variance vector on filter PE 0 ***
    ! ******************************************************

    WRITE (*,*) 'TEMPLATE prepoststep_ens_pdaf.F90: Initialize variance, either directly or with MPI'
    if (filterpe) then
        call MPI_Barrier(comm_filter, ierror)
    end if

    IF (allocated(dim_state_p_stride)) deallocate(dim_state_p_stride)
    allocate(dim_state_p_stride(npes_model))
    do i = 1, npes_model
      dim_state_p_stride(i) = 0
      do j = 1, i - 1
        dim_state_p_stride(i) = dim_state_p_count(j) + dim_state_p_stride(i)
      end do
    end do
#ifdef PDAF_DEBUG
    ! Debug output: summed until index local state dimension array
    if (mype_model == 0 ) WRITE(*, '(a,x,a,i5,x,a,x,i9)') "TSMP-PDAF-debug", "mype(w)=", mype_world, "init_pdaf: dim_state_p_stride in modified:", dim_state_p_stride
#endif

    !!!!!!!!!!!!!!!!!!!!!!!! case below contains dummy CLM component  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (filterpe) then
        print *, "prepoststep: gathering variance"
        !        call MPI_Gather(variance_p, dim_p, MPI_DOUBLE_PRECISION, &
        !            variance, dim_p, MPI_DOUBLE_PRECISION, &
        !            0, comm_filter, MPIerr)
    call MPI_Gatherv(variance_p, dim_p, MPI_DOUBLE_PRECISION, variance, dim_state_p_count, &
        dim_state_p_stride, MPI_DOUBLE_PRECISION, 0, comm_filter, MPIerr);
        if (MPIerr .ne. MPI_SUCCESS) then
            print *,"mpi gather failed"
            call MPI_Abort(MPI_COMM_WORLD, MPIerr)
        end if
        print *, "prepoststep: gathering variance succeeded"
    end if
    DEALLOCATE(variance_p)
    IF (allocated(dim_state_p_stride)) deallocate(dim_state_p_stride)


    ! ************************************************************
    ! *** Compute RMS errors according to sampled covar matrix ***
    ! ************************************************************

    ! total estimated RMS error
    IF (mype_filter == 0) THEN
        DO i = 1, dim_state
            rmserror_est = rmserror_est + variance(i)
        ENDDO
        rmserror_est = SQRT(rmserror_est / dim_state)
    END IF
    DEALLOCATE(variance)


    ! *****************
    ! *** Screen IO ***
    ! *****************

    ! Output RMS errors given by sampled covar matrix
    !    if (model == tag_model_parflow) then
    IF (mype_filter == 0) THEN
        WRITE (*, '(12x, a, es12.4)') &
            'RMS error according to sampled variance: ', rmserror_est
    END IF
    !    end if

    ! *******************
    ! *** File output ***
    ! *******************
    notfirst: IF (.not. firsttime) THEN
        WRITE (*,*) 'TEMPLATE prepoststep_ens_pdaf.F90: Implement writing of output files here!'
    END IF notfirst


    ! for testing purpose, print the state vector
    print *, "prepoststep: first element in state vector: ", state_p(1)
    ! ********************
    ! *** finishing up ***
    ! ********************

    firsttime = .FALSE.
    end if
END SUBROUTINE prepoststep_ens_pdaf
