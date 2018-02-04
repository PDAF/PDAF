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
!$Id: PDAF-D_lestkf_update.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lestkf_update --- Control analysis update of the LESTKF
!
! !INTERFACE:
SUBROUTINE  PDAF_lestkf_update(step, dim_p, dim_obs_f, dim_ens, rank, &
     state_p, Ainv, ens_p, state_inc_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prodRinvA_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, U_prepoststep, screen, &
     subtype, incremental, type_forget, type_sqrt, dim_lag, &
     sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the LESTKF filter.
!
! The analysis is performed by first preparing several
! global quantities on the PE-local domain, like the
! observed part of the state ensemble for all local
! analysis domains on the PE-local state domain.
! Then the analysis and ensemble tranformations are 
! performed within a loop over all local analysis domains
! in the PE-local state domain in the subroutine
! (PDAF\_lestkf\_analysis). In this loop, the local state and 
! observation dimensions are initialized and the global 
! state ensemble is restricted to the local analysis domain.
! In addition, the routine U\_prepoststep is called prior
! to the analysis and after the resampling outside of
! the loop over the local domains to allow the user
! to access the ensemble information.
!
! Variant for domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: type_trans, filterstr, obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l, npes_filter, COMM_filter, MPIerr, &
       MPI_SUM, MPI_MAX, MPI_INTEGER

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
! !   suffix _f: Denotes a full variable of all observations required for the
! !              analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  INTEGER, INTENT(in) :: rank        ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Ainv(rank, rank)      ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: subtype     ! Filter subtype
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
  INTEGER, INTENT(in) :: dim_lag     ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_obs_op, &    ! Observation operator
       U_init_n_domains_p, & ! Provide number of local analysis domains
       U_init_dim_l, &       ! Init state dimension for local ana. domain
       U_init_dim_obs, &     ! Initialize dimension of observation vector
       U_init_dim_obs_l, &   ! Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &         ! Initialize observation vector
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_obsvar, &      ! Initialize mean observation error variance
       U_init_obsvar_l, &    ! Initialize local mean observation error variance
       U_g2l_state, &        ! Get state on local ana. domain from global state
       U_l2g_state, &        ! Init full state from state on local analysis domain
       U_g2l_obs, &          ! Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &      ! Compute product of R^(-1) with HV
       U_prepoststep         ! User supplied pre/poststep routine

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_lestkf
! Calls: U_prepoststep
! Calls: U_init_n_domains_p
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_init_dim_l
! Calls: U_init_dim_obs_l
! Calls: U_g2l_state
! Calls: U_l2g_state
! Calls: PDAF_set_forget
! Calls: PDAF_lestkf_analysis
! Calls: PDAF_seik_omega
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: MPI_reduce
!EOP

! *** local variables ***
  INTEGER :: i, j, member, row     ! Counters
  INTEGER :: domain_p              ! Counter for local analysis domain
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  REAL    :: invdimens             ! Inverse global ensemble size
  INTEGER :: minusStep             ! Time step counter
  INTEGER :: n_domains_p           ! number of PE-local analysis domains
  REAL    :: forget_ana_l          ! forgetting factor supplied to analysis routine
  REAL    :: forget_ana            ! Possibly globally adaptive forgetting factor
  LOGICAL :: storeOmega = .FALSE.  ! Store matrix Omega instead of recomputing it
  REAL, ALLOCATABLE :: HX_f(:,:)   ! HX for PE-local ensemble
  REAL, ALLOCATABLE :: HXbar_f(:)  ! PE-local observed mean state
  REAL, ALLOCATABLE :: obs_f(:)    ! PE-local observation vector
  REAL, ALLOCATABLE :: Omega(:,:)  ! Transformation matrix Omega
  REAL, ALLOCATABLE :: OmegaT(:,:) ! Transpose of transformation matrix Omeg
  REAL, SAVE, ALLOCATABLE :: OmegaT_save(:,:) ! Stored OmegaT
  ! Variables on local analysis domain
  INTEGER :: dim_l                ! State dimension on local analysis domain
  INTEGER :: dim_obs_l            ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: ens_l(:,:) ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:) ! Mean state on local analysis domain
  REAL, ALLOCATABLE :: stateinc_l(:)  ! State increment on local analysis domain
  REAL, ALLOCATABLE :: TA_l(:,:)  ! Local ensemble transform matrix
  ! Variables for statistical information on local analysis
  INTEGER :: obsstats(4)           ! PE-local statistics
  INTEGER :: obsstats_g(4)         ! Global statistics
  ! obsstats(1): Local domains with observations
  ! obsstats(2): Local domains without observations
  ! obsstats(3): Sum of all available observations for all domains
  ! obsstats(4): Maximum number of observations over all domains
  REAL, ALLOCATABLE :: Ainv_l(:,:)  ! thread-local matrix Ainv


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  CALL PDAF_timeit(5, 'new')
  minusStep = - step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF


! **************************************
! *** Preparation for local analysis ***
! **************************************

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(4, 'new')

  ! Query number of analysis domains for the local analysis
  ! in the PE-local domain
  CALL U_init_n_domains_p(step, n_domains_p)
  
  IF (screen > 0) THEN
     IF (mype == 0) THEN
        IF (subtype /= 3) THEN
           WRITE (*, '(a, i7, 3x, a)') 'PDAF ', step, 'Local ESTKF analysis'
        ELSE
           WRITE (*, '(a, i7, 3x, a)') 'PDAF ', step, 'LESTKF analysis for fixed covariance matrix'
        END IF
     END IF
     WRITE (*, '(a, 5x, a, i6, a, i10)') &
          'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
  END IF


! *** Local analysis: initialize global quantities ***

  ! Get observation dimension for all observations required 
  ! for the loop of local analyses on the PE-local domain.
  CALL U_init_dim_obs(step, dim_obs_f)

  IF (screen > 0) THEN
     WRITE (*, '(a, 5x, a, i6, a, i10)') &
          'PDAF', '--- PE-Domain:', mype, &
          ' dimension of PE-local full obs. vector', &
          dim_obs_f
  END IF

  ! HX = [Hx_1 ... Hx_(r+1)] for full DIM_OBS_F region on PE-local domain
  ALLOCATE(HX_f(dim_obs_f, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f * dim_ens)

  ENS: DO member = 1,dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! Call observation operator
     CALL U_obs_op(step, dim_p, dim_obs_f, ens_p(:, member), HX_f(:, member))
  END DO ENS

  CALL PDAF_timeit(11, 'new')

  ! *** Compute mean forecast state
  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')

  ! *** Compute mean state of ensemble on PE-local observation space 
  ALLOCATE(HXbar_f(dim_obs_f))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f)

  HXbar_f = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_obs_f
        HXbar_f(row) = HXbar_f(row) + invdimens * HX_f(row, member)
     END DO
  END DO

  ! Set forgetting factor globally
  forget_ana = forget
  IF (type_forget == 1) THEN
     ALLOCATE(obs_f(dim_obs_f))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f)

     ! get observation vector
     CALL U_init_obs(step, dim_obs_f, obs_f)

     ! Set FORGET
     CALL PDAF_set_forget(step, filterstr, dim_obs_f, dim_ens, HX_f, &
          HXbar_f, obs_f, U_init_obsvar, forget, forget_ana)
     
     DEALLOCATE(obs_f)
  ENDIF

  ! *** Initialize OmegaT
  CALL PDAF_timeit(33, 'new')
  ALLOCATE(omega(dim_ens, rank))
  ALLOCATE(omegaT(rank, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank * dim_ens)

  O_store: IF (.NOT. storeOmega .OR. (storeOmega .AND. allocflag == 0)) THEN

     CALL PDAF_seik_omega(rank, Omega, type_trans, screen)
     omegaT = TRANSPOSE(Omega)

     IF (storeOmega) THEN
        ALLOCATE(omegaT_save(rank, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank * dim_ens)
        omegaT_save = omegaT
     END IF

  ELSE O_store
     ! Re-use stored Omega
     if (mype == 0 .AND. screen > 0) &
          write (*,'(a, 5x, a)') 'PDAF', '--- Use stored Omega'
     omegaT = omegaT_save
  END IF O_store

  DEALLOCATE(omega)
  CALL PDAF_timeit(33, 'old')

  CALL PDAF_timeit(4, 'old')


! ************************************
! *** Perform analysis and re_init ***
! ************************************


  ! Initialize counters for statistics on local observations
  obsstats = 0

  CALL PDAF_timeit(6, 'new')

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, stateinc_l, TA_l, Ainv_l, flag) firstprivate(forget_ana_l)

  forget_ana_l = forget_ana

  ! Allocate ensemble transform matrix
  ALLOCATE(TA_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  TA_l = 0.0

  ! Allocate ensemble transform matrix
  ALLOCATE(Ainv_l(dim_ens-1, dim_ens-1))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', (dim_ens-1)**2)
  Ainv_l = 0.0

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag)
  localanalysis: DO domain_p = 1, n_domains_p

     ! local state dimension
     CALL U_init_dim_l(step, domain_p, dim_l)

     ! Get observation dimension for local domain
     CALL PDAF_timeit(9, 'new')
     CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
     CALL PDAF_timeit(9, 'old')

     ! Gather statistical information on local observations
!$OMP CRITICAL
     IF (dim_obs_l > obsstats(4)) obsstats(4) = dim_obs_l
     IF (dim_obs_l > 0) THEN
        obsstats(3) = obsstats(3) + dim_obs_l
        obsstats(1) = obsstats(1) + 1
     ELSE
        obsstats(2) = obsstats(2) + 1
     END IF
!$OMP END CRITICAL
     
     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     ALLOCATE(stateinc_l(dim_l))

     CALL PDAF_timeit(15, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))
     END DO
     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     CALL PDAF_timeit(15, 'old')

     ! Check whether we have observations for the current local domain
     ! Perform analysis only if we have observations
     CALL PDAF_timeit(7, 'new')

     havelocalobs: IF (dim_obs_l > 0) THEN
        IF (subtype /= 3) THEN
           ! LESTKF analysis for current domain
           CALL PDAF_lestkf_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
                dim_ens, rank, state_l, Ainv_l, ens_l, HX_f, &
                HXbar_f, stateinc_l, OmegaT, forget_ana_l, U_g2l_obs, &
                U_init_obs_l, U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, screen, &
                incremental, type_forget, type_sqrt, TA_l, flag)
        ELSE
           ! LESTKF analysis with state update but no ensemble transformation
           CALL PDAF_lestkf_analysis_fixed(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
                dim_ens, rank, state_l, Ainv_l, ens_l, HX_f, &
                HXbar_f, stateinc_l, OmegaT, forget_ana_l, U_g2l_obs, &
                U_init_obs_l, U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, screen, &
                incremental, type_forget, type_sqrt, flag)
        END IF
     END IF havelocalobs

     CALL PDAF_timeit(7, 'old')
    
     CALL PDAF_timeit(16, 'new')

     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
     IF (subtype /= 4) THEN
        CALL U_l2g_state(step, domain_p, dim_l, state_l, dim_p, state_p)
     END IF
    
     ! Initialize global state increment
     IF (incremental == 1) THEN
        CALL U_l2g_state(step, domain_p, dim_l, stateinc_l, dim_p, state_inc_p)
     END IF

     CALL PDAF_timeit(16, 'old')
     CALL PDAF_timeit(17, 'new')

     ! *** Perform smoothing of past ensembles ***
     CALL PDAF_smoother_local(domain_p, step, dim_p, dim_l, dim_ens, &
          dim_lag, TA_l, ens_l, sens_p, cnt_maxlag, &
          U_g2l_state, U_l2g_state, forget_ana_l, screen)

     CALL PDAF_timeit(17, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l, stateinc_l)

  END DO localanalysis

!$OMP CRITICAL
  ! Set Ainv - required for subtype=3
  Ainv = Ainv_l
!$OMP END CRITICAL

  DEALLOCATE(TA_l, Ainv_l)
!$OMP END PARALLEL

  CALL PDAF_timeit(6, 'old')

  ! *** Print statistics for local analysis to the screen ***
  IF (npes_filter>1) THEN
     CALL MPI_Reduce(obsstats, obsstats_g, 3, MPI_INTEGER, MPI_SUM, &
          0, COMM_filter, MPIerr)
     CALL MPI_Reduce(obsstats(4), obsstats_g(4), 1, MPI_INTEGER, MPI_MAX, &
          0, COMM_filter, MPIerr)
  ELSE
     ! This is a work around for working with nullmpi.F90
     obsstats_g = obsstats
  END IF

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', '--- Global statistics for local analysis:'
     WRITE (*, '(a, 8x, a, i10)') &
          'PDAF', 'Local domains with observations:       ', obsstats_g(1)
     WRITE (*, '(a, 8x, a, i10)') &
          'PDAF', 'Local domains without observations:    ', obsstats_g(2)
     WRITE (*, '(a, 8x, a, i10)') &
          'PDAF', 'Maximum local observation dimension:   ', obsstats_g(4)
     WRITE (*, '(a, 8x, a, f9.1)') &
          'PDAF', 'Total avg. local observation dimension:', &
          REAL(obsstats_g(3)) / REAL(obsstats_g(1) + obsstats_g(2))
     IF (obsstats_g(2) > 0) THEN
        WRITE (*, '(a, 8x, a, f9.1)') &
             'PDAF', 'Avg. for domains with observations:    ', &
             REAL(obsstats_g(3)) / REAL(obsstats_g(1))
     END IF
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF

! *** Clean up from local analysis update ***
  DEALLOCATE(HX_f, HXbar_f, OmegaT)
#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

! *** Poststep for analysis ensemble ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_lestkf_update
