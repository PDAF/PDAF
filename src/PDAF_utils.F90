! Copyright (c) 2004-2026 Lars Nerger
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
!
!> Module collecting differrent utlility routines for PDAF
MODULE PDAF_utils
!> Utility routines for PDAF
!!
!! !  These are core routines of PDAF and
!!    should not be changed by the user   !
!!
CONTAINS

!-------------------------------------------------------------------------------
!> General allocation routine of PDAF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code from restructuring
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_alloc(dim_p, dim_ens, dim_ens_task, dim_es, &
     statetask, outflag)

  USE mpi
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: screen, state, ens, &
       Ainv, sens, bias
  USE PDAF_mod_parallel, &
       ONLY: mype, mype_model, mype_couple, filterpe, task_id, &
       COMM_couple

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< Size of state vector
  INTEGER, INTENT(in) :: dim_ens         !< Ensemble size
  INTEGER, INTENT(in) :: dim_ens_task    !< Ensemble size handled by a model task
  INTEGER, INTENT(in) :: dim_es          !< Dimension of error space (size of Ainv)
  INTEGER, INTENT(in) :: statetask       !< Task ID forecasting a single state
  INTEGER, INTENT(inout):: outflag       !< Status flag

! *** local variables ***
  INTEGER :: allocstat                   ! Status for allocate


! ****************************
! *** Allocate PDAF arrays ***
! ****************************

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_p)

     ! Allocate full ensemble on filter-PEs
     ALLOCATE(ens(dim_p, dim_ens), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_ens)
     
     IF (dim_es > 1) THEN
        ALLOCATE(Ainv(dim_es, dim_es), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of Ainv'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(1, 'r', dim_es**2)
     ELSE
        ALLOCATE(Ainv(1, 1), stat = allocstat)
     END IF

     IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate ens of size ', &
          dim_ens, ' on pe(f) ', mype

     ! SENS and BIAS have their own allocation routines.
     ! If they are not alled in PDAF_X_SET_IPARAM, we allocate
     ! these arrays to their minimum size.
     IF (.NOT.ALLOCATED(sens)) ALLOCATE(sens(1,1,1))
     IF (.NOT.ALLOCATED(bias)) ALLOCATE(bias(1))
     
  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        ALLOCATE(ens(dim_p, dim_ens_task), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of ens on model-pe'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_ens_task)

        IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate ens of size ', &
             dim_ens_task, ' on pe(m) ', mype_model, ' of model task ',task_id

        ! Some of the model-PEs may integrate a central state
        IF (mype_couple+1 == statetask) THEN
           ALLOCATE(state(dim_p), stat = allocstat)
           IF (allocstat /= 0) THEN
              WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
              outflag = 20
           END IF
              
           IF (screen > 2) WRITE (*,*) 'PDAF: alloc - allocate state of size ', &
                dim_p, ' on pe(m) ', mype_model, ' of model task ',task_id

           ! count allocated memory
           CALL PDAF_memcount(2, 'r', dim_p)
        END IF
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_alloc

!-------------------------------------------------------------------------------
!> Allocation routine of PDAF for smoother array
!!
!! The smoother array is allocated here separate from the
!! general allocation in PDAF_alloc. The routine is 
!! called from the PDAF_X_set_iparam of each filter method
!! that supports smoothing. The separation is needed to 
!! allow setting dim_lag using PDAF_set_iparam.
!!
!! __Revision history:__
!! * 2025-06 - Lars Nerger - Initial code from restructuring
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_alloc_sens(dim_p, dim_ens, dim_lag, outflag)

  USE mpi
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: sens
  USE PDAF_mod_parallel, &
       ONLY: filterpe 

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< Size of state vector
  INTEGER, INTENT(in) :: dim_ens         !< Ensemble size
  INTEGER, INTENT(in) :: dim_lag         !< Smoother lag
  INTEGER, INTENT(inout):: outflag       !< Status flag

! *** local variables ***
  INTEGER :: allocstat                   ! Status for allocate


! ****************************
! *** Allocate PDAF arrays ***
! ****************************

  on_filterpe: IF (filterpe) THEN
     
     IF (ALLOCATED(sens)) DEALLOCATE(sens)

     ! Allocate array for past ensembles for smoothing on filter-PEs
     IF (dim_lag > 0) THEN
        ALLOCATE(sens(dim_p, dim_ens, dim_lag), stat = allocstat)

        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_ens * dim_lag)
     ELSE
        ALLOCATE(sens(1, 1, 1), stat = allocstat)
     END IF
     IF (allocstat /= 0) THEN
        WRITE (*,*) 'PDAF-ERROR(20): error in allocation of sens'
        outflag = 20
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_alloc_sens

!-------------------------------------------------------------------------------
!> Allocation routine of PDAF for bias array
!!
!! The bias array is allocated here separate from the
!! general allocation in PDAF_alloc. The routine should be
!! called from the PDAF_X_set_iparam of each filter method
!! that supports a bias estimation. The separation is needed
!! to allow setting dim_bias using PDAF_set_iparam.
!!
!! __Revision history:__
!! * 2025-06 - Lars Nerger - Initial code from restructuring
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_alloc_bias(dim_bias_p, outflag)

  USE mpi
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: bias
  USE PDAF_mod_parallel, &
       ONLY: filterpe 

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_bias_p      !< Size of bias vector
  INTEGER, INTENT(inout):: outflag       !< Status flag

! *** local variables ***
  INTEGER :: allocstat                   ! Status for allocate


! ****************************
! *** Allocate PDAF arrays ***
! ****************************

  on_filterpe: IF (filterpe) THEN
     
     IF (ALLOCATED(bias)) DEALLOCATE(bias)

     IF (dim_bias_p > 0) THEN
        ALLOCATE(bias(dim_bias_p), stat = allocstat)
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_bias_p)
        
        ! initialize bias field
        bias = 0.0
     ELSE
        ALLOCATE(bias(1), stat = allocstat)
     ENDIF
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of BIAS'
        outflag = 20
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_alloc_bias


!-------------------------------------------------------------------------------
!> Deallocate PDAF-internal arrays
!!
!! Perform deallocation of PDAF-internal arrays
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2017-08 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_deallocate()

  USE mpi
  USE PDAF_mod_core, &
       ONLY: dim_bias_p, state, Ainv, ens, &
       sens, bias, dim_lag
  USE PDAF_mod_parallel, &
       ONLY: filterpe, COMM_couple, mpi_init_by_pdaf, MPIerr
  USE PDAF_iau, &
       ONLY: PDAF_iau_dealloc

  IMPLICIT NONE


! ******************************
! *** Deallocate PDAF fields ***
! ******************************

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full ensemble matrix on Filter-PEs

     DEALLOCATE(state)
     DEALLOCATE(Ainv)

     ! Allocate full ensemble on filter-PEs
     DEALLOCATE(ens)

     ! Allocate array for past ensembles for smoothing on filter-PEs
     IF (dim_lag > 0) THEN
        DEALLOCATE(sens)
     END IF

     IF (dim_bias_p > 0) THEN
        DEALLOCATE(bias)
     ENDIF

     ! Deallocate IAU arrays
     CALL PDAF_iau_dealloc()
     
  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local ensemble
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs that do coupling communication
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        DEALLOCATE(ens)
     END IF

  END IF on_filterpe

  ! Finalize MPI, if initialized by PDAF
  IF (mpi_init_by_pdaf) CALL MPI_finalize(MPIerr)

END SUBROUTINE PDAF_deallocate


!-------------------------------------------------------------------------------
!> Get value of a correlation function
!!
!! This routine returns the value of the chosen correlation
!! function according to the specified length scale.
!!
!!  This is a core routine of PDAF and
!!  should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_correlation_function(ctype, length, distance, value)

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: ctype          !< Type of correlation function
                                        !< (1): Gaussian with f(0)=1.0
                                        !< (2): 5th-order polynomial (Gaspace/Cohn, 1999)
  REAL, INTENT(in) :: length            !< Length scale of function
                                        !< (1): standard deviation
                                        !< (2): support length f=0 for distance>length 
  REAL, INTENT(in) :: distance          !< Distance at which the function is evaluated
  REAL, INTENT(out) :: value            !< Value of the function
  

! *** Local variables ***
  REAL :: shalf                         ! Half of support distance for Gaspari-Cohn


! **************************
! *** Set function value ***
! **************************

  IF (ctype == 1) THEN
     ! *********************************************
     ! *** Gaussian function scaled for f(0)=1.0 ***
     ! *********************************************

     ! Compute weight
     IF (length > 0.0) THEN

        value = exp(-distance*distance/ (2.0*length*length))

     ELSE

        IF (distance > 0.0) THEN
           value = 0.0
        ELSE
           value = 1.0
        END IF

     END IF

  ELSEIF (ctype == 2) THEN
     ! ************************************************************************
     ! *** 5th-order polynomial mimicking Gaussian but with compact support ***
     ! *** Equation (4.10) of Gaspari&Cohn, QJRMS125, 723 (1999)            ***
     ! ************************************************************************

     shalf = REAL(length) / 2.0

     ! Evaluate function
     cradnull: IF (length > 0.0) THEN

        cutoff: IF (distance <= length) THEN
           IF (distance <= length / 2.0) THEN
              value = -0.25 * (distance / shalf)**5 &
                   + 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   - 5.0 / 3.0 * (distance / shalf)**2 + 1.0
           ELSEIF (distance > length / 2.0 .AND. distance < length * 0.9) THEN
              value = 1.0 / 12.0 * (distance / shalf)**5 &
                   - 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   + 5.0 / 3.0 * (distance / shalf)**2 &
                   - 5.0 * (distance / shalf) &
                   + 4.0 - 2.0 / 3.0 * shalf / distance
           ELSEIF (distance >= length * 0.9 .AND. distance < length) THEN
              value = MAX(1.0 / 12.0 * (distance / shalf)**5 &
                   - 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   + 5.0 / 3.0 * (distance / shalf)**2 &
                   - 5.0 * (distance / shalf) &
                   + 4.0 - 2.0 / 3.0 * shalf / distance, 0.0)
           ELSE
              value = 0.0
           ENDIF
        ELSE cutoff
           value = 0.0
        END IF cutoff

     ELSE cradnull

        IF (distance > 0.0) THEN
           value = 0.0
        ELSE
           value = 1.0
        END IF

     END IF cradnull

  END IF


END SUBROUTINE PDAF_correlation_function


!-------------------------------------------------------------------------------
!> Set ensemble index to force an analysis step
!!
!! Helper routine for PDAF.
!! The routine overwrite member index of the ensemble 
!! state by local_dim_ens and the counter cnt_steps
!! by nsteps-1. This forces that the analysis
!! step is executed at the next call to PDAF_put_state
!! or PDAF_assimilate.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_force_analysis()

  USE PDAF_mod_core, &
       ONLY: member, local_dim_ens, nsteps, cnt_steps, step_obs, &
       step, screen, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world

  IMPLICIT NONE

! *** Set ensemble member ***

  member = local_dim_ens

  nsteps = cnt_steps + 1

  ! Only in case of using PDAF_assimilate, we need to reset the step counting
  IF (use_PDAF_assim) step_obs = step + nsteps - 1

  IF (screen>0 .AND. mype_world==0) THEN
     WRITE (*,'(a,5x,a,i8)') 'PDAF','!! Force analysis at step', step_obs
  END IF

END SUBROUTINE PDAF_force_analysis


!-------------------------------------------------------------------------------
!> Generate vector of random values
!!
!! Helper routine for PDAF.
!! The routine returns a vector of random number
!! of specified distribution.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code combining previous codes
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_generate_rndvec(len, vec, stddev, dist, iseed)

  implicit none

! *** Arguments ***
  integer, intent(in) :: len           ! Length of vector to process
  real, intent(inout) :: vec(:)        ! value to be perturbed
  real, intent(in)    :: stddev        ! Standard deviation of lognormal distribution
  integer, intent(in) :: dist          ! Distribution: 
                                       !   (1) Normal, N(vec, stddev^2)
                                       !   (2) Log-normal, logN(vec, stddev^2)
                                       !   (3) uniform: vec + stddev*U(0,1))
                                       !   (4) uniform: vec + stddev*U(-1,1))
                                       !   (5) Laplace
  integer, intent(in) :: iseed(4)      ! Seed for dlarnv (last entry has to be odd)

! *** Local variables ***
  integer :: i                         ! Counter
  real :: sigma2                       ! Variance for lognormal
  real, allocatable :: logval(:)       ! Logrithmic value of input value
  real, allocatable :: rndval(:)       ! Normal random value


  ALLOCATE(rndval(len))

  IF (dist==1 .OR. dist==2) THEN

     ! Generate Gaussian-distributed random number
     CALL dlarnv(3, iseed, len, rndval)

     IF (dist==1) THEN

        ! Normal distribution

        vec = vec + stddev*rndval

     ELSEIF (dist==2) THEN

        ! Log-normal distribution

        ALLOCATE(logval(len))

        sigma2 = log(1.0 + stddev*stddev)
        logval = log(vec) -0.5 * sigma2

        vec = exp(logval + sqrt(sigma2) * rndval) ! Eq. (A10) from Ciavatta et al. (2016)

        DEALLOCATE(logval)
     END IF

  ELSEIF (dist==3) THEN

     ! Generate uniformly distributed random vector (0,1)
     CALL dlarnv(1, iseed, len, rndval)

     vec = vec + stddev*rndval

  ELSEIF (dist==4) THEN

     ! Generate uniformly distributed random vector (-1,1)
     CALL dlarnv(2, iseed, len, rndval)

     vec = vec + stddev*rndval

  ELSEIF (dist==5) THEN

     ! Generate random Laplace distributed noise

     CALL dlarnv(1, iseed, len, rndval)

     DO i = 1, len
        rndval(i) = sign(1.0, (rndval(i)-0.5)) * log(1.0 - 2.0*ABS(rndval(i)-0.5))
     END DO

     ! Add noise to generate observations
     vec = vec + sqrt(stddev/2.0) * rndval

  END IF

  DEALLOCATE(rndval)

END SUBROUTINE PDAF_generate_rndvec


END MODULE PDAF_utils
