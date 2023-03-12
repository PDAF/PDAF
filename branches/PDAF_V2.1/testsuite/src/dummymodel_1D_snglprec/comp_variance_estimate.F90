!$Id: comp_variance_estimate.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: comp_variance_estimate --- Compute estimated variance for subtype=3
!
! !INTERFACE:
SUBROUTINE comp_variance_estimate(step, dim_p, dim_ens, dim_u, state_p, Uinv, &
     ens_p, variance_p, initialstep)

! !DESCRIPTION:
! Helper routine for for pre/poststep routine.
! Used in the filters: SEIK/ETKF/ESTKF
! It is also suited for the local filter variants. However, as
! each local analysis domain has in general a different
! transform matrix, and only the transform matrix of the last
! domain is available after the analysis, this computation
! will most likely be incorrect.
! 
! The routine is called by prepoststep_ens_pdaf if the estimated
! error covariance should be computed in the case that a filter
! is used with a fixed ensemble-represented error covariance 
! matrix (subtype=3). In this case the ensemble spread is not
! changed during the analysis step. To compute the analysis
! error estimate, one has to mimic the update of the ensemble 
! spread to compute the analysis error. This routine performs
! this update to compute the analysis variance estimate.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2012-03 - Lars Nerger - Initial code extractd from prepoststep routine
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_model, &
       ONLY: step_null
  USE mod_assimilation, &
       ONLY: covartype, filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_u       ! Size of Uinv
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! For subtype=3, the array 'state_p' is initialized after the first
  ! forecast phase. It should not be changed any more.
  REAL, INTENT(inout) :: Uinv(dim_u, dim_u)    ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL, INTENT(inout) :: variance_p(dim_p)     ! PE-local forecast/analysis state
  LOGICAL, INTENT(in) :: initialstep  ! Whether routine is called at the initial time step

! !CALLING SEQUENCE:
! Called by: prepoststep_ens_pdaf  
! Called by: prepoststep_etkf_pdaf  

! *** local variables ***
  INTEGER :: i, j, row, member         ! counters
  INTEGER :: dgesv_info                ! output flag for DGESV
  REAL :: rdim_ens                     ! dim_ens in real format
  REAL :: fac                          ! Temporary factor
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL, ALLOCATABLE :: Ttrans(:,:)     ! matrix T^T
  REAL, ALLOCATABLE :: TUT(:,:)        ! temporary matrix TUT^T
  REAL, ALLOCATABLE :: tempUinv(:,:)   ! temporary matrix Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for SGESV
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting


  ! Initialize numbers
  invdim_ens    = 1.0 / REAL(dim_ens)  

  ! Compute ensemble mean state
  ! We compute it only at the initial step, because later state_p 
  ! holds the forecast or analysis step
  IF (initialstep) THEN
     ! local 
     state_p = 0.0
     DO member = 1, dim_ens
        DO i = 1, dim_p
           state_p(i) = state_p(i) + ens_p(i, member)
        END DO
     END DO
     state_p(:) = invdim_ens * state_p(:)
  END IF

  haveSEIKorESTKF: IF (filtertype==1 .OR. filtertype==3) THEN
! *** Variant for SEIK

     ! Allocate fields
     ALLOCATE(tempUinv(dim_ens - 1, dim_ens - 1))
     ALLOCATE(Ttrans(dim_ens - 1, dim_ens))
     ALLOCATE(TUT(dim_ens, dim_ens))
     ALLOCATE(ipiv(dim_ens - 1)) 
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(3, 'r', (dim_ens-1)**2 + (dim_ens-1)*dim_ens + dim_ens**2)
        CALL memcount(3, 'i', dim_ens-1)
     END IF

     ! Initialize matrix T^T
     DO i = 1, dim_ens - 1
        DO j = 1, dim_ens
           Ttrans(i, j) = -invdim_ens
        END DO
     END DO
     DO i = 1, dim_ens - 1
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO

     IF (step > 0 .AND. (step - step_null /= 0)) THEN
        ! Initialize temporary Uinv (We must not change Uinv here!)
        tempUinv(1:dim_ens-1, 1:dim_ens-1) = Uinv(1:dim_ens-1, 1:dim_ens-1)
     ELSE IF (step < 0 .OR. (step - step_null == 0)) THEN
        ! Initialize invariant Uinv (dim_ens T T^T)
        IF (covartype == 1) THEN
           ! For covariance matrix with factor r^-1 (new SEIK - real ensemble)
           rdim_ens = REAL(dim_ens - 1)
        ELSE
           ! For covariance matrix with factor (r+1)^-1 (old SEIK)
           rdim_ens = REAL(dim_ens)
        END IF
        CALL sgemm('n', 't', dim_ens - 1, dim_ens - 1, dim_ens, &
             rdim_ens, Ttrans, dim_ens - 1, Ttrans, dim_ens - 1, &
             0.0, tempUinv, dim_ens - 1)
     END IF

     ! call solver - compute W = U T^T
     CALL sgesv(dim_ens - 1, dim_ens, tempUinv, dim_ens - 1, ipiv, &
          Ttrans, dim_ens - 1, dgesv_info) 

     ! Compute T W = T U T^T using operation in subroutine
     CALL PDAF_seik_TtimesA(dim_ens - 1, dim_ens, Ttrans, TUT)

  ELSE IF (filtertype==4 .OR. filtertype==5) THEN
! *** Variant for ETKF and LETKF

     ! Allocate fields
     ALLOCATE(tempUinv(dim_ens, dim_ens))
     ALLOCATE(Ttrans(dim_ens, dim_ens))
     ALLOCATE(TUT(dim_ens, dim_ens))
     ALLOCATE(ipiv(dim_ens)) 
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(3, 'r', 3*dim_ens**2)
        CALL memcount(3, 'i', dim_ens)
     END IF

     ! Initialize matrix T^T
     DO i = 1, dim_ens
        DO j = 1, dim_ens
           Ttrans(i, j) = -invdim_ens
        END DO
     END DO
     DO i = 1, dim_ens
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO

     IF (step > 0 .AND. (step - step_null /= 0)) THEN
        ! Initialize temporary Uinv (We must not change Uinv here!)
        tempUinv(1:dim_ens, 1:dim_ens) = Uinv(1:dim_ens, 1:dim_ens)

        ! call solver - compute W = U T^T
        CALL sgesv(dim_ens, dim_ens, tempUinv, dim_ens, ipiv, &
             Ttrans, dim_ens, dgesv_info) 

        ! Compute T W = T U T^T using operation in subroutine     
        CALL PDAF_etkf_Tleft(dim_ens, dim_ens, Ttrans)
        TUT = Ttrans

     ELSE IF (step < 0 .OR. (step - step_null == 0)) THEN
        ! After the forecast or at the initial time
        TUT = Ttrans/REAL(dim_ens-1)
     END IF
  ELSE IF (filtertype==6 .OR. filtertype==7) THEN
! *** Variant for ESTKF
     
     ! Allocate fields
     ALLOCATE(tempUinv(dim_ens - 1, dim_ens - 1))
     ALLOCATE(Ttrans(dim_ens - 1, dim_ens))
     ALLOCATE(TUT(dim_ens, dim_ens))
     ALLOCATE(ipiv(dim_ens - 1)) 
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(3, 'r', (dim_ens-1)**2 + (dim_ens-1)*dim_ens + dim_ens**2)
        CALL memcount(3, 'i', dim_ens-1)
     END IF

     ! Initialize matrix T^T
     rdim_ens = REAL(dim_ens)
     fac = 1.0 / (rdim_ens * (1.0/SQRT(rdim_ens)+1.0))
     DO i = 1, dim_ens - 1
        DO j = 1, dim_ens - 1
           Ttrans(i, j) = - fac
        END DO
     END DO
     DO i = 1, dim_ens - 1
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO
     DO i = 1, dim_ens - 1
        Ttrans(i, dim_ens) = - 1.0 / SQRT(rdim_ens)
     END DO

     IF (step > 0 .AND. (step - step_null /= 0)) THEN
        ! Initialize temporary Uinv (We must not change Uinv here!)
        tempUinv(1:dim_ens-1, 1:dim_ens-1) = Uinv(1:dim_ens-1, 1:dim_ens-1)
     ELSE IF (step < 0 .OR. (step - step_null == 0)) THEN
        ! Initialize invariant Uinv (dim_ens T T^T)
        rdim_ens = REAL(dim_ens - 1)

        CALL sgemm('n', 't', dim_ens - 1, dim_ens - 1, dim_ens, &
             rdim_ens, Ttrans, dim_ens - 1, Ttrans, dim_ens - 1, &
             0.0, tempUinv, dim_ens - 1)
     END IF

     ! call solver - compute W = U T^T
     CALL sgesv(dim_ens - 1, dim_ens, tempUinv, dim_ens - 1, ipiv, &
          Ttrans, dim_ens - 1, dgesv_info) 

     ! Compute Omega W = T U T^T using operation in subroutine
     CALL PDAF_estkf_OmegaA(dim_ens - 1, dim_ens, Ttrans, TUT)

  END IF haveSEIKorESTKF

  ! Compute local sampled variances
  variance_p(:) = 0.0
  DO i = 1, dim_ens
     DO j = 1, dim_ens
        DO row = 1, dim_p
           variance_p(row) = variance_p(row) &
                + ens_p(row, j) * ens_p(row, i) * TUT(i, j)
        END DO
     END DO
  END DO
  
  DEALLOCATE(tempUinv, ipiv, Ttrans, TUT)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE comp_variance_estimate
