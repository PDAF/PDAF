!$Id: init_dim_obs_local.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_local --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_local(domain, step, dim_obs, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analsis domain.
!
! This variant is for the Lorenz96 model without
! parallelization. A local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state
  USE mod_assimilation, &
       ONLY: local_range, local_range2, use_obs_mask, obsindx, obsindx_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain     ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs    ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
!EOP

! *** local variables ***
  INTEGER :: i, k, ilow, iup            ! Counters


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  obsgaps: IF (.NOT. use_obs_mask) THEN
     ! Variant if full state is observed

     dim_obs_l = 1 + local_range + local_range2

     ! local dimension is larger than state dimension:
     ! reset to state_dimension
     IF (dim_obs_l > dim_state) dim_obs_l = dim_state

  ELSE
     ! Variant for gappy observations
     ! Here we initialize dim_obs_l and also the index array obsindx_l that
     ! points to the required elements in the global observation vector.

     dim_obs_l = 0
     obsindx_l = -1

     ilow = domain - local_range
     iup = domain + local_range2

     ! Perform localization
     IF (ilow >= 1 .AND. iup <= dim_state) THEN
        ! Observed region completely within observed region
        k = 1
        DO i = 1, dim_obs
           IF (obsindx(i) >= ilow .AND. obsindx(i) <= iup) THEN
              dim_obs_l = dim_obs_l + 1
              obsindx_l(k) = i
              k = k + 1
           END IF
        END DO
     ELSE IF (ilow < 1) THEN
        ! Use lower periodic BC
        k = 1
        DO i = 1, dim_obs
           IF (obsindx(i) >= ilow + dim_state .AND. obsindx(i) <= dim_state) THEN
              dim_obs_l = dim_obs_l + 1
              obsindx_l(k) = i
              k = k + 1
           END IF
        END DO
        DO i = 1, dim_obs
           IF (obsindx(i) >= 1 .AND. obsindx(i) <= iup) THEN
              dim_obs_l = dim_obs_l + 1
              obsindx_l(k) = i
              k = k + 1
           END IF
        END DO
     ELSE IF (iup > dim_state) THEN
        ! Use upper periodic BC
        k = 1
        DO i = 1, dim_obs
           IF (obsindx(i) >= ilow .AND. obsindx(i) <= dim_state) THEN
              dim_obs_l = dim_obs_l + 1
              obsindx_l(k) = i
              k = k + 1
           END IF
        END DO
        DO i = 1, dim_obs
           IF (obsindx(i) >= 1 .AND. obsindx(i) <= iup - dim_state) THEN
              dim_obs_l = dim_obs_l + 1
              obsindx_l(k) = i
              k = k + 1
           END IF
        END DO
     END IF
  END IF obsgaps

END SUBROUTINE init_dim_obs_local

