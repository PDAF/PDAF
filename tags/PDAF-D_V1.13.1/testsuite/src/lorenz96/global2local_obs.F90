!$Id: global2local_obs.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: global2local_obs --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
SUBROUTINE global2local_obs(domain, step, dim_obs, dim_obs_l, mstate, &
     mstate_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the analysis step
! on each of the local analysis domains.
! It has to restrict the full vector of all 
! observations required for the loop of localized 
! analyses on the PE-local domain to the current 
! local analysis domain.
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
       ONLY: local_range, local_range2, use_obs_mask, obsindx_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain            ! Current local analysis domain
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: dim_obs           ! Dimension of full PE-local obs. vector
  INTEGER, INTENT(in) :: dim_obs_l         ! Local dimension of observation vector
  REAL, INTENT(in)    :: mstate(dim_obs)   ! Full PE-local obs. vector
  REAL, INTENT(out)   :: mstate_l(dim_obs_l)   ! Obs. vector on local domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_g2l_obs)
!EOP


! *** local variables ***
  INTEGER :: i, ilow, iup  ! Counters


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  obsgaps: IF (.NOT. use_obs_mask) THEN
     ! This is the case that the full state is observed
     ilow = domain - local_range
     iup = domain + local_range2

     ! Perform localization
     IF (ilow >= 1 .AND. iup <= dim_state) THEN
        ! Observed region completely within observed region
        DO i = ilow, iup
           mstate_l(i - ilow + 1) = mstate(i)
        END DO
     ELSE IF (ilow < 1) THEN
        ! Use lower periodic BC
        DO i = ilow + dim_state, dim_state
           mstate_l(i-ilow-dim_state+1) = mstate(i)
        END DO
        DO i = 1, iup
           mstate_l(i-ilow+1) = mstate(i)
        END DO
     ELSE IF (iup > dim_state) THEN
        ! Use upper periodic BC
        DO i = ilow, dim_state
           mstate_l(i-ilow+1) = mstate(i)
        END DO
        DO i = 1, iup - dim_state
           mstate_l(i+dim_state-ilow+1) = mstate(i)
        END DO
     END IF

  ELSE obsgaps
     ! Gappy observations

     DO i = 1, dim_obs_l
        mstate_l(i) = mstate(obsindx_l(i))
     END DO

  END IF obsgaps

END SUBROUTINE global2local_obs
