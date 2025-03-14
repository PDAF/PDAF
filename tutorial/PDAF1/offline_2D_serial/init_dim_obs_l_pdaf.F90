!$Id: init_dim_obs_l_pdaf.F90 1860 2017-12-18 13:04:35Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analysis domain.
!
! Implementation for the 2D offline example.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: nx, ny, &
       cradius, coords_obs_f, id_lobs_in_fobs, coords_l, distance_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
!EOP


! *** local variables ***
  INTEGER :: i, cnt                   ! Counters
  REAL :: limits_x(2), limits_y(2)    ! Coordinate limits for observation domain
  REAL :: distance                    ! Distance between observation and analysis domain


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  !Determine coordinate limits for observation domain
  limits_x(1) = coords_l(1) - cradius
  IF (limits_x(1) < 1.0) limits_x(1) = 1.0
  limits_x(2) = coords_l(1) + cradius
  IF (limits_x(2) > REAL(nx)) limits_x(2) = REAL(nx)

  limits_y(1) = coords_l(2) - cradius
  IF (limits_y(1) < 1.0) limits_y(1) = 1.0
  limits_y(2) = coords_l(2) + cradius
  IF (limits_y(2) > REAL(ny)) limits_y(2) = REAL(ny)

  ! Count observations within cradius
  dim_obs_l = 0
  DO i = 1, dim_obs_f
     IF (coords_obs_f(1, i) >= limits_x(1) .AND. coords_obs_f(1, i) <= limits_x(2) .AND. &
          coords_obs_f(2, i) >= limits_y(1) .AND. coords_obs_f(2, i) <= limits_y(2)) THEN
        
        distance = SQRT((coords_l(1) - coords_obs_f(1,i))**2 + &
             (coords_l(2) - coords_obs_f(2,i))**2)
        IF (distance <= cradius) dim_obs_l = dim_obs_l + 1

     END IF
  END DO

  ! Initialize index array for local observations in full observed vector
  ! and array of distances for local observations
  IF (ALLOCATED(id_lobs_in_fobs)) DEALLOCATE(id_lobs_in_fobs)
  ALLOCATE(id_lobs_in_fobs(dim_obs_l))
  IF (ALLOCATED(distance_l)) DEALLOCATE(distance_l)
  ALLOCATE(distance_l(dim_obs_l))

  cnt = 0
  DO i = 1, dim_obs_f
     IF (coords_obs_f(1, i) >= limits_x(1) .AND. coords_obs_f(1, i) <= limits_x(2) .AND. &
          coords_obs_f(2, i) >= limits_y(1) .AND. coords_obs_f(2, i) <= limits_y(2)) THEN

        distance = SQRT((coords_l(1) - coords_obs_f(1,i))**2 + &
             (coords_l(2) - coords_obs_f(2,i))**2)
        IF (distance <= cradius) THEN
           cnt = cnt + 1
           id_lobs_in_fobs(cnt) = i
           distance_l(cnt) = distance
        END IF
     END IF
  END DO

END SUBROUTINE init_dim_obs_l_pdaf
