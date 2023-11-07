!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!prodrinva_l_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                      'prodrinva_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: prodrinva_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_l_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the product of the inverse of the local
! observation error covariance matrix with
! the matrix of locally observed ensemble 
! perturbations.
! Next to computing the product,  a localizing 
! weighting (similar to covariance localization 
! often used in EnKF) can be applied to matrix A.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange, obs_index_p, &
        rms_obs, distance 
  USE mod_parallel_pdaf, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  ! Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis    (as U_prodRinvA_l)
! Called by: PDAF_lestkf_analysis   (as U_prodRinvA_l)
! Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)
!EOP

! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
  INTEGER :: domain        ! Global domain index
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)     ! Localization weights
 ! REAL, ALLOCATABLE :: distance(:)   ! Localization distance
  REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL    :: meanvar                 ! Mean variance in observation domain
  REAL    :: svarpovar               ! Mean state plus observation variance
  REAL    :: var_obs                 ! Variance of observation error

! *** NO CHANGES REQUIRED BELOW IF OBSERVATION ERRORS ARE CONSTANT ***

! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p

  ! Screen output
  IF (verbose == 1) THEN
     WRITE (*, '(8x, a, f12.3)') &
           '--- Use global rms for observations of ', rms_obs
     WRITE (*, '(8x, a, 1x)') &
          '--- Domain localization'
     WRITE (*, '(12x, a, 1x, f12.2)') &
          '--- Local influence radius', local_range

     IF (locweight > 0) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observation errors'

        IF (locweight == 3) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 4) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF
  
  ! *** initialize numbers (this is for constant observation errors)
  ! Set observation variance and inverse here
  ivariance_obs = 1.0 / rms_obs**2
  var_obs = rms_obs**2

! ********************************
! *** Initialize weight array. ***
! ********************************

  ! Allocate weight array
  ALLOCATE(weight(dim_obs_l))

  if (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  else if (locweight == 1) THEN
     ! Exponential weighting
     wtype = 1
     rtype = 0
  ELSE IF (locweight == 2 .OR. locweight == 3 .OR. locweight == 4) THEN
     ! 5th-order polynomial (Gaspari&Cohn, 1999)
     wtype = 2

     IF (locweight < 3) THEN
        ! No regulated weight
        rtype = 0
     ELSE
        ! Use regulated weight
        rtype = 1
     END IF

  end if

  IF (locweight == 4) THEN
     ! Allocate array for single observation point
     ALLOCATE(A_obs(1, rank))
  END IF

  DO i=1, dim_obs_l
     ! Control verbosity of PDAF_local_weight
     IF (verbose==1 .AND. i==1) THEN
        verbose_w = 1
     ELSE
        verbose_w = 0
     END IF

     IF (locweight /= 4) THEN
        ! All localizations except regulated weight based on variance at 
        ! single observation point
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             dim_obs_l, rank, A_l, var_obs, weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        A_obs(1,:) = A_l(i,:)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance(i), &
             1, rank, A_obs, var_obs, weight(i), verbose_w)
     END IF
  END DO

  IF (locweight == 4) DEALLOCATE(A_obs)


! ********************
! *** Apply weight ***
! ********************

  DO j = 1, rank
     DO i = 1, dim_obs_l
         C_l(i, j) =  ivariance_obs * weight(i) * A_l(i, j)
     END DO
  END DO

! *** Clean up ***
  DEALLOCATE(weight)
  
END SUBROUTINE prodRinvA_l_pdaf
