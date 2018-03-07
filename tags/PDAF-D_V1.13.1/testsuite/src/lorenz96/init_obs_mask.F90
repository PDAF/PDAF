!$Id: init_obs_mask.F90 1633 2016-08-15 08:02:27Z lnerger $
!BOP
!
! !ROUTINE: init_obs_mask - Initialize mask for gapps observations
!
! !INTERFACE:
SUBROUTINE init_obs_mask(dim)

! !DESCRIPTION:
! User-supplied routine for PDAF:
! Here a mask for observations is initialized. It
! has the size fo the state. The observation is used, 
! if the corresponding entry in the array obs_mask
! is one.
! The mask is read from the file file_obs_mask. This 
! file has to contain as many lines as the state dimension.
! A value of 1 indicates on observation at the corresponding
! grid point. A value of 0 is an observational gap.
!
! !REVISION HISTORY:
! 2010-05 - Lars Nerger
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: file_obs_mask, obs_mask, use_maskfile, numobs, dx_obs
  USE mod_parallel, &
       ONLY: abort_parallel
  USE parser, &
       ONLY: parse
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim            ! State dimension

! *** local variables ***
  INTEGER :: i                  ! Index of observation component
  INTEGER :: ioerr              ! Error flag for file handling
  CHARACTER(len=dim) :: maskstr ! String showing the observation availability
  CHARACTER(len=dim) :: maskleg ! Lagend for maskstr
  CHARACTER(len=1) :: maskstr1  ! String showing the observation availability


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(obs_mask(dim))

  WRITE (*,'(/1x, a)') 'Initialize mask for observations'

  maskf: IF (use_maskfile) THEN
     WRITE (*,'(5x, a, a)') 'Read mask from file: ', TRIM(file_obs_mask)

     OPEN(unit = 20, file = file_obs_mask, status = 'old', iostat=ioerr)

     IF (ioerr /= 0) THEN
        WRITE(*,*) 'ERROR in reading observation mask from file: ', &
             TRIM(file_obs_mask)
        CALL abort_parallel()
     END IF

     DO i = 1, dim
        READ(20,*) obs_mask(i)
     END DO

     CLOSE(20)

  ELSE maskf
     ! Initialize observation mask according to numobs

     IF (dx_obs==1) THEN
        WRITE (*,'(5x, a, i4)') 'Initialize mask for observations 1 to ', numobs
     ELSE IF (dx_obs>1) THEN
        WRITE (*,'(5x, a, i3, a)') 'Initialize mask for observations at each ', &
             dx_obs, '. grid point'
     END IF

     ! Reset numobs if larger than state
     if (numobs > dim) numobs = dim

     obs_mask(:) = 0
     IF (numobs>0) obs_mask(1) = 1

     DO i=1, dim
        IF (i*dx_obs+1 <= dim) THEN
           obs_mask(i*dx_obs+1) = 1
        END IF
     END DO

  END IF maskf

! *** Print observation mask to screen ***

  maskstr = ''
  maskleg = ''
  DO i = 1, dim
     WRITE (maskstr1,'(i1)') obs_mask(i)
     maskstr = TRIM(maskstr) // maskstr1
     IF (MOD(i,10)==0) THEN
        maskleg = TRIM(maskleg) // '|'
     ELSE
        maskleg = TRIM(maskleg) // '_'
     END IF
  END DO
  WRITE (*,'(5x, a)') 'Mask (mark | at each 10th entry)'
  WRITE (*,'(5x, a)') maskleg
  WRITE (*,'(5x, a)') maskstr

  
END SUBROUTINE init_obs_mask
