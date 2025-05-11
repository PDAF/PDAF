!$Id: initialize.F90 1426 2013-09-25 15:21:35Z lnerger $
!BOP
!
! !ROUTINE: filtering --- Program part to run exclusively on filter processes
!
! !INTERFACE:
SUBROUTINE filtering()

! !DESCRIPTION:
! Initialization routine for the simple 2D model with
! parallelization of the model. 
!
! In this routine only the dimensions are initialized
! which are required on the filter processes. No state
! information is read. 
!
! !REVISION HISTORY:
! 2014-04 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: nx, ny, nx_p, total_steps
  USE mod_parallel_model, &
       ONLY: abort_parallel
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
!EOP

! *** local variables ***
  INTEGER :: step


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
  nx = 36          ! Extent of grid in x-direction
  ny = 18          ! Extent of grid in y-direction
  total_steps = 18 ! Number of time steps to perform

! *** Screen output ***
  IF (mype_filter == 0) THEN
     WRITE (*, '(1x, a)') 'PDAF-side: INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
     WRITE (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
     WRITE (*, '(10x,a,i4)') 'Time steps', total_steps
  END IF

! *** Initialize size of local nx for parallelization ***
  IF (npes_filter==1 .OR. npes_filter==2 .OR. npes_filter==3 .OR. npes_filter==4 .OR. &
       npes_filter==6 .OR.npes_filter==9) THEN
     ! Split x-direction in chunks of equal size
     nx_p = nx / npes_filter
  ELSE
     WRITE (*,*) 'ERROR: Invalid number of processes'
     CALL abort_parallel()
  END IF

  IF (mype_filter == 0) THEN
     WRITE (*, '(/2x, a, i3, a)') &
          '-- Domain decomposition over', npes_filter, ' PEs'
     WRITE (*, '(2x,a,i3,a,i3)') &
          '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
  END IF

  ! Initialize PDAF
  CALL init_pdaf()

  IF (mype_filter==0) WRITE (*, '(1x, a)') 'PDAF-side: START INTEGRATION'

  stepping: DO step = 1 , total_steps

     CALL assimilate_pdaf()

  END DO stepping


END SUBROUTINE filtering
