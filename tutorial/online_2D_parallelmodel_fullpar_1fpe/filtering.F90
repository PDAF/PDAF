!> Program part to run exclusively on filter processes
!!
!! Filtering routine running the data assimilation
!! part for the case that the data assimilation and
!! model integrations are run on distinct sets of 
!! processes.
!!
!! __Revision history:__
!! * 2014-04 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE filtering()

  USE mod_model, &               ! Model variables
       ONLY: nx, ny, nx_p, total_steps
  USE mod_parallel_model, &      ! Model parallelization variables
       ONLY: abort_parallel
  USE mod_parallel_pdaf, &       ! PDAF parallelization variables
       ONLY: mype_filter, npes_filter

  IMPLICIT NONE

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
