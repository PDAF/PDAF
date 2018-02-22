!$Id$
!BOP
!
! !ROUTINE: get_modelerror --- initialize random model error
!
! !INTERFACE:
SUBROUTINE get_modelerror(dim_state, modelerror, rnd_type)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Get a single realization of modelerror.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner
! Later revisions - see svn log

  USE mod_parallel_pdaf, & 
       ONLY: mype_world
  USE mod_assimilation, &
       ONLY: modelerr_amp

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT(in) :: dim_state
  REAL,INTENT(out)    :: modelerror(dim_state)
  INTEGER, INTENT(in) :: rnd_type   ! 0 - normal random numbers
                                    ! 1 - uniform random numbers

! Local variables
  INTEGER  :: i, j, k
  INTEGER, SAVE :: seed(4)
  INTEGER, SAVE :: firsttime = 1


! Initialize seed
  IF (firsttime==1) THEN
     seed(1) = 1
     seed(2) = 2
     seed(3) = 7+mype_world
     seed(4) = 17

     firsttime = 0
  ENDIF

! Initialize error
  IF (rnd_type==0) THEN
     CALL dlarnv(3,seed,dim_state,modelerror)
  ELSEIF (rnd_type==1) THEN
     CALL dlarnv(2,seed,dim_state,modelerror)
  END IF

  modelerr_amp = 0.01

  modelerror = modelerror * modelerr_amp


END SUBROUTINE get_modelerror
 

