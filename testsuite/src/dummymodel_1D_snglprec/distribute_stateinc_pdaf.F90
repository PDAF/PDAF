!$Id: distribute_stateinc_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: distribute_stateinc_pdaf --- Add analysis increment to model fields
!
! !INTERFACE:
SUBROUTINE distribute_stateinc_pdaf(dim_p, state_inc_p, first, steps)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase of the filter from PDAF\_incremental
! supplying the analysis state increment.
! The routine has to compute the fraction of 
! the increment to be added to the model state 
! at each time step. Further, it has to transform 
! the increment vector into increments of the 
! fields of the model (typically available 
! trough a module).
!
! The routine is executed by each process that 
! is participating in the model integrations.
!
! For the dummy model and PDAF with domain 
! decomposition the state vector and the model 
! field are identical. The increment vector 
! divided by the number of time steps in the 
! forecast phase is added to the model field.
!
! !REVISION HISTORY:
! 2006-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_model
  USE mod_model, &
       ONLY: field

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! Dimension of PE-local state
  REAL, INTENT(in) :: state_inc_p(dim_p) ! PE-local state vector
  INTEGER, INTENT(in) :: first           ! Flag for first call of each forecast
  INTEGER, INTENT(in) :: steps           ! number of time steps in forecast

! !CALLING SEQUENCE:
! Called by: PDAF_incremental   (as U_dist_stateinc)
!EOP  

! local variables
  INTEGER, SAVE :: allocfirst = 1        ! Flag for allocation
  REAL, SAVE, ALLOCATABLE :: inc_save(:) ! Stored state increment


! *******************************
! *** prepare state increment ***
!********************************

  IF (first > 0) THEN
     ! At begin of each forecast phase compute increment
     ! per update step. (E.g., at each time step)

     IF (allocfirst == 1) THEN
        ALLOCATE(inc_save(dim_p))
        allocfirst = 0
     ENDIF

     IF (mype_model == 0) THEN
        WRITE (*, '(8x, a)') '--- Initialize Incremental updating'
     END IF

     ! compute increment per time step
     inc_save = state_inc_p / REAL(steps)
  ENDIF


! *************************************
! *** Add increment to model fields ***
!**************************************

  field(:) = field(:) + inc_save(:)


END SUBROUTINE distribute_stateinc_pdaf
