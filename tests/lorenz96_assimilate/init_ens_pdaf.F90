!$Id$
!BOP
!
! !ROUTINE: init_ens_pdaf --- Initialize ensemble for SEIK
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim, dim_ens, state, Uinv, &
     ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens model 
! states. Here, we supply three methods to 
! initialize the ensemble:
! eof - By second-order exact sampling. 
!       This follows Pham (2001) and was used
!       in most of our papers.
! rnd - By random sampling form a long state
!       trajectory. This method is often described
!       In papers on the EnKF and the ensemble
!       square-root filters.
! tru - This variant is only used when synthetic
!       observation are generated. In this case
!       all ensemble states are equal.
!
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &
       ONLY: timeit
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: type_ensinit

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype          ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
!  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(inout) :: Uinv(1,1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_init_ens)
! Calls: init_ens_eof
! Calls: init_ens_rnd
! Calls: collect_state
!EOP

! *** local variables ***
   INTEGER :: i       ! counter


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** Generate full ensemble ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'

  CALL timeit(6, 'new')

  IF (TRIM(type_ensinit) == 'eof') THEN
     ! Initialize by 2nd-order exact sampling from EOFs
     CALL init_ens_eof(dim, dim_ens, state, ens, flag)
  ELSE IF (TRIM(type_ensinit) == 'rnd') THEN
     ! Initialize by random sampling from state trajectory
     CALL init_ens_rnd(dim, dim_ens, state, ens, flag)
  ELSE IF (TRIM(type_ensinit) == 'tru') THEN
     ! Initialize from true initial condition
     WRITE (*, '(9x, a)') '--- generate from model initial state'

     DO i=1, dim_ens
        CALL collect_state_pdaf(dim, ens(:,i))
     END DO
     
  END IF

  CALL timeit(6, 'old')


END SUBROUTINE init_ens_pdaf
