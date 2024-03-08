C     $Id$
C
C     !ROUTINE: init_ens_pdaf --- Initialize ensemble
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The routine is called when the filter is
C     initialized in PDAF\_filter\_init.  It has
C     to initialize an ensemble of dim\_ens states.
C     Typically, the ensemble will be directly read from files.
C
C     The routine is called by all filter processes and
C     initializes the ensemble for the PE-local domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p,
     &     Uinv, ens_p, flag)

C     !USES:
      USE mod_parallel_pdaf,
     &     ONLY: mype_filter

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: filtertype                  ! Type of filter to initialize
      INTEGER, INTENT(in) :: dim_p                       ! PE-local state dimension
      INTEGER, INTENT(in) :: dim_ens                     ! Size of ensemble
      REAL*8, INTENT(inout) :: state_p(dim_p)            ! PE-local model state
C     It is not necessary to initialize the array 'state_p' for SEIK. 
C     It is available here only for convenience and can be used freely.
      REAL*8, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
      REAL*8, INTENT(out)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble
      INTEGER, INTENT(inout) :: flag                     ! PDAF status flag

C     !CALLING SEQUENCE:
C     Called by: PDAF_filter_init    (as U_init_ens)


#include "CPP_OPTIONS.h"
#include "SIZE.h"
#include "DYNVARS.h"

C     *** local variables ***
      INTEGER :: member  ! Counters


C     **********************
C     *** INITIALIZATION ***
C     **********************

C     *** Generate full ensemble on filter-PE 0 ***
      IF (mype_filter==0) THEN
         WRITE (*, '(/a, 9x, a)') 'PDAFuser: ',
     &        'Initialize state ensemble'
         WRITE (*, '(a, 9x, a)') 'PDAFuser: ',
     &        '--- read ensemble from files'
         WRITE (*, '(a, 9x, a, i5)') 'PDAFuser: ',
     &        '--- Ensemble size:  ', dim_ens
      END IF


C     ********************************
C     *** Read ensemble from files ***
C     ********************************

      ! This is a very simplified ensemble initialization
      ! just using a scale according to the ensemble member
      ! index. This will not work with real cases.
      CALL collect_state_pdaf(dim_p, state_p)

      DO member = 1, dim_ens
         ens_p(:, member)  = state_p(:) + float(member-1)*1.d-2
      ENDDO


C     ****************
C     *** clean up ***
C     ****************

      END SUBROUTINE init_ens_pdaf
