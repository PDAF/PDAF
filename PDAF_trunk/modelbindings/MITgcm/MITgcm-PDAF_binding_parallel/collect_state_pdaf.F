C     $Id$
#include "CPP_OPTIONS.h"
C!BOP
C
C     !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
C
C     !INTERFACE:
      SUBROUTINE collect_state_pdaf(dim_p, state_p)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
C
C     This subroutine is called during the forecast 
C     phase from PDAF\_put\_state\_X after the 
C     propagation of each ensemble member. 
C     The supplied state vector has to be initialized
C     from the model fields (typically via a module). 
C     With parallelization, MPI communication might be 
C     required to initialize state vectors for all 
C     subdomains on the model PEs. 
C
C     The routine is executed by each process that is
C     participating in the model integrations.
C
C     !REVISION HISTORY:
C     2004-11 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: myIter, myThid

      IMPLICIT NONE
  
C     !ARGUMENTS:
      INTEGER, INTENT(in) :: dim_p ! PE-local state dimension
      _RL, INTENT(inout) :: state_p(dim_p) ! local state vector

C     !Local variables:
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

      INTEGER :: i, j, k, bi, bj, koffset
      CHARACTER*(MAX_LEN_MBUF) ::  msgBuf
  
C     !CALLING SEQUENCE:
C     Called by: PDAF_put_state_X   (as U_coll_state)
C!EOP
  
C     *************************************************
C     *** Initialize state vector from model fields ***
C     *************************************************

      koffset = 0

C     *** Sea surface height (etaN) ***
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j=1,sNy
               DO i=1,sNx
                  koffset = koffset + 1
                  state_p(koffset) = etaN(i,j,bi,bj)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      write (*,*) 'KOFFSET ETA: ', koffset
C     *** TEMP (theta) ***
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     koffset = koffset + 1
                     state_p(koffset) = theta(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      write (*,*) 'KOFFSET THETA: ', koffset
C     *** SALT (salt) ***
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     koffset = koffset + 1
                     state_p(koffset) = salt(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      write (*,*) 'KOFFSET SALT: ', koffset
C     *** U velocity (uVel) ***
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     koffset = koffset + 1
                     state_p(koffset) = uVel(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      write (*,*) 'KOFFSET UVEL: ', koffset      
C     *** V velocity (vVel) ***
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     koffset = koffset + 1
                     state_p(koffset) = vVel(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      write (*,*) 'KOFFSET VVEL: ', koffset

      IF ( koffset < dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'COLLECT_STATE_PDAF: koffset < dim_p:',koffset,' < ',dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF
      IF ( koffset > dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'COLLECT_STATE_PDAF: koffset > dim_p:',koffset,' > ',dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF

      END
