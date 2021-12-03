C     $Id$
C     !BOP
C
C     !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in all filters
C
C     This subroutine is called during the forecast
C     phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
C     after the propagation of each ensemble member.
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
C     !INTERFACE:
      SUBROUTINE collect_state_pdaf(dim_p, state_p)

C     !USES:
      USE mod_assimilation,
     &     ONLY: myThid, id, off_fields

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: dim_p ! PE-local state dimension
      REAL(kind=8), INTENT(inout) :: state_p(dim_p) ! local state vector

C     !Local variables:
#include "CPP_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

      INTEGER :: i, j, k, bi, bj, cnt
      CHARACTER*(MAX_LEN_MBUF) ::  msgBuf

C     !CALLING SEQUENCE:
C     Called by: PDAF_put_state_X   (as U_coll_state)
C!EOP

C     *************************************************
C     *** Initialize state vector from model fields ***
C     *************************************************

C     *** Sea surface height (etaN) ***
      cnt = off_fields(id%etan)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j=1,sNy
               DO i=1,sNx
                  cnt = cnt + 1
                  state_p(cnt) = etaN(i,j,bi,bj)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     *** TEMP (theta) ***
      cnt = off_fields(id%theta)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     cnt = cnt + 1
                     state_p(cnt) = theta(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     *** SALT (salt) ***
      cnt = off_fields(id%salt)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     cnt = cnt + 1
                     state_p(cnt) = salt(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     *** U velocity (uVel) ***
      cnt = off_fields(id%uvel)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     cnt = cnt + 1
                     state_p(cnt) = uVel(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     *** V velocity (vVel) ***
      cnt = off_fields(id%vvel)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO k=1,Nr
               DO j=1,sNy
                  DO i=1,sNx
                     cnt = cnt + 1
                     state_p(cnt) = vVel(i,j,k,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF ( cnt < dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'COLLECT_STATE_PDAF: cnt < dim_p:',cnt,' < ',dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF
      IF ( cnt > dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'COLLECT_STATE_PDAF: cnt > dim_p:',cnt,' > ',dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF

      END
