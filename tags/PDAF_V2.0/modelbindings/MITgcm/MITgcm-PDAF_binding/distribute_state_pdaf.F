C     $Id$

C     !BOP
C
C     !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in all filters.
C
C     During the forecast phase of the filter this
C     subroutine is called from PDAF\_get\_state
C     supplying a model state which has to be evolved.
C     The routine has to initialize the fields of the
C     model (typically available through a module) from
C     the state vector of PDAF. With parallelization,
C     MPI communication might be required to
C     initialize all subdomains on the model PEs.
C
C     The routine is executed by each process that is
C     participating in the model integrations.
C
C     !REVISION HISTORY:
C     2004-10 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
#include "CPP_OPTIONS.h"
      
C     !INTERFACE:
      SUBROUTINE distribute_state_pdaf(dim_p, state_p)

C     !USES:
      USE mod_assimilation,
     &     ONLY: myThid, id, off_fields

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: dim_p ! PE-local state dimension
      REAL(kind=8), INTENT(inout) :: state_p(dim_p) ! PE-local state vector

C     !Local variables:
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     Local variables
      INTEGER i,j,k,bi,bj,cnt
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     !CALLING SEQUENCE:
C     Called by: PDAF_get_state   (as U_dist_state)
C!EOP


C     *******************************************
C     *** Initialize model fields from state  ***
C     *** Each model PE knows his sub-state   ***
C     *******************************************

C     *** Sea surface height (etaN) ***
      cnt = off_fields(id%etan)
      DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
            DO j=1,sNy
               DO i=1,sNx
                  cnt = cnt + 1
                  etaN(i,j,bi,bj) = state_p(cnt)
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
                     theta(i,j,k,bi,bj)= state_p(cnt)
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
                     salt(i,j,k,bi,bj)= state_p(cnt)
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
                     uVel(i,j,k,bi,bj) = state_p(cnt)
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
                     vVel(i,j,k,bi,bj)= state_p(cnt)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL DO_FIELDS_BLOCKING_EXCHANGES(myThid)


      IF ( cnt < dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'DISTRIBUTE_STATE_PDAF: cnt < dim_p:', cnt,
     &        ' > ', dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF

      IF ( cnt > dim_p ) THEN
         WRITE(msgBuf,'(A,I7,A,I7)')
     &        'DISTRIBUTE_STATE_PDAF: cnt > dim_p:', cnt,
     &        ' > ', dim_p
         CALL PRINT_ERROR( msgBuf, myThid )
         CALL ALL_PROC_DIE( 0 )
         STOP 'ABNORMAL END: S/R CONFIG_CHECK'
      ENDIF

      END
