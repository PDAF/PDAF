C     $Id$
C!BOP
C
C     !ROUTINE: init_pdaf_info - Screen output on assimilation configuration
C
C     !INTERFACE:
      SUBROUTINE init_pdaf_info()

C     !DESCRIPTION:
C     This routine performs a model-sided screen output about
C     the coniguration of the data assimilation system.
C     Using this output is optional. Most of the information
C     is also displayed by PDAF itself when it is initialized
C     in PDAF_init. Not displayed by PDAF is the assimilation
C     interval (nsteps_between_obs), which is unknown to PDAF.
C
C     !REVISION HISTORY:
C     2011-05 - Lars Nerger - Initial code extracted from init_pdaf
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,     ! Variables for assimilation       
     &     ONLY: filtertype, subtype, dim_ens,
     &     dim_state_p, nsteps_between_obs, forget

      IMPLICIT NONE

C     !CALLING SEQUENCE:
C     Called by: init_pdaf
C!EOP


C     *****************************
C     *** Initial Screen output ***
C     *****************************

      write (*,'(a,5x,a,i9)')
     &     'PDAFuser: ', 'State dimension: ', dim_state_p

      IF (filtertype == 4) THEN
         WRITE (*, '(a, 21x, a)') 'PDAFuser: ', 'Filter: ETKF'
         IF (subtype == 0) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Variant using T-matrix'
         ELSE IF (subtype == 1) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ', 
     &           '-- Variant following Hunt et al. (2007)'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &       '-- Offline mode'
         END IF

      ELSE IF (filtertype == 5) THEN
         WRITE (*, '(a, 21x, a)') 'PDAFuser: ', 'Filter: LETKF'
         IF (subtype == 0) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Variant using T-matrix'
         ELSE IF (subtype == 1) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ', 
     &           '-- Variant following Hunt et al. (2007)'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Offline mode'
         END IF

      ELSE IF (filtertype == 6) THEN
         WRITE (*, '(a, 21x, a)') 'PDAFuser: ',
     &        'Filter: ESTKF'
         IF (subtype == 0) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &       '-- Standard mode'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Offline mode'
         END IF

      ELSE IF (filtertype == 7) THEN
         WRITE (*, '(a, 21x, a)') 'PDAFuser: ',
     &        'Filter: LESTKF'
         IF (subtype == 0) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Standard mode'
         ELSE IF (subtype == 5) THEN
            WRITE (*, '(a, 6x, a)') 'PDAFuser: ',
     &           '-- Offline mode'
         END IF
      END IF

      WRITE (*, '(a, 14x, a, i5)') 'PDAFuser: ',
     &     'ensemble size:', dim_ens
      IF (subtype /= 5)
     &     WRITE (*, '(a, 6x, a, i5)') 'PDAFuser: ', 
     &     'Assimilation interval:', nsteps_between_obs
      WRITE (*, '(a, 10x, a, f5.2)') 'PDAFuser: ',
     &     'forgetting factor:', forget
      
      END SUBROUTINE init_pdaf_info
