C     $Id$
C     !BOP
C
CC    !DESCRIPTION:
C  This module provides routine to parse command line
C  arguments of different types. This version is for
C  use with MPI parallelization.
C  By default, this routine uses the intrinsics
C  'get\_command\_count' and 'get\_command\_argument'
C  that are define by the Fortran 2003 standard.
C  If a compiler does not support these functions, you
C  can use '-DF77' as a definition for the preprocessor.
C  In this case the Fortran77 standard 'iargc()' and
C  'getarg()' are used.
C
C  The module provides a generic subroutine to parse
C  variables of type INTEGER, REAL, or CHARACTER
C  (with length up to 100) from the command line.
C
C  Usage:                      \begin{verbatim}
C  SUBROUTINE PARSE(char(len=32) handle, variable)
C    The string 'handle' determines the name of
C    the parsed variable.
C    Example: handle='iters' parses a variable
C             specified on the command line by
C             '-iters value'
C
C     Usage:
C     CALL PARSE(handle, int_variable)
C          Parses a variable of type integer
C          whose name is given by the string
C          handle.
C
C     CALL PARSE(handle, real_variable)
C          Parses a variable of type real
C          whose name is given by the string
C          handle.
C
C     CALL PARSE(handle, character_variable)
C          Parses a string variable of maxmimal
C          length of 100 characters whose name is
C          given by the string handle.
C
C     CALL PARSE(handle, logical_variable)
C          Parses a variable of type logical
C          whose name is given by the string
C          handle. In the command line it has
C          to be specified as 'T' or 'F'.
C                                \end{verbatim}
C
C  !REVISION HISTORY:
C  2003-02 - Stephan Frickenhaus, Lars Nerger - Initial code
C  Later revisions - see svn log
C
C     !MODULE:
      MODULE parser
      
C     !USES:

#ifdef ALLOW_USE_MPI
      USE mpi
#else
These lines are here on purpose to STOP the compilation,
because online-PDAF needs MPI
#endif
      
      IMPLICIT NONE
      SAVE


C     !PUBLIC MEMBER FUNCTIONS:
      PUBLIC :: parse
      CHARACTER(len=32), PUBLIC :: handle ! handle for command line parser

      PRIVATE
      CHARACTER(len=100) :: str1, str2
      INTEGER :: i
      INTEGER :: mype, MPIerr
C     INTEGER,EXTERNAL :: iargc


C     *** define interface ***
      INTERFACE parse
        MODULE PROCEDURE parse_int
        MODULE PROCEDURE parse_real
        MODULE PROCEDURE parse_string
        MODULE PROCEDURE parse_logical
      END INTERFACE

      CONTAINS
      SUBROUTINE parse_int(handle, intvalue)
C     ******************************
C     *** Parse an integer value ***
C     ******************************

C     *** subroutine arguments ***
      CHARACTER(len=32), INTENT(in) :: handle
      INTEGER,INTENT(inout) :: intvalue

C     *** local variables ***
      CHARACTER(len=32) :: string
      INTEGER :: parsed_int
      LOGICAL :: modified

C     *** Initialization ***
      CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

      string = '-' // TRIM(handle)
      modified = .FALSE.

C     *** Parsing ***
#ifdef F77
      write (*,*) 'PARSE for F77!!!!!!!!!!!!!!!'
      IF (iargc() > 0) THEN
        DO i = 1, iargc() - 1
          CALL getarg(i, str1)
          CALL getarg(i + 1, str2)
#else
      IF (command_argument_count() > 0) THEN
        DO i = 1, command_argument_count() - 1
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
#endif
          IF (str1 == TRIM(string)) THEN
            READ(str2, *) parsed_int
            modified = .TRUE.
          END IF
        ENDDO
      ENDIF

C     *** Finalize ***
      IF (modified) THEN
        intvalue = parsed_int
C     IF (mype == 0) WRITE (*, '(2x, a, a, a, i)') &
        IF (mype == 0) WRITE (*, '(2x, a, a, a, i10)')
     &            'PARSER: ', TRIM(handle), '=', parsed_int
      END IF

      END SUBROUTINE parse_int

      SUBROUTINE parse_real(handle, realvalue)
C     **************************
C     *** Parse a real value ***
C     **************************

C     *** function arguments ***
      CHARACTER(len=32), INTENT(in) :: handle
      REAL(kind=8), INTENT(inout) :: realvalue

C     *** local variables ***
      CHARACTER(len=32) :: string
      REAL(kind=8) :: parsed_real
      LOGICAL :: modified

C     *** Initialize ***
      CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

      string = '-' // TRIM(handle)
      modified = .FALSE.

C     *** Parsing ***
#ifdef F77
      IF (iargc() > 0) THEN
        DO i = 1, iargc() - 1
          CALL getarg(i, str1)
          CALL getarg(i + 1, str2)
#else
      IF (command_argument_count() > 0) THEN
        DO i = 1, command_argument_count() - 1
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
#endif
          IF (str1 == TRIM(string)) THEN
            READ(str2, *) parsed_real
            modified = .TRUE.
          END IF
        ENDDO
      ENDIF

C     *** Finalize ***
      IF (modified) THEN
        realvalue = parsed_real
        IF (mype == 0) WRITE (*, '(2x, a, a, a, es12.4)')
     &        'PARSER: ', TRIM(handle), '=', parsed_real
      END IF
      END SUBROUTINE parse_real


      SUBROUTINE parse_string(handle, charvalue)
C     **********************
C     *** Parse a string ***
C     **********************

C     *** function arguments ***
      CHARACTER(len=32), INTENT(in) :: handle
      CHARACTER(len=*), INTENT(inout) :: charvalue

C     *** local variables ***
      CHARACTER(len=100) :: string
      CHARACTER(len=100) :: parsed_string
      LOGICAL :: modified

C     *** Initialize ***
      CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

      string = '-' // TRIM(handle)
      modified = .FALSE.

C     *** Parsing ***
#ifdef F77
      IF (iargc() > 0) THEN
        DO i = 1, iargc() - 1
          CALL getarg(i, str1)
          CALL getarg(i + 1, str2)
#else
      IF (command_argument_count() > 0) THEN
        DO i = 1, command_argument_count() - 1
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
#endif
          IF (str1 == TRIM(string)) THEN
            READ(str2, *) parsed_string
            modified = .TRUE.
          END IF
        ENDDO
      ENDIF

C     *** Finalize ***
      IF (modified) THEN
        charvalue = parsed_string
        IF (mype == 0) WRITE (*, '(2x, a, a, a, a)')
     &        'PARSER: ', TRIM(handle), '= ', TRIM(parsed_string)
      END IF

      END SUBROUTINE parse_string

      SUBROUTINE parse_logical(handle, logvalue)
C     ******************************
C     *** Parse an logical value ***
C     ******************************

C     *** subroutine arguments ***
      CHARACTER(len=32), INTENT(in) :: handle
      LOGICAL, INTENT(inout) :: logvalue

C     *** local variables ***
      CHARACTER(len=32) :: string
      LOGICAL :: parsed_log
      LOGICAL :: modified

C     *** Initialization ***
      CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

      string = '-' // TRIM(handle)
      modified = .FALSE.

C  *** Parsing ***
#ifdef F77
      IF (iargc() > 0) THEN
        DO i = 1, iargc() - 1
          CALL getarg(i, str1)
          CALL getarg(i + 1, str2)
#else
      IF (command_argument_count() > 0) THEN
        DO i = 1, command_argument_count() - 1
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
#endif
          IF (str1 == TRIM(string)) THEN
            READ(str2, *) parsed_log
               modified = .TRUE.
          END IF
        ENDDO
      ENDIF

C  *** Finalize ***
      IF (modified) THEN
        logvalue = parsed_log
        IF (mype == 0) WRITE (*, '(2x, a, a, a, l)')
     &        'PARSER: ', TRIM(handle), '=', parsed_log
      END IF
      END SUBROUTINE parse_logical

      END MODULE parser
