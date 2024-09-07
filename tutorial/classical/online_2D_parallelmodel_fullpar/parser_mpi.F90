!> Command line parser
!!
!! This module provides routine to parse command line
!! arguments of different types. This version is for 
!! use with MPI parallelization.
!!
!! By default, this routine uses the intrinsics 
!! 'get_command_count' and 'get_command_argument' 
!! that are defined by the Fortran 2003 standard.
!! If a compiler does not support these functions, you
!! can use '-DF77' as a definition for the preprocessor.
!! In this case the Fortran77 standard 'iargc()' and
!! 'getarg()' are used.
!!
!! The module provides a generic subroutine to parse
!! variables of type INTEGER, REAL, or CHARACTER
!! (with length up to 100) from the command line.
!!
!! Usage: 
!! SUBROUTINE PARSE(char(len=32) handle, variable)
!!   The string 'handle' determines the name of    
!!   the parsed variable.                          
!!   Example: handle='iters' parses a variable     
!!            specified on the command line by     
!!            '-iters value'
!!                                                 
!!    Usage:                                       
!!    CALL PARSE(handle, int_variable)             
!!         Parses a variable of type integer       
!!         whose name is given by the string       
!!         handle.                                 
!!                                                 
!!    CALL PARSE(handle, real_variable)            
!!         Parses a variable of type real          
!!         whose name is given by the string       
!!         handle.                                 
!!                                                 
!!    CALL PARSE(handle, character_variable)       
!!         Parses a string variable of maxmimal    
!!         length of 100 characters whose name is  
!!         given by the string handle.             
!!                                                 
!!    CALL PARSE(handle, logical_variable)         
!!         Parses a variable of type logical       
!!         whose name is given by the string       
!!         handle. In the command line it has      
!!         to be specified as 'T' or 'F'.          
!!
!! __Revision history:__
!! * 2003-02 - Stephan Frickenhaus, Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE parser

  use mpi 
  USE mod_parallel_model, &
    ONLY: abort_parallel

  IMPLICIT NONE
  SAVE
  
! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: parse
  CHARACTER(len=32), PUBLIC :: handle  ! handle for command line parser
!EOP

  PRIVATE
  CHARACTER(len=100) :: str1, str2 
  INTEGER :: i   
  INTEGER :: mype, MPIerr
!   INTEGER,EXTERNAL :: iargc


  ! *** define interface ***
  INTERFACE parse
    MODULE PROCEDURE parse_int
    MODULE PROCEDURE parse_real
    MODULE PROCEDURE parse_string
    MODULE PROCEDURE parse_logical
  END INTERFACE

CONTAINS
  SUBROUTINE parse_int(handle, intvalue)
! ******************************
! *** Parse an integer value ***
! ******************************

! *** subroutine arguments ***    
    CHARACTER(len=32), INTENT(in) :: handle
    INTEGER,INTENT(inout) :: intvalue

! *** local variables ***
    CHARACTER(len=32) :: string
    INTEGER :: parsed_int
    LOGICAL :: modified

! *** Initialization ***
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // TRIM(handle)
    modified = .FALSE.
    
! *** Parsing ***
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

! *** Finalize ***
    IF (modified) THEN
       intvalue = parsed_int
!        IF (mype == 0) WRITE (*, '(2x, a, a, a, i)') &
       IF (mype == 0) WRITE (*, '(2x, a, a, a, i10)') &
            'PARSER: ', TRIM(handle), '=', parsed_int
    END IF
  END SUBROUTINE parse_int

  SUBROUTINE parse_real(handle, realvalue)
! **************************
! *** Parse a real value ***
! **************************

! *** function arguments ***    
    CHARACTER(len=32), INTENT(in) :: handle
    REAL, INTENT(inout) :: realvalue

! *** local variables ***
    CHARACTER(len=32) :: string
    REAL :: parsed_real
    LOGICAL :: modified

! *** Initialize ***
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // TRIM(handle)
    modified = .FALSE.

! *** Parsing ***
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

! *** Finalize ***
    IF (modified) THEN
       realvalue = parsed_real
       IF (mype == 0) WRITE (*, '(2x, a, a, a, es12.4)') &
            'PARSER: ', TRIM(handle), '=', parsed_real
    END IF
  END SUBROUTINE parse_real


  SUBROUTINE parse_string(handle, charvalue)
! **********************
! *** Parse a string ***
! **********************

! *** function arguments ***    
    CHARACTER(len=32), INTENT(in) :: handle
    CHARACTER(len=*), INTENT(inout) :: charvalue

! *** local variables ***
    CHARACTER(len=100) :: string
    CHARACTER(len=100) :: parsed_string
    CHARACTER(len=110) :: str1_check
    CHARACTER(len=110) :: str2_check
    LOGICAL :: modified

! *** Initialize ***
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // TRIM(handle)
    modified = .FALSE.
    
! *** Parsing ***
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

          ! Add check for inadmissible strings longer than 100
          ! characters
          CALL get_command_argument(i, str1_check)
          CALL get_command_argument(i+1, str2_check)
          IF (mype == 0) THEN
             IF (.NOT. TRIM(str2_check) == TRIM(str2)) THEN
                WRITE (*,'(2x, a)') "PARSER: ERROR, command line input too long."
                WRITE (*,'(2x, a, 1x, a)') "called handle=", TRIM(string)
                WRITE (*,'(2x, a, 1x, a)') "parsed handle=", TRIM(str1)
                WRITE (*,'(2x, a, 1x, a)') "parsed input(cut)=", TRIM(str2)
                call abort_parallel()
             END IF
          END IF


#endif
          IF (str1 == TRIM(string)) THEN
             ! Format specifier is needed for reading paths.  Using
             ! `*` as format specifier, reading stops at a `/`
             READ(str2, '(a)') parsed_string
             modified = .TRUE.
          END IF
       ENDDO
    ENDIF

! *** Finalize ***
    IF (modified) THEN
       charvalue = parsed_string
       IF (mype == 0) WRITE (*, '(2x, a, a, a, a)') &
           'PARSER: ', TRIM(handle), '= ', TRIM(parsed_string)
    END IF

  END SUBROUTINE parse_string

  SUBROUTINE parse_logical(handle, logvalue)
! ******************************
! *** Parse an logical value ***
! ******************************

! *** subroutine arguments ***    
    CHARACTER(len=32), INTENT(in) :: handle
    LOGICAL, INTENT(inout) :: logvalue

! *** local variables ***
    CHARACTER(len=32) :: string
    LOGICAL :: parsed_log
    LOGICAL :: modified

! *** Initialization ***
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // TRIM(handle)
    modified = .FALSE.
    
! *** Parsing ***
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

! *** Finalize ***
    IF (modified) THEN
       logvalue = parsed_log
       IF (mype == 0) WRITE (*, '(2x, a, a, a, l1)') &
            'PARSER: ', TRIM(handle), '=', parsed_log
    END IF
  END SUBROUTINE parse_logical

END MODULE parser
