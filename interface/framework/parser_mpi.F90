!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!parser_mpi.F90: TSMP-PDAF implementation of routine
!                'parser_mpi' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: parser_mpi.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !MODULE:
MODULE parser

! !DESCRIPTION:
! This module provides routine to parse command line
! arguments of different types. This version is for 
! use with MPI parallelization.
! By default, this routine uses the intrinsics 
! 'get\_command\_count' and 'get\_command\_argument' 
! that are define by the Fortran 2003 standard.
! If a compiler does not support these functions, you
! can use '-DF77' as a definition for the preprocessor.
! In this case the Fortran77 standard 'iargc()' and
! 'getarg()' are used.
!
! The module provides a generic subroutine to parse
! variables of type INTEGER, REAL, or CHARACTER
! (with length up to 100) from the command line.
!
! Usage:                      \begin{verbatim}
! SUBROUTINE PARSE(char(len=32) handle, variable)
!   The string 'handle' determines the name of    
!   the parsed variable.                          
!   Example: handle='iters' parses a variable     
!            specified on the command line by     
!            '-iters value'
!                                                 
!    Usage:                                       
!    CALL PARSE(handle, int_variable)             
!         Parses a variable of type integer       
!         whose name is given by the string       
!         handle.                                 
!                                                 
!    CALL PARSE(handle, real_variable)            
!         Parses a variable of type real          
!         whose name is given by the string       
!         handle.                                 
!                                                 
!    CALL PARSE(handle, character_variable)       
!         Parses a string variable of maxmimal    
!         length of 100 characters whose name is  
!         given by the string handle.             
!                                                 
!    CALL PARSE(handle, logical_variable)         
!         Parses a variable of type logical       
!         whose name is given by the string       
!         handle. In the command line it has      
!         to be specified as 'T' or 'F'.          
!                               \end{verbatim}
!
! !REVISION HISTORY:
! 2003-02 - Stephan Frickenhaus, Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mpi
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
