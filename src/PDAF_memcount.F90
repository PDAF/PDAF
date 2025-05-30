! Copyright (c) 2004-2025 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!
!> Module for memory alllocation
!!
!! This Module provides methods to count allocated memory.
!! 
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_memcounting

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE
  SAVE

  PUBLIC :: PDAF_memcount_ini, PDAF_memcount_define
  PUBLIC :: PDAF_memcount, PDAF_memcount_get, PDAF_memcount_get_global
  
  PRIVATE
  
  REAL, ALLOCATABLE :: mcounts(:)
  INTEGER :: wlength_i = 1
  INTEGER :: wlength_r = WORDLENGTH_REAL
  INTEGER :: wlength_d = 2
  INTEGER :: wlength_c = 4
  INTEGER :: bytespword = 4
  INTEGER :: ncnt = 0

CONTAINS
!-------------------------------------------------------------------------------
!> Initialize counters
!!
!! Subroutine to allocate and initialize 'ncounters' counters
!!
  SUBROUTINE PDAF_memcount_ini(ncounters)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ncounters  !< Number of memory counters

    ! Allocate and initialize counters
    IF (.NOT. (ALLOCATED(mcounts))) ALLOCATE(mcounts(ncounters))

    mcounts = 0.0

    ! Store number of available counters
    ncnt = ncounters

  END SUBROUTINE PDAF_memcount_ini

!-------------------------------------------------------------------------------
!> Define word length of variables
!!
!! Subroutine to define the word length of variables with type 'stortype'. 
!! In addition the length of one word in bytes can be set.
!!
!! Default lengths are (with 4 bytes per word):
!! * Integer: 1 word
!! * Real: 2 words
!! * Double: 2 words
!! * Complex: 4 words
!!
  SUBROUTINE PDAF_memcount_define(stortype, wordlength)

    IMPLICIT NONE

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: stortype  !< Type of variable
    !<    Supported are: 
    !<    (i) Integer, (r) Real, (d) Double, (c) Complex, (w) Word
    INTEGER, INTENT(IN) :: wordlength         !< Word length for chosen type


    IF (stortype == 'i') THEN
       wlength_i = wordlength
    ELSE IF (stortype == 'r') THEN
       wlength_r = wordlength
    ELSE IF (stortype == 'd') THEN
       wlength_d = wordlength
    ELSE IF (stortype == 'c') THEN
       wlength_c = wordlength
    ELSE IF (stortype == 'w') THEN
       bytespword = wordlength
    ELSE
       WRITE (*,'(a)') 'PDAF-ERROR: Storage type not supported in PDAF_MEMCOUNT!'
    END IF

  END SUBROUTINE PDAF_memcount_define

!-------------------------------------------------------------------------------
!> Count memory 
!!
!! Subroutine to count memory for the counter with index 'ID'. 
!! The allocated variable has type 'stortype' and dimension 'dim'.
!!
  SUBROUTINE PDAF_memcount(ID, stortype, dim)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ID                !< Id of the counter
    CHARACTER(len=1), INTENT(IN) :: stortype !< Type of variable
    !<    Supported are: 
    !<    (i) Integer, (r) Real, (d) Double, (c) Complex, (w) Word
    INTEGER, INTENT(in) :: dim               !< Dimension of allocated variable

!$OMP CRITICAL
    IF (stortype == 'i') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_i) * REAL(dim)
    ELSE IF (stortype == 'r') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_r) * REAL(dim)
    ELSE IF (stortype == 'd') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_d) * REAL(dim)
    ELSE IF (stortype == 'c') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_c) * REAL(dim)
    END IF
!$OMP END CRITICAL

  END SUBROUTINE PDAF_memcount

!-------------------------------------------------------------------------------
!> Reading out a memory counter
!!
!! Read out the memory count with index 'ID'. 
!! Provide size in unit 'munit'.
!!
  REAL FUNCTION PDAF_memcount_get(ID, munit)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ID             !< Id of the counter
    CHARACTER(len=1), INTENT(in) :: munit !< Unit of output
    !<    Supported are: 
    !<    (B) bytes, (K) kilo-bytes, (M) mega-bytes, (G) giga-bytes

    IF (munit == 'B' .OR. munit == 'b') THEN
       PDAF_memcount_get = REAL(bytespword) * mcounts(ID)
    ELSE IF (munit == 'k' .OR. munit == 'K') THEN
       PDAF_memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0
    ELSE IF (munit == 'm' .OR. munit == 'M') THEN
       PDAF_memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0**2
    ELSE IF (munit == 'g' .OR. munit == 'G') THEN
       PDAF_memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0**3
    ELSE
       PDAF_memcount_get = 0.0
    END IF

  END FUNCTION PDAF_memcount_get

!-------------------------------------------------------------------------------
!> Reading out a memory counter with parallelization
!!
!! This routine reads out the memory count with index 'ID'. 
!! Provide size in unit 'munit'. To get the globally counted
!! memory MPI_Allreduce is executd for the specified communicator.
!!
  REAL FUNCTION PDAF_memcount_get_global(ID, munit, comm)

    use mpi

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ID             !< Id of the counter
    CHARACTER(len=1), INTENT(in) :: munit !< Unit of output
    !<    Supported are: 
    !<    (B) bytes, (K) kilo-bytes, (M) mega-bytes, (G) giga-bytes
    INTEGER, INTENT(in) :: comm           !< Communicator

! *** Local variables
    INTEGER :: MPIerr
    REAL :: memcount_get

    ! Get Process-local memory xount
    IF (munit == 'B' .OR. munit == 'b') THEN
       memcount_get = REAL(bytespword) * mcounts(ID)
    ELSE IF (munit == 'k' .OR. munit == 'K') THEN
       memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0
    ELSE IF (munit == 'm' .OR. munit == 'M') THEN
       memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0**2
    ELSE IF (munit == 'g' .OR. munit == 'G') THEN
       memcount_get = REAL(bytespword) * mcounts(ID) / 1024.0**3
    ELSE
       memcount_get = 0.0
    END IF

    ! Get global sum of memory count
    CALL MPI_allreduce(memcount_get, PDAF_memcount_get_global, 1, &
         MPI_REALTYPE, MPI_SUM, COMM, MPIerr)

  END FUNCTION PDAF_memcount_get_global

END MODULE PDAF_memcounting
