!$Id: mod_memcount.F90 332 2019-12-30 09:37:03Z lnerger $
!> Methods to count allocated memory
!!
!! This Module provides methods to count allocated memory.
!!
!! * _memcount_ini_:
!!     Initialize memory counters
!! * _memcount_define_: 
!!     Define length of variable types
!! * _memcount_:
!!     Count allocated memory
!! * _memcount_get_:
!!     Read out memory counter
!! 
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_memcount

  IMPLICIT NONE
  SAVE

! *** Public functions
  PUBLIC :: memcount_ini, memcount_define
  PUBLIC :: memcount, memcount_get
  
  PRIVATE

  INTEGER, ALLOCATABLE :: mcounts(:)
  INTEGER :: wlength_i = 1
  INTEGER :: wlength_r = 2
  INTEGER :: wlength_d = 2
  INTEGER :: wlength_c = 4
  INTEGER :: bytespword = 4

CONTAINS
!-------------------------------------------------------------------------------
!> Initialize counters
!!
!! Subroutine to allocate and initialize 'ncounters' counters.
!!
  SUBROUTINE memcount_ini(ncounters)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ncounters  !> Number of memory counters
    
    IF (.NOT. (ALLOCATED(mcounts))) ALLOCATE(mcounts(ncounters))

    mcounts = 0

  END SUBROUTINE memcount_ini

!-------------------------------------------------------------------------------
!> Define word length of variables
!!
!! Subroutine to define the word length of variables with type 'stortype'. 
!! In addition the length of one word in bytes can be set.
!! The types and their default lengths are:
!! * (i) Integer: 1 word
!! * (r) Real: 2 words
!! * (d) Double: 2 words
!! * (c) Complex: 4 words
!! * (w) WordL 1 word
!!
!! Bytes per word: 4
!!
  SUBROUTINE memcount_define(stortype, wordlength)

    IMPLICIT NONE

! *** Arguments ***
    CHARACTER(len=1), INTENT(IN) :: stortype  !< Type of variable;
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
       WRITE (*,*) 'Storage type not supported in MEMCOUNT!'
    END IF

  END SUBROUTINE memcount_define

!-------------------------------------------------------------------------------
!> Count memory 
!!
!! Subroutine to count memory for the counter with index 'ID'. 
!! The allocated variable has type 'stortype' and dimension 'dim'.
!!
  SUBROUTINE memcount(ID, stortype, dim)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ID             !< Id of the counter
    CHARACTER(len=1), INTENT(IN) :: stortype !< Type of variable
    !<    Supported are: 
    !<    (i) Integer, (r) Real, (d) Double, (c) Complex, (w) Word
    INTEGER, INTENT(in) :: dim            !< Dimension of allocated variable


    IF (stortype == 'i') THEN
       mcounts(ID) = mcounts(ID) + wlength_i * dim
    ELSE IF (stortype == 'r') THEN
       mcounts(ID) = mcounts(ID) + wlength_r * dim
    ELSE IF (stortype == 'd') THEN
       mcounts(ID) = mcounts(ID) + wlength_d * dim
    ELSE IF (stortype == 'c') THEN
       mcounts(ID) = mcounts(ID) + wlength_c * dim
    END IF

  END SUBROUTINE memcount

!-------------------------------------------------------------------------------
!> Reading out a memory counter
!!
!! Read out the memory count with index 'ID'. 
!! Provide size in unit 'munit'.
!!
  REAL FUNCTION memcount_get(ID, munit)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ID             !< Id of the counter
    CHARACTER(len=1), INTENT(in) :: munit !< Unit of output
    !<    Supported are: 
    !<    (B) bytes, (K) kilo-bytes, (M) mega-bytes, (G) giga-bytes


    IF (munit == 'B' .OR. munit == 'b') THEN
       memcount_get = REAL(bytespword * mcounts(ID))
    ELSE IF (munit == 'k' .OR. munit == 'K') THEN
       memcount_get = REAL(bytespword * mcounts(ID)) / 1024.0
    ELSE IF (munit == 'm' .OR. munit == 'M') THEN
       memcount_get = REAL(bytespword * mcounts(ID)) / 1024.0**2
    ELSE IF (munit == 'g' .OR. munit == 'G') THEN
       memcount_get = REAL(bytespword * mcounts(ID)) / 1024.0**3
    END IF

  END FUNCTION memcount_get

END MODULE mod_memcount
