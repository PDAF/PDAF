!$Id: mod_memcount.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_memcount

! !DESCRIPTION:
! This Module provides methods to count allocated memory.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: memcount_ini, memcount_define
  PUBLIC :: memcount, memcount_get
!EOP
  
  PRIVATE
  
  REAL, ALLOCATABLE :: mcounts(:)
  INTEGER :: wlength_i = 1
  INTEGER :: wlength_r = 2
  INTEGER :: wlength_d = 2
  INTEGER :: wlength_c = 4
  INTEGER :: bytespword = 4

CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: memcount_ini - Initialize counters
!
! !INTERFACE: memcount_ini()
  SUBROUTINE memcount_ini(ncounters)

! !DESCRIPTION:
! Subroutine to allocate and initialize 'ncounters' counters.\\
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: ncounters  ! Number of memory counters
!EOP
    
    IF (.NOT. (ALLOCATED(mcounts))) ALLOCATE(mcounts(ncounters))

    mcounts = 0.0

  END SUBROUTINE memcount_ini

!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: memcount_define - Define word length of variables
!
! !INTERFACE: memcount_define()
  SUBROUTINE memcount_define(stortype, wordlength)

! !DESCRIPTION:
! Subroutine to define the word length of variables with type 'stortype'. 
! In addition the length of one word in bytes can be set.
! Default lengths are:\\
! Integer: 1 word
! - Real: 2 words
! - Double: 2 words
! - Complex: 4 words\\
! Bytes per word: 4
!
! !USES:      
    IMPLICIT NONE

! !ARGUMENTS:
    CHARACTER(len=1), INTENT(IN) :: stortype  ! Type of variable
    !    Supported are: 
    !    (i) Integer, (r) Real, (d) Double, (c) Complex, (w) Word
    INTEGER, INTENT(IN) :: wordlength         ! Word length for chosen type
!EOP

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
!BOP
!
! !ROUTINE: memcount - Count memory 
!
! !INTERFACE: memcount()
  SUBROUTINE memcount(ID, stortype, dim)

! !DESCRIPTION:
! Subroutine to count memory for the counter with index 'ID'. 
! The allocated variable has type 'stortype' and dimension 'dim'.

! !USES:
    IMPLICIT NONE

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: ID             ! Id of the counter
    CHARACTER(len=1), INTENT(IN) :: stortype ! Type of variable
    !    Supported are: 
    !    (i) Integer, (r) Real, (d) Double, (c) Complex, (w) Word
    INTEGER, INTENT(in) :: dim            ! Dimension of allocated variable
!EOP

    IF (stortype == 'i') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_i) * REAL(dim)
    ELSE IF (stortype == 'r') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_r) * REAL(dim)
    ELSE IF (stortype == 'd') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_d) * REAL(dim)
    ELSE IF (stortype == 'c') THEN
       mcounts(ID) = mcounts(ID) + REAL(wlength_c) * REAL(dim)
    END IF

  END SUBROUTINE memcount

!-------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: memcount_get - Reading out a memory counter
!
! !INTERFACE: memcount_get()
  REAL FUNCTION memcount_get(ID, munit)

! !DESCRIPTION:
! Read out the memory count with index 'ID'. 
! Provide size in unit 'munit'.

! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: ID             ! Id of the counter
    CHARACTER(len=1), INTENT(in) :: munit ! Unit of output
    !    Supported are: 
    !    (B) bytes, (K) kilo-bytes, (M) mega-bytes, (G) giga-bytes
!EOP

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

  END FUNCTION memcount_get

END MODULE mod_memcount
