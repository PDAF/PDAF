!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case.
!!
!! __Revision history__
!! * 2026-02 - Lars Nerger - Initial code from restructuring
!! * Later revisions - see repository log
!!
MODULE mod_statevector_pdaf

  IMPLICIT NONE
  SAVE

! *** Variables to handle multiple fields in the state vector ***

!+++ Specific part for 2D tutorial model

  !< Fortran type holding the indices of model fields in the state vector
  !< This should be adapted to the fields in the state vector - it serves to give each field a name
  TYPE field_ids
     INTEGER :: NAME_OF_FIELD_1 
!     INTEGER :: NAME_OF_FIELD_1 
!     INTEGER :: ...
  END TYPE field_ids

  !< Fortran type storing size and offset of each model field in the state vector
  !< This is generic, but one could extend this type to more variables
  TYPE state_field
     INTEGER :: dim               ! size of field in state vector
     INTEGER :: off               ! offset of field in state vector
     CHARACTER(len=10) :: name    ! Name of field variable
  END TYPE state_field

!+++ End of specific part

  !---- The next variables usually do not need editing -----

  !< Type variable holding field IDs in state vector
  TYPE(field_ids) :: id

  !< number of fields in state vector
  INTEGER :: n_fields                   

  !< Vector of type variable holding dimension and offset of each field
  TYPE(state_field), ALLOCATABLE :: sfields(:)

CONTAINS


! -----------------------------------------------------------------
!> This routine initializes the array `id`
!!
!! The initialization of n_fields and of the IDs id%X
!! in this routine should be adapted to the particular state vector.
!!
  SUBROUTINE init_id(n_fields)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: n_fields

  ! Template reminder - delete when implementing functionality
  write (*,*) 'TEMPLATE mod_statevector_pdaf.F90: Define n_fields and IDs of fields in state vector!'


! Set total number of fields
!    n_fields = ?

  ! Dummy setting to ensure that code can be run
    n_fields = 1

! Set field IDs
    id%NAME_OF_FIELD_1 = 1
!    id%NAME_OF_FIELD_2 = 2


  END SUBROUTINE init_id


! -----------------------------------------------------------------
!> This routine initializes the array `sfields`
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector. It has to be adapted to
!! the particular fields used in the state vector.
!!
  SUBROUTINE init_sfields()

    ! Specific for model
!    USE mod_model, &       ! Model variables
!         ONLY: nx, ny

    IMPLICIT NONE

! *** Local variables ***
    INTEGER :: i           ! Counter


! *** Allocate ***

    ALLOCATE(sfields(n_fields))

! *****************************************************
! *** Specify sfields entry for each field variable ***
! *****************************************************

    ! Template reminder - delete when implementing functionality
    write (*,*) 'TEMPLATE mod_statevector_pdaf.F90: define field properties in variables SFIELDS!'

    ! field NAME_OF_FIELD_1
    sfields(id%NAME_OF_FIELD_1)%name = 'A'
!    sfields(id%NAME_OF_FIELD_1)%dim = ??

!+++ TEMPLATE: Dummy initialization for compilation
    sfields(id%NAME_OF_FIELD_1)%dim = 10

    ! field NAME_OF_FIELD_2
!    sfields(id%NAME_OF_FIELD_2)%name = ??
!    sfields(id%NAME_OF_FIELD_2)%dim = ??


! **************************************
! ***   Set offsets                  ***
! **************************************

! +++ The following is generic

    ! Define field offsets in state vector
    sfields(1)%off = 0
    DO i = 2, n_fields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    END DO

  END SUBROUTINE init_sfields


! -----------------------------------------------------------------
!> Initialize the state vector
!!
!! This routine is generic. Case-specific adaptions should only
!! be done in the routines init_id and init_sfields.
!!
  SUBROUTINE setup_statevector(dim_state, dim_state_p, screen)

    USE mod_parallel_pdaf, &
         ONLY: mype_model, npes_model, task_id, &
         comm_model, MPI_SUM, MPI_INTEGER, MPI_COMM_WORLD

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: dim_state    !< Global dimension of state vector
    INTEGER, INTENT(out) :: dim_state_p  !< Local dimension of state vector
    INTEGER, INTENT(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    INTEGER :: i                         ! Counters
    INTEGER :: MPIerr                    ! Error flag for MPI


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    CALL init_id(n_fields)

! *** Initialize array `sfields` ***

    CALL init_sfields()

! *** Set state vector dimension ***

    dim_state_p = SUM(sfields(:)%dim)


! *** Get global state dimension ***
    CALL MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)

! *** Write information about the state vector ***

    IF (task_id==1) THEN
       IF (mype_model==0) THEN
          WRITE (*,'(/a,2x,a)') 'model-PDAF', '*** Setup of state vector ***'
          WRITE (*,'(a,3x,a,i5)') 'model-PDAF', '--- Number of fields in state vector:', n_fields
          WRITE (*,'(a,a7,3x,a2,4x,a8,5x,a9,6x,a6)') &
               'model-PDAF','proc.','ID', 'variable', 'dimension', 'offset'
       END IF

       IF ((mype_model==0 .AND. screen<=2) .OR. screen>2) THEN
          DO i = 1, n_fields
             WRITE (*,'(a, i6,2x,i4,4x,a10,2x,i10,2x,i10)') &
                  'model-PDAF', mype_model, i, sfields(i)%name, sfields(i)%dim, sfields(i)%off
          END DO
       END IF

       IF (npes_model>1) THEN
          IF (screen>2 .OR. mype_model==0) WRITE (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'model-PDAF', 'PE', mype_model, 'process-local full state dimension: ',dim_state_p
       END IF
       IF (mype_model==0) &
            WRITE (*,'(a,2x,a,1x,i10)') 'model-PDAF', 'Global state dimension: ',dim_state
    END IF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

  END SUBROUTINE setup_statevector

END MODULE mod_statevector_pdaf
