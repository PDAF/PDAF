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
module mod_statevector_pdaf

  implicit none
  save

! *** Variables to handle multiple fields in the state vector ***

!+++ Specific part for 2D tutorial model

  !< Fortran type holding the indices of model fields in the state vector
  !< This should be adapted to the fields in the state vector - it serves to give each field a name
  type field_ids
     integer :: fieldA 
     integer :: fieldB
  end type field_ids

  !< Fortran type storing size and offset of each model field in the state vector
  !< This is generic, but one could extend this type to more variables
  type state_field
     integer :: dim               ! size of field in state vector
     integer :: off               ! offset of field in state vector
     character(len=10) :: name    ! Name of field variable
  end type state_field

!+++ End of specific part

  !---- The next variables usually do not need editing -----

  !< Type variable holding field IDs in state vector
  type(field_ids) :: id

  !< number of fields in state vector
  integer :: n_fields                   

  !< Vector of type variable holding dimension and offset of each field
  type(state_field), allocatable :: sfields(:)

contains


! -----------------------------------------------------------------
!> This routine initializes the array `id`
!!
!! The initialization of n_fields and of the IDs id%X
!! in this routine should be adapted to the particular state vector.
!!
  subroutine init_id(n_fields)

    implicit none

! *** Arguments ***
    integer, intent(out) :: n_fields

!+++ Specific part for 2D tutorial model

! Set total number of fields
    n_fields = 2

! Set field IDs
    id%fieldA = 1
    id%fieldB = 2

!+++ End of specific part

  end subroutine init_id


! -----------------------------------------------------------------
!> This routine initializes the array `sfields`
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector. It has to be adapted to
!! the particular fields used in the state vector.
!!
  subroutine init_sfields()

    ! Specific for model
    use mod_model, &       ! Model variables
         only: nx, ny

    implicit none

! *** Local variables ***
    integer :: i           ! Counter


! *** Allocate ***

    allocate(sfields(n_fields))

! *****************************************************
! *** Specify sfields entry for each field variable ***
! *****************************************************

!+++ Specific part for 2D tutorial model

    ! fieldA
    sfields(id%fieldA)%name = 'A'

    ! fieldB
    sfields(id%fieldB)%name = 'B'


! **************************************
! ***   Set dimensions and offsets   ***
! **************************************

    ! Set field dimensions
    do i = 1, n_fields
       sfields(i)%dim = nx * ny
    end do

! +++ The following is generic

    ! Define field offsets in state vector
    sfields(1)%off = 0
    do i = 2, n_fields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do

  end subroutine init_sfields


! -----------------------------------------------------------------
!> Initialize the state vector
!!
!! This routine is generic. Case-specific adaptions should only
!! be done in the routines init_id and init_sfields.
!!
  subroutine setup_statevector(dim_state, dim_state_p, screen)

    use mod_parallel_pdaf, &
         only: mype_model, npes_model, task_id, &
         comm_model, MPI_SUM, MPI_INTEGER, MPI_COMM_WORLD

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state    !< Global dimension of state vector
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector
    integer, intent(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    integer :: i                         ! Counters
    integer :: MPIerr                    ! Error flag for MPI


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    call init_id(n_fields)

! *** Initialize array `sfields` ***

    call init_sfields()

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)


! *** Get global state dimension ***
    call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)

! *** Write information about the state vector ***

    if (task_id==1) then
       if (mype_model==0) then
          write (*,'(/a,2x,a)') 'model-PDAF', '*** Setup of state vector ***'
          write (*,'(a,3x,a,i5)') 'model-PDAF', '--- Number of fields in state vector:', n_fields
          write (*,'(a,a7,3x,a2,4x,a8,5x,a9,6x,a6)') &
               'model-PDAF','proc.','ID', 'variable', 'dimension', 'offset'
       end if

       if ((mype_model==0 .and. screen<=2) .or. screen>2) then
          do i = 1, n_fields
             write (*,'(a, i6,2x,i4,4x,a10,2x,i10,2x,i10)') &
                  'model-PDAF', mype_model, i, sfields(i)%name, sfields(i)%dim, sfields(i)%off
          end do
       end if

       if (npes_model>1) then
          if (screen>2 .or. mype_model==0) write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'model-PDAF', 'PE', mype_model, 'process-local full state dimension: ',dim_state_p
       end if
       if (mype_model==0) &
            write (*,'(a,2x,a,1x,i10)') 'model-PDAF', 'Global state dimension: ',dim_state
    end if
    call MPI_Barrier(MPI_COMM_WORLD, MPIerr)

  end subroutine setup_statevector

end module mod_statevector_pdaf
