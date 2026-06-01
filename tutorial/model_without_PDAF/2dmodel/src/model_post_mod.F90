!>  Postprocessing of model
!!
!! The routine does the postprocessing of the model run. Here we simply deallocate
!! model arrays.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_post_mod

contains

  subroutine postprocess()

    use model_mod, &              ! Model variables
         only: fieldA, fieldB, coords_x, coords_y

    implicit none


! **********************
! *** Finalization   ***
! **********************

    ! Dellocate model arrays
    deallocate(fieldA, fieldB)
    deallocate(coords_x, coords_y)

    write (*, '(/10x, a)') '+++++ PDAF tutorial model completed +++++'

  end subroutine postprocess

end module model_post_mod
