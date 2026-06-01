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
         only: fieldA_p, fieldB_p, coords_x_p, coords_y_p
    use model_parallel_mod, &     ! Model parallelzation variables
         only: myproc_world

    implicit none

! **********************
! *** Finalization   ***
! **********************

    ! Dellocate model arrays
    deallocate(fieldA_p, fieldB_p)
    deallocate(coords_x_p, coords_y_p)

    if (myproc_world==0) then
       write (*, '(/10x, a)') '+++++ PDAF tutorial model completed +++++'
    end if

  end subroutine postprocess

end module model_post_mod
