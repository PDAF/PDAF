!> Analysis routine for DA method GLOBALTEMPLATE
!!
!! This routine computes the analysis update of the DA method.
!! Thus, for ensemble DA, it transforms the forecast ensemble
!! into the analysis ensemble.
!!
!! ADAPTING THE TEMPLATE:
!! This template contains a few typical steps of ensemble filters.
!! On this basis one can implement another DA method. Below we 
!! describe the steps that are included in this code template.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
MODULE PDAF_GLOBALTEMPLATE_analysis

CONTAINS
  SUBROUTINE PDAF_GLOBALTEMPLATE_ana(step, dim_p, dim_obs_p, dim_ens, &
       state_p, Ainv, ens_p, &
       HZ_p, HXbar_p, obs_p, &
       forget, U_prodRinvA, type_trans, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi                          ! MPI library module
    USE PDAF_timer, &                ! Routines for timings
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &          ! Routine for memory counting
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &        ! Variables for parallelization
         ONLY: mype, MPIerr, COMM_filter
    USE PDAF_analysis_utils, &       ! Utility functions
         ONLY: PDAF_subtract_rowmean

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step         !< Current time step
    INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
    INTEGER, INTENT(in ) :: dim_obs_p   !< PE-local dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
    REAL, INTENT(out)   :: state_p(dim_p)           !< on exit: PE-local forecast state
    REAL, INTENT(out)   :: Ainv(dim_ens, dim_ens)   !< on exit: weight matrix for ensemble transformation
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
    REAL, INTENT(inout) :: HZ_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
    REAL, INTENT(in) :: HXbar_p(dim_obs_p)          !< PE-local observed state
    REAL, INTENT(in) :: obs_p(dim_obs_p)            !< PE-local observation vector
    REAL, INTENT(in)    :: forget       !< Forgetting factor
    INTEGER, INTENT(in) :: type_trans   !< Type of ensemble transformation
    INTEGER, INTENT(in) :: screen       !< Verbosity flag
    INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
    INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_prodRinvA             !< Provide product R^-1 A
       
! *** local variables ***
    INTEGER :: i, col, row              ! counters
    INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
    INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
    REAL, ALLOCATABLE :: innov_p(:)     ! PE-local observation innovation
    REAL, ALLOCATABLE :: RiHZ_p(:,:)    ! Temporary matrices for analysis
    REAL, ALLOCATABLE :: tmp_Ainv(:,:)  ! Temporary storage of Ainv
    REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
    REAL, ALLOCATABLE :: rndmat(:,:)    ! Temporary random matrix


! **********************
! *** INITIALIZATION ***
! **********************

    CALL PDAF_timeit(51, 'new')

    ! Allocate arrays and count allocated memory
    IF (dim_obs_p > 0) THEN
       ALLOCATE(innov_p(dim_obs_p))
       ALLOCATE(RiHZ_p(dim_obs_p, dim_ens))
       IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p + dim_obs_p * dim_ens)
    END IF
    ALLOCATE(tmp_Ainv(dim_ens, dim_ens))
    ALLOCATE(ens_blk(maxblksize, dim_ens))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2 + maxblksize * dim_ens)


! **************************
! *** Compute innovation ***
! ***     d = y - H x    ***
! **************************

! +++ TEMPLATE:
! +++ We include this part in the template because it is generic
! +++ and will likely be needed for most ensemble DA methods

    IF (dim_obs_p > 0) THEN
       ! The innovatiopn only exists for domains with observations
     
       CALL PDAF_timeit(10, 'new')

       innov_p = obs_p - HXbar_p

       CALL PDAF_timeit(10, 'old')
    END IF

    CALL PDAF_timeit(51, 'old')


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                    +++
! +++ At this point the observations, the observed ensemble states +++
! +++ and the innovation are initialized. Now one can use these    +++
! +++ to e.g. calculate a transformation matrix for the ensemble,  +++
! +++ which is specific for each DA method.                        +++
! +++ Below we just describe some routines that can be used in     +++
! +++ the calculations.                                            +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! **********************************************
! ***  Calculate transformation matrix Ainv  ***
! **********************************************

    CALL PDAF_timeit(11, 'new')

    haveobsA: IF (dim_obs_p > 0) THEN
       ! *** The contribution of observation matrix ist only ***
       ! *** computed for domains with observations          ***

! +++ TEMPLATE: To compute the array of ensemble perturbations
! +++ one can use PDAF_subtract_rowmean which replaces the values of 
! +++ the input array HX_p by the perturbations.

       CALL PDAF_timeit(51, 'new')

       ! Subtract mean from observed ensemble: HZ = [Hx_1 ... Hx_N] T
       CALL PDAF_subtract_rowmean(dim_obs_p, dim_ens, HZ_p)

       CALL PDAF_timeit(51, 'old')

       ! ***                RiHZ = Rinv HZ                
       ! *** This is implemented as a subroutine thus that
       ! *** Rinv does not need to be allocated explicitly.

! +++ TEMPLATE: All DA methods take observation error into account.
! +++ One variant, e.g. used in ESTKF and ETKF is to multiply with
! +++ the inverse observation error covariance matrix. The routine
! +++ U_prodRinvA compute this product. This is a call-back routine.
! +++ If PDAF-OMI is used this routine is provided by OMI.

       CALL PDAF_timeit(48, 'new')
       CALL U_prodRinvA(step, dim_obs_p, dim_ens, obs_p, HZ_p, RiHZ_p)
       CALL PDAF_timeit(48, 'old')

! +++ TEMPLATE: Other possibilities to take the observation error 
! +++ into account are:
! +++ U_likelihood - compute the likelihood of an ensemble member 
! +++                (used in PF and NETF, see PDAF_pf_analysis.F90)
! +++ U_add_obs_err - add the observation error covariance matrix to some matrix
! +++                (used in EnKF and LEnKF, see PDAF_enkf_analysis_rsm.F90)

    END IF haveobsA


! +++ TEMPLATE: Note that the template does not compute the full
! +++ matrix Ainv here because the exact operation might depend on 
! +++ the DA method. 

    CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE: Note that this implementation is what is done
! +++ in the ETKF. It can be specific for the DA method.

    ! *** Initialize Ainv = (N-1) I ***
    Ainv = 0.0
    DO i = 1, dim_ens
       Ainv(i, i) = REAL(dim_ens - 1)
    END DO

! +++ TEMPLATE: The calculation of Ainv before is local for each
! +++ MPI process. One needs a global sum to get the overall value

    ! get total sum on all filter PEs
    CALL MPI_allreduce(Ainv, tmp_Ainv, dim_ens**2, &
         MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

    ! *** Complete computation of Ainv ***
    Ainv = tmp_Ainv

    CALL PDAF_timeit(51, 'old')
    CALL PDAF_timeit(11, 'old')


! +++ TEMPLATE: One might want to apply a random rotation to Ainv
! +++ PDAF_generate_rndmat provides a matrix with such random rotation
! +++ which is used in this example

    ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
    multrnd: IF (type_trans == 2) THEN
       WRITE (*,'(a, 5x,a)') 'PDAF', '--- Apply random rotation to ensemble'

       ALLOCATE(rndmat(dim_ens, dim_ens))
       IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        
       ! Initialize random matrix
       CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

       CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
            1.0, Ainv, dim_ens, rndmat, dim_ens, &
            0.0, tmp_Ainv, dim_ens)

       DEALLOCATE(rndmat)
    ELSE
       ! Non-random case
       tmp_Ainv = Ainv
    END IF multrnd


! ************************************************
! ***     Transform state ensemble             ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv.   ***
! ************************************************

    CALL PDAF_timeit(51, 'new')

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
    END IF

! +++ TEMPLATE: Finally one multiplies the ensemble (or ensemble perturbations)
! +++ with the transformation matrix Ainv. Since this overwrites the ensemble 
! +++ matrix ens_p we use here a blocked variant. A block of 'blocksize' rows
! +++ of ens_p is updated at a time so that only a small temporary matrix
! +++ (ens_blk) is required. The value of 'maxblksize' is hardcoded, but could
! +++ be changed.
! +++ The example below computes the product and adds the forecast mean state.
! +++ Thus, the update of the ensemble mean and perturbations is done together.

    CALL PDAF_timeit(21, 'new')

    ! Use block formulation for transformation
    maxblksize = 200
    IF (mype == 0 .AND. screen > 0) &
         WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize

    blocking: DO blklower = 1, dim_p, maxblksize

       blkupper = MIN(blklower + maxblksize - 1, dim_p)

       ! Store forecast ensemble
       DO col = 1, dim_ens
          ens_blk(1 : blkupper - blklower + 1, col) &
               = ens_p(blklower : blkupper, col)
       END DO

       ! Store mean forecast in ensemble matrix
       DO col = 1,dim_ens
          ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
       END DO

       !                        a  _f   f
       ! Transform ensemble:   X = X + X  Ainv
       CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
            1.0, ens_blk(1, 1), maxblksize, Ainv(1, 1), dim_ens, &
            1.0, ens_p(blklower, 1), dim_p)

    END DO blocking

    CALL PDAF_timeit(21, 'old')
    CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

! +++ TEMPLATE: One should be careful to deallocate all allocated arrays
! +++ We collect most deallocates here

    IF (dim_obs_p > 0) THEN
       DEALLOCATE(innov_p)
       DEALLOCATE(RiHZ_p)
    END IF
    DEALLOCATE(tmp_Ainv)
    DEALLOCATE(ens_blk)

! +++ TEMPLATE: Below is generic operation that is required
! +++ memory counting work

    ! Set flag that allocation was already done once (used for memory counting)
    IF (allocflag == 0) allocflag = 1

  END SUBROUTINE PDAF_GLOBALTEMPLATE_ana

END MODULE PDAF_GLOBALTEMPLATE_analysis
