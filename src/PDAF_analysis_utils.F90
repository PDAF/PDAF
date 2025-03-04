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
!> Utility routines for analysis step
!!
MODULE PDAF_analysis_utils

  ! Array for observation statistics in local analysis
  INTEGER :: obsstats(4)           ! PE-local statistics
  ! obsstats(1): Local domains with observations
  ! obsstats(2): Local domains without observations
  ! obsstats(3): Sum of all available observations for all domains
  ! obsstats(4): Maximum number of observations over all domains


CONTAINS

!-------------------------------------------------------------------------------
!> Generate random matrix with special properties
!!
!! Generate a transformation matrix OMEGA for
!! the generation and transformation of the 
!! ensemble in the SEIK and LSEIK filter.
!! Generated is a uniform orthogonal matrix OMEGA
!! with R columns orthonormal in $R^{r+1}$
!! and orthogonal to (1,...,1)' by iteratively 
!! applying the Householder matrix onto random 
!! vectors distributed uniformly on the unit sphere.
!!
!! This version initializes at each iteration step
!! the whole Householder matrix and subsequently
!! computes Omega using GEMM from BLAS. All fields are 
!! allocated once at their maximum required size.
!! (On SGI O2K this is about a factor of 2.5 faster
!! than the version applying BLAS DDOT, but requires
!! more memory.)
!!
!! For Omegatype=0 a deterministic Omega is computed
!! where the Housholder matrix of (1,...,1)' is operated
!! on an identity matrix.
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2002-01 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_generate_rndmat(dim, rndmat, mattype)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: seedset, new_seedset

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim         ! Size of matrix rndmat
  REAL, INTENT(out)   :: rndmat(dim, dim) ! Matrix
  INTEGER, INTENT(in) :: mattype     ! Select type of random matrix:
                                     !   (1) orthonormal random matrix
                                     !   (2) orthonormal with eigenvector (1,...,1)^T

!  *** local variables ***
  INTEGER :: iter, col, row          ! counters
  INTEGER :: i, j, k                 ! counters
  INTEGER :: dimrnd                  ! Size of random matrix to be generation at first part
  INTEGER, SAVE :: iseed(4)          ! seed array for random number routine
  REAL :: norm                       ! norm of random vector
  INTEGER :: pflag                   ! pointer flag
  INTEGER, SAVE :: first = 1         ! flag for init of random number seed
  REAL :: rndval                     ! temporary value for init of Householder matrix
  REAL, ALLOCATABLE :: rndvec(:)     ! vector of random numbers
  REAL, ALLOCATABLE :: h_rndvec(:)   ! vector of random numbers
  REAL, ALLOCATABLE :: house(:,:)    ! Householder matrix
  REAL, ALLOCATABLE :: matUBB(:,:)   ! Temporary matrix
  REAL, POINTER :: mat_iter(:,:)     ! Pointer to temporary random array
  REAL, POINTER :: mat_itermin1(:,:) ! Pointer to temporary random array
  REAL, POINTER :: matU(:,:)         ! Pointer to temporary array
  REAL, POINTER :: matUB(:,:)        ! Pointer to temporary array
  REAL, POINTER :: matB(:,:)         ! Pointer to temporary array
  REAL, ALLOCATABLE, TARGET :: temp1(:,:)  ! Target array
  REAL, ALLOCATABLE, TARGET :: temp2(:,:)  ! Target array


! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine size of matrix build through householder reflections
  randOmega: IF (mattype == 1) THEN
     ! Random orthonormal matrix
     dimrnd = dim
  ELSE
     ! Random orthonormal matrix with eigenvector (1,...,1)^T
     dimrnd = dim - 1
  END IF randOmega


! ******************************************
! *** Generate orthonormal random matrix ***
! ******************************************

  ! allocate fields
  ALLOCATE(rndvec(dim))
  ALLOCATE(house(dim + 1, dim))
  ALLOCATE(temp1(dim, dim), temp2(dim, dim))

  ! set pointers
  mat_itermin1 => temp1
  mat_iter     => temp2
  pflag = 0

  ! Initialized seed for random number routine
  IF (first == 1 .OR. new_seedset) THEN
     IF (seedset == 2) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 3) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 4) THEN
        iseed(1)=1
        iseed(2)=6
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 5) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 6) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 7) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 8) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 9) THEN
        iseed(1)=3
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 10) THEN
        iseed(1)=3
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 11) THEN
        iseed(1)=13
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 12) THEN
        iseed(1)=13
        iseed(2)=11
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 13) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 14) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=17
        iseed(4)=11
     ELSE IF (seedset == 15) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 16) THEN
        iseed(1)=15
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 17) THEN
        iseed(1)=15
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 18) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 19) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=17
     ELSE IF (seedset == 20) THEN
        iseed(1)=15
        iseed(2)=15
        iseed(3)=47
        iseed(4)=17
     ELSE
        ! Standard seed
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
     END IF
     first = 2
     new_seedset = .FALSE.
  END IF


! *** First step of iteration       ***  
! *** Determine mat_iter for iter=1 ***

  ! Get random number [-1,1]
  CALL larnvTYPE(2, iseed, 1, rndvec(1))
  
  IF (rndvec(1) >= 0.0) THEN
     mat_itermin1(1, 1) = +1.0
  ELSE
     mat_itermin1(1, 1) = -1.0
  END IF

! *** Iteration ***

  iteration: DO iter = 2, dimrnd

! Initialize new random vector
      
     ! Get random vector of dimension DIM (elements in [-1,1])
     CALL larnvTYPE(2, iseed, iter, rndvec(1:iter))

     ! Normalize random vector
     norm = 0.0
     DO col = 1, iter
        norm = norm + rndvec(col)**2
     END DO
     norm = SQRT(norm)
        
     DO col = 1, iter
        rndvec(col) = rndvec(col) / norm
     END DO

! Compute Householder matrix

     ! First ITER-1 rows
     rndval = 1.0 / (ABS(rndvec(iter)) + 1.0)
     housecol: DO col = 1, iter - 1
        houserow: DO row = 1,iter - 1
           house(row, col) = - rndvec(row) * rndvec(col) * rndval
        END DO houserow
     END DO housecol
        
     DO col = 1, iter - 1
        house(col, col) = house(col, col) + 1.0
     END DO

     ! Last row
     housecol2: DO col = 1, iter - 1
        house(iter, col) = - (rndvec(iter) + SIGN(1.0, rndvec(iter))) &
             * rndvec(col) * rndval
     END DO housecol2

! Compute matrix on this iteration stage

     ! First iter-1 columns
     CALL gemmTYPE ('n', 'n', iter, iter - 1, iter - 1, &
          1.0, house, dim + 1, mat_itermin1, dim, &
          0.0, mat_iter, dim)

     ! Final column
     DO row = 1, iter
        mat_iter(row, iter) = rndvec(row)
     END DO

! Adjust pointers to temporal OMEGA fields

     IF (pflag == 0) THEN
        mat_itermin1 => temp2
        mat_iter     => temp1
        pflag = 1
     ELSE IF (pflag == 1) THEN
        mat_itermin1 => temp1
        mat_iter     => temp2
        pflag = 0
     END IF

  END DO iteration


! ****************************************************
! *** Ensure eigenvector (1,...1,)^T for mattype=2 ***
! ****************************************************

  mattype2: IF (mattype == 1) THEN

     ! *** Generation of random matrix completed for mattype=1
     rndmat = mat_itermin1

  ELSE mattype2

     ! *** Complete generation of random matrix with eigenvector
     ! *** (1,...,1)^T by transformation with a basis that
     ! *** includes (1,...,1)^T. (We follow the description 
     ! *** Sakov and Oke, MWR 136, 1042 (2008)).

     NULLIFY(mat_iter, mat_itermin1)

     ALLOCATE(h_rndvec(dim))

! *** Complete initialization of random matrix with eigenvector ***
! *** (1,...,1)^T in the basis that includes (1,...,1)^T        ***

     IF (pflag == 0) THEN
        matU   => temp1
        matUB => temp2
     ELSE
        matU   => temp2
        matUB => temp1
     END if

     matUB(:,:) = 0.0
     matUB(1,1) = 1.0
     DO col = 2, dim
        DO row = 2, dim
           matUB(row, col) = matU(row - 1, col - 1)
        END DO
     END DO
     NULLIFY(matU)

! *** Generate orthonormal basis including (1,...,1)^T as leading vector ***
! *** We again use houesholder reflections.                              ***

     IF (pflag == 0) THEN
        matB => temp1
     ELSE
        matB => temp2
     END IF

     ! First column
     DO row = 1, dim
        matB(row, 1) = 1.0 / SQRT(REAL(dim))
     END DO

     ! columns 2 to dim
     buildB: DO col = 2, dim

        ! Get random vector of dimension DIM (elements in [0,1])
        CALL larnvTYPE(1, iseed, dim, rndvec)

        loopcols: DO i = 1, col - 1
           DO j = 1, dim
              DO k = 1, dim
                 house(k, j) = - matB(k,i) * matB(j,i)
              END DO
           END DO
           DO j = 1, dim
              house(j, j) = house(j, j) + 1.0
           END DO

           ! Apply house to random vector
           CALL gemvTYPE ('n', dim, dim, &
                1.0, house, dim+1, rndvec, 1, &
                0.0, h_rndvec, 1)
           rndvec = h_rndvec

        END DO loopcols

        ! Normalize vector
        norm = 0.0
        DO i = 1, iter
           norm = norm + h_rndvec(i)**2
        END DO
        norm = SQRT(norm)
        
        DO i = 1, iter
           h_rndvec(i) = h_rndvec(i) / norm
        END DO

        ! Inialize column of matB
        matB(:, col) = h_rndvec

     END DO buildB


! *** Final step: Transform random matrix  ***
! *** rndmat = matB matUB matB^T  ***

     ALLOCATE(matUBB(dim, dim))

     ! matUB * matB^T
     CALL gemmTYPE ('n', 't', dim, dim, dim, &
          1.0, matUB, dim, matB, dim, &
          0.0, matUBB, dim)

     ! matB * matUB * matB^T
     CALL gemmTYPE ('n', 'n', dim, dim, dim, &
          1.0, matB, dim, matUBB, dim, &
          0.0, rndmat, dim)

! *** CLEAN UP ***

     NULLIFY(matUB, matB)
     DEALLOCATE(matUBB)
     DEALLOCATE(h_rndvec)

  END IF mattype2


! ************************
! *** General clean up ***
! ************************

  DEALLOCATE(temp1, temp2)
  DEALLOCATE(rndvec, house)

END SUBROUTINE PDAF_generate_rndmat


!-------------------------------------------------------------------------------
!> Print statistics on local analysis domains
!!
  SUBROUTINE PDAF_print_domain_stats(n_domains_p)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, COMM_filter, MPIerr

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: n_domains_p   ! Number of PE-local analysis domains

! *** Local Variables ***
    INTEGER :: n_domains_stats(4)        ! Gobal statistics for number of analysis domains


    IF (npes_filter>1) THEN
       CALL MPI_Reduce(n_domains_p, n_domains_stats(1), 1, MPI_INTEGER, MPI_MIN, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(n_domains_p, n_domains_stats(2), 1, MPI_INTEGER, MPI_MAX, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(n_domains_p, n_domains_stats(3), 1, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
       IF (mype == 0) THEN
          WRITE (*, '(a, 5x, a, i9, 1x, i9, 1x, f11.1)') &
               'PDAF', '--- local analysis domains (min/max/avg):', n_domains_stats(1:2), &
               REAL(n_domains_stats(3)) / REAL(npes_filter)
       END IF
    ELSE
       ! This is a work around for working with nullmpi.F90
       IF (mype == 0) THEN
          WRITE (*, '(a, 5x, a, i9)') &
               'PDAF', '--- local analysis domains:', n_domains_p
       END IF
    END IF

  END SUBROUTINE PDAF_print_domain_stats


!-------------------------------------------------------------------------------

!> Initialize local observation statistics
!!
  SUBROUTINE PDAF_init_local_obsstats()

    IMPLICIT NONE

! *** Set obsstats to zero ***

    obsstats = 0

  END SUBROUTINE PDAF_init_local_obsstats


!> Update local observation statistics
!!
  SUBROUTINE PDAF_incr_local_obsstats(dim_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l   ! Number of locally assimilated observations

! *** Update observation statistics ***

!$OMP CRITICAL
    IF (dim_obs_l > obsstats(4)) obsstats(4) = dim_obs_l
    IF (dim_obs_l > 0) THEN
       obsstats(3) = obsstats(3) + dim_obs_l
       obsstats(1) = obsstats(1) + 1
    ELSE
       obsstats(2) = obsstats(2) + 1
    END IF
!$OMP END CRITICAL

  END SUBROUTINE PDAF_incr_local_obsstats



!-------------------------------------------------------------------------------
!> Print observation statistics
!!
  SUBROUTINE PDAF_print_local_obsstats(screen, n_domains_with_obs)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, COMM_filter, MPIerr
    USE PDAFomi, &
         ONLY: omi_n_obstypes => n_obstypes, PDAFomi_obsstats_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: screen      ! Verbosity flag
    INTEGER, OPTIONAL, INTENT(out) :: n_domains_with_obs

! *** Local variables ***
    INTEGER :: obsstats_g(4)           ! Global statistics

    ! *** Print statistics for local analysis to the screen ***
    IF ( npes_filter>1) THEN
       CALL MPI_Reduce(obsstats, obsstats_g, 3, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(obsstats(4), obsstats_g(4), 1, MPI_INTEGER, MPI_MAX, &
            0, COMM_filter, MPIerr)
    ELSE
       ! This is a work around for working with nullmpi.F90
       obsstats_g = obsstats
    END IF

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- Global statistics for local analysis:'
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Local domains with observations:       ', obsstats_g(1)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Local domains without observations:    ', obsstats_g(2)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAF', 'Maximum local observation dimension:   ', obsstats_g(4)
       IF (obsstats_g(2) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Total avg. local observation dimension:', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1) + obsstats_g(2))
       END IF
       IF (obsstats_g(2) > 0 .AND. obsstats_g(1) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAF', 'Avg. for domains with observations:    ', &
               REAL(obsstats_g(3)) / REAL(obsstats_g(1))
       END IF
    END IF

    IF (omi_n_obstypes > 0) CALL PDAFomi_obsstats_l(screen)

    IF (PRESENT(n_domains_with_obs)) n_domains_with_obs = obsstats_g(1)

  END SUBROUTINE PDAF_print_local_obsstats


!-------------------------------------------------------------------------------
!> Operate SEIK matrix T on A as AT
!!
!! Operate matrix T on another matrix as
!!         $A = A T$.
!!
!! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
!! There are two proposed forms of T (ensemble size N):\\
!! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!!          last row= -1/N\\
!! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
!! We typically use TypeT=0, but both variants are implemented.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2002-01 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - moved into PDAF_analysis_utils in general revision
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF_seik_matrixT(dim, dim_ens, A)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim               !< dimension of states
    INTEGER, INTENT(in) :: dim_ens           !< Size of ensemble
    REAL, INTENT(inout) :: A(dim, dim_ens)   !< Input/output matrix
  
! *** local variables ***
    INTEGER :: row, col             ! counters
    INTEGER :: typeT = 0            ! Which type of T
    REAL :: invdimens               ! Inverse of ensemble size
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: rowmean(:) ! Mean values of rows of A

!$OMP THREADPRIVATE(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(rowmean(dim))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim)
       allocflag = 1
    END IF
    rowmean   = 0.0
    invdimens = 1.0 / REAL(dim_ens)

    IF (typeT == 0) THEN

       ! *** Compute row means of A ***
       DO col = 1, dim_ens
          DO row = 1, dim
             rowmean(row) = rowmean(row) + A(row, col)
          END DO
       END DO
       rowmean = invdimens * rowmean

    ELSE

       ! *** Get last column of A ***
       DO row = 1, dim
          rowmean(row) = A(row, dim_ens)
       END DO
     
    END IF


! **********************************************
! ***  Operate T on A                        ***
! ***                                        ***
! *** v^TT = (v_1-mean(v), ... ,v_r-mean(v)) ***
! *** with v = (v_1,v_2, ... ,r_(r+1))       ***
! **********************************************

    DO col = 1, dim_ens - 1
       DO row = 1, dim
          A(row, col) = A(row, col) - rowmean(row)
       END DO
    END DO
  
    DO row = 1, dim
       A(row, dim_ens) = 0.0
    END DO


! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(rowmean)

  END SUBROUTINE PDAF_seik_matrixT


!-------------------------------------------------------------------------------
!> Operate SEIK matrix T from left side on some matrix
!!
!! Operate matrix T on another matrix as
!!                 B = T A\\
!! \\
!! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
!! There are two proposed forms of T (ensemble size N):\\
!! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!!          last row= -1/N\\
!! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2002-01 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - moved into PDAF_analysis_utils in general revision
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF_seik_TtimesA(rank, dim_col, A, B)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: rank               !< Rank of initial covariance matrix
    INTEGER, INTENT(in) :: dim_col            !< Number of columns in A and B
    REAL, INTENT(in)    :: A(rank, dim_col)   !< Input matrix
    REAL, INTENT(out)   :: B(rank+1, dim_col) !< Output matrix (TA)

! *** local variables ***
    INTEGER :: row, col  ! counters
    INTEGER :: typeT = 0 ! Which type of T
    REAL :: invdimens    ! Inversize of ensemble size
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: colmean(:) ! Mean values of columns of A

!$OMP THREADPRIVATE(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(colmean(dim_col))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim_col)
       allocflag = 1
    END IF
    colmean = 0.0
    invdimens = -1.0 / REAL(rank + 1)

    whichT: IF (typeT == 0) THEN

       ! *** Compute column means of A ***
       DO col = 1, dim_col
          DO row = 1, rank
             colmean(col) = colmean(col) + invdimens * A(row, col)
          END DO
       END DO

    END IF whichT


! ****************************************************
! ***  Operate T on A                              ***
! ***                                              ***
! *** Tv_1 = (v_11-mean(v_1), ... ,v_r1-mean(v_1)) ***
! *** with v_1 = (v_11,v_21, ... ,v_N1  )          ***
! ****************************************************

    ! first DIM rows
    DO col = 1, dim_col
       DO row = 1, rank
          B(row, col) = A(row, col) + colmean(col)
       END DO
    END DO

    ! row RANK+1
    DO col = 1, dim_col
       B(rank + 1, col) = colmean(col)
    END DO


! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(colmean)

  END SUBROUTINE PDAF_seik_TtimesA


!------------------------------------------------------------------------------
!> Generate random matrix with special properties
!!
!! Generate a transformation matrix OMEGA for
!! the generation and transformation of the 
!! ensemble in the SEIK and LSEIK filter.
!! Generated is a uniform orthogonal matrix OMEGA
!! with R columns orthonormal in $R^{r+1}$
!! and orthogonal to (1,...,1)' by iteratively 
!! applying the Householder matrix onto random 
!! vectors distributed uniformly on the unit sphere.
!!
!! This version initializes at each iteration step
!! the whole Householder matrix and subsequently
!! computes Omega using GEMM from BLAS. All fields are 
!! allocated once at their maximum required size.
!! (On SGI O2K this is about a factor of 2.5 faster
!! than the version applying BLAS DDOT, but requires
!! more memory.)
!!
!! For Omegatype=0 a deterministic Omega is computed
!! where the Housholder matrix of (1,...,1)' is operated
!! on an identity matrix.
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2002-01 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_seik_Omega(rank, Omega, Omegatype, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_mod_filtermpi, &
         ONLY: mype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: rank                !< Approximated rank of covar matrix
    REAL, INTENT(inout) :: Omega(rank+1, rank) !< Matrix Omega
    INTEGER, INTENT(in) :: Omegatype           !< Select type of Omega:
                                               !<   (1) generated from random vectors
                                               !<   (0) generated from deterministic vectors
                                               !< (other) product of matrix from (2) with
                                               !<      orthonormal random matrix orthogonal (1....1)T
    INTEGER, INTENT(in) :: screen              !< Verbosity flag

!  *** local variables ***
    INTEGER :: col, row                  ! counters
    REAL :: rndval                       ! temporary value for init of Householder matrix
    REAL :: rndnum                       ! Value of randum entry
    REAL, ALLOCATABLE :: house(:,:)      ! Householder matrix
    REAL, POINTER :: rndmat(:,:)         ! Pointer to temporary Omega field


    randOmega: IF (Omegatype == 0) THEN 
! *************************************************
! *** Generate deterministic Omega as           ***
! *** Householder matrix associated with the    ***
! *** vector  1/sqrt(rank) (1,...,1)^T          ***
! *************************************************

       IF (mype == 0 .AND. screen > 0) &
            WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute deterministic Omega'

       rndnum = 1.0 / SQRT(REAL(rank + 1))

       ! First r rows
       rndval = - rndnum * rndnum / (rndnum + 1.0)
       Omegacolb: DO col = 1, rank
          Omegarowb: DO row = 1, rank
             Omega(row, col) = rndval
          END DO Omegarowb
       END DO Omegacolb

       DO col = 1, rank
          Omega(col, col) = Omega(col, col) + 1.0
       END DO

       ! Last row
       rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
       Omegacolc: DO col = 1, rank
          Omega(rank + 1, col) = rndval
       END DO Omegacolc

    ELSEIF (Omegatype == 1) THEN randOmega
! ****************************************
! *** Generate Omega by random vectors ***
! ****************************************

       IF (mype == 0 .AND. screen > 0) &
            WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute random Omega'

! *** Initialization ***

       ! Allocate fields
       ALLOCATE(house(rank + 1, rank))
       ALLOCATE(rndmat(rank, rank))

! *** Initialize orthonormal random matrix of size rank*rank ***

       CALL PDAF_generate_rndmat(rank, rndmat, 1)

! *** Project rndmat orthogonal to (1,...,1)^T ***

       ! *** Compute Householder matrix ***

       rndnum = 1.0 / SQRT(REAL(rank + 1))

       ! First r rows
       rndval = - rndnum * rndnum / (rndnum + 1.0)
       housecol: DO col = 1, rank
          houserow: DO row = 1, rank
             house(row, col) = rndval
          END DO houserow
       END DO housecol

       DO col = 1, rank
          house(col, col) = house(col, col) + 1.0
       END DO

       ! Last row
       rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
       housecolb: DO col = 1, rank
          house(rank + 1, col) = rndval
       END DO housecolb

       ! *** Complete Omega: house * rndmat ***

       CALL gemmTYPE ('n', 'n', rank + 1, rank, rank, &
            1.0, house, rank + 1, rndmat, rank, &
            0.0, Omega, rank + 1)

! *** CLEAN UP ***

       DEALLOCATE(house)
       DEALLOCATE(rndmat)

    ELSE randOmega
! *** Generate Omega as a product of a deterministic  ***
! *** transformation with an orthonormal random       ***
! *** matrix that preserves the mean.                 ***
! *** 1. The deterministic matrix matrix given by the ***
! *** householder matrix from Omegatype=0.            ***
! *** 2. The random matrix is generated analogously   ***
! *** to Omegatype=1 followed by a transformation to  ***
! *** ensure the (1,....,1)^T is an eigenvector of    ***
! *** the matrix.                                     ***

       IF (mype == 0 .AND. screen > 0) &
            WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute product Omega'

! *** Initialization ***

       ! Allocate fields
       ALLOCATE(house(rank + 1, rank))
       ALLOCATE(rndmat(rank, rank))

! *** 1. Deterministic part:                            ***
! *** Compute Householder matrix associated with the    ***
! *** vector  1/sqrt(rank) (1,...,1)^T                  ***
! *** (this is the transformation used for Omegatype=0) ***

       rndnum = 1.0 / SQRT(REAL(rank + 1))

       ! First r rows
       rndval = - rndnum * rndnum / (rndnum + 1.0)
       housecolc: DO col = 1, rank
          houserowc: DO row = 1, rank
             house(row, col) = rndval
          END DO houserowc
       END DO housecolc

       DO col = 1, rank
          house(col, col) = house(col, col) + 1.0
       END DO

       ! Last row
       rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
       housecalc: DO col = 1, rank
          house(rank + 1, col) = rndval
       END DO housecalc

! *** 2. Random part: 
! *** Initialize orthonormal random matrix of size rank*rank 
! *** with eigenvector (1,...,1)^T

       CALL PDAF_generate_rndmat(rank, rndmat, 2)

! *** 3. Multiply deterministic and random parts: 

       CALL gemmTYPE ('n', 'n', rank+1, rank, rank, &
            1.0, house, rank+1, rndmat, rank, &
            0.0, Omega, rank+1)

! *** CLEAN UP ***

       DEALLOCATE(house, rndmat)

    END IF randOmega

  END SUBROUTINE PDAF_seik_Omega


!------------------------------------------------------------------------------
!> Initialize matrix Uinv from matrix T
!!
!! Initialize matrix Uinv by
!! $U^{-1} = FAC\ T^T T$
!! where $FAC$ = rank+1 for a covariance matrix with factor
!! (rank+1)$^{-1}$ and $FAC$ = rank for a covariance matrix
!! with factor rank$^{-1}$.
!!
!! There are two proposed forms of T (ensemble size N):\\
!! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!!          last row= -1/N\\
!! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
!! We typically use TypeT=0, but both variants are implemented.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2002-01 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF_seik_Uinv(rank, Uinv)

    USE PDAF_seik, &
         ONLY: Nm1vsN

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: rank             !< Rank of initial covariance matrix
    REAL, INTENT(inout) :: Uinv(rank, rank) !< Inverse of matrix U
  
! *** local variables ***
    INTEGER :: row, col       ! counters
    INTEGER :: typeT = 0      ! Choose type of T
    REAL :: rdivrp1, r2divrp1 ! scaling factors for Uinv


! ***********************
! *** Initialize Uinv ***
! ***********************

    ttype: IF (typeT == 0) THEN

       ! Scaling factors
       IF (Nm1vsN == 1) THEN
          ! For ensemble covariance matrix with factor (N-1)^-1
          rdivrp1 = REAL(rank) / REAL(rank + 1)
          r2divrp1 = REAL(rank)**2 / REAL(rank + 1)
       ELSE
          ! For ensemble covariance matrix with factor N^-1
          rdivrp1 = 1
          r2divrp1 = REAL(rank)
       END IF

       DO col = 1, rank
          ! non-diagonal elements - upper triangle
          DO row = 1, col - 1
             Uinv(row, col) = - rdivrp1
          END DO
          ! diagonal
          Uinv(col, col) = r2divrp1
          ! non-diagonal elements - lower triangle
          DO row = col + 1, rank
             Uinv(row, col) = -rdivrp1
          END DO
       END DO

    ELSE ttype

       ! Scaling factors
       IF (Nm1vsN == 1) THEN
          ! For ensemble covariance matrix with factor (N-1)^-1
          rdivrp1 = REAL(rank) / REAL(rank + 1)
       ELSE
          ! For ensemble covariance matrix with factor N^-1
          rdivrp1 = 1
       END IF

       DO col = 1, rank
          ! non-diagonal elements - upper triangle
          DO row = 1, col - 1
             Uinv(row, col) = rdivrp1
          END DO
          ! diagonal
          Uinv(col, col) = 2.0 * rdivrp1
          ! non-diagonal elements - lower triangle
          DO row = col + 1, rank
             Uinv(row, col) = rdivrp1
          END DO
       END DO

    END IF ttype

  END SUBROUTINE PDAF_seik_Uinv

!------------------------------------------------------------------------------
!> Generate random matrix with special properties
!!
!! Generate a random matrix OMEGA for the initilization
!! of ensembles for the EnKF.
!!
!! The following properties can be set:\\
!! 1. Simply fill the matrix with random numbers from a 
!!   Gaussian distribution with mean zero and unit 
!!   variance. (This corresponds to the simple random 
!!   initialization of EnKF.)\\
!! 2. Constrain the columns of OMEGA to be of unit norm
!!   (This corrects error in the variance estimates caused
!!   by the finite number of random numbers.)\\
!! 3. Constrain the columns of OMEGA to be of norm dim\_ens$^{-1/2}$
!!   (This corrects variance errors as in 2)\\
!! 4. Project columns of OMEGA to be orthogonal to the vector
!!   $(1,....,1)^T$ by Householder reflections. (This assures
!!   that the mean of the generated ensemble equals the
!!   prescribed state estimate.)\\
!! Property 4 can be combined with either property 2 or 3.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2005-04 - Lars Nerger - Initial code PDAF_enkf_omega for EnKF
!! * 2025-02 - Lars Nerger - REname to general name and include in PDAF_analysis utils
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_ens_Omega(seed, r, dim_ens, Omega, norm, &
       otype, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: seed(4)  !< Seed for random number generation
    INTEGER, INTENT(in) :: r        !< Approximated rank of covar matrix
    INTEGER, INTENT(in) :: dim_ens  !< Ensemble size
    REAL, INTENT(inout) :: Omega(dim_ens,r)  !< Random matrix
    REAL, INTENT(inout) :: norm     !< Norm for ensemble transformation
    INTEGER, INTENT(in) :: otype    !< Type of Omega:
                                    !< (1) Simple Gaussian random matrix
                                    !< (2) Columns of unit norm
                                    !< (3) Columns of norm dim_ens^(-1/2)
                                    !< (4) Projection orthogonal (1,..,1)^T
                                    !< (6) Combination of 2 and 4
                                    !< (7) Combination of 3 and 4
                                    !< (8) Rows of sum 0 and variance 1
    INTEGER, INTENT(in) :: screen   !< Control verbosity

!  *** local variables ***
    INTEGER :: i, j                      ! Counters
    INTEGER, SAVE :: iseed(4)            ! seed array for random number routine
    REAL, ALLOCATABLE :: rndvec(:)       ! Vector of random numbers
    REAL, ALLOCATABLE :: house(:,:)      ! Householder matrix
    REAL, ALLOCATABLE :: Omega_tmp(:,:)  ! Temporary OMEGA for projection
    REAL :: colnorm                      ! Norm of matrix column
    REAL :: rownorm                      ! Norm of matrix row


! **********************
! *** INITIALIZATION ***
! **********************

    IF (seed(1) >= 0) THEN
       ! Use given seed
       iseed=seed
    END IF

    ALLOCATE(rndvec(r))

! *** Initialize random matrix ***

    DO i = 1, dim_ens
       ! Fill row-wise to be consistent with old sampling formulation
       CALL larnvTYPE(3, seed, r, rndvec)
       Omega(i, :) = rndvec
    END DO
  
! *** Normalize columns of Omega ***
    normcols: IF (otype == 2 .OR. otype == 3 .OR. otype == 6 .OR. otype == 7) THEN
       IF ((screen > 0) .AND. (otype == 2 .OR. otype == 6)) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- EnKF_Omega: Normalize columns of random matrix'
       ELSE IF (screen > 0) THEN
          WRITE (*, '(a, 5x,a)') &
               'PDAF', '--- EnKF_Omega: Normalize columns of random matrix to dim_ens^(-1/2)'
       END IF

       DO j = 1, r
          ! Compute norm
          colnorm = 0.0
          DO i = 1, dim_ens
             colnorm = colnorm + Omega(i, j)**2
          END DO
          IF (otype == 3 .OR. otype == 7) THEN
             ! Set column norm to 1/sqrt(dim_ens)
             colnorm = colnorm / REAL(dim_ens)
          END IF
          colnorm = SQRT(colnorm)

          ! Perform normalization
          DO i = 1, dim_ens
             Omega(i, j) = Omega(i, j) / colnorm
          END DO
       END DO
    END IF normcols


! *** Project columns orthogonal to (1,1,...,1)^T ***
    doproject: IF (otype == 4 .OR. otype == 6 .OR. otype == 7) THEN
       IF (screen > 0) &
            WRITE (*, '(a, 5x, a)') &
            'PDAF', '--- EnKF_Omega: Project columns orthogonal to (1,...,1)^T'

       ALLOCATE(Omega_tmp(dim_ens, r))
       ALLOCATE(house(dim_ens, dim_ens))

       ! Store Omega
       Omega_tmp = Omega

       ! Initialize Householder matrix
       housecolb: DO j = 1, dim_ens
          houserowb: DO i = 1, dim_ens
             house(i, j) = -1.0 / REAL(dim_ens)
          END DO houserowb
       END DO housecolb
       DO j = 1, dim_ens
          house(j, j) = house(j, j) + 1.0
       END DO

       ! Perform reflection
       CALL gemmTYPE ('n', 'n', dim_ens, r, dim_ens, &
            1.0, house, dim_ens, Omega_tmp, dim_ens, &
            0.0, Omega, dim_ens)

       DEALLOCATE(Omega_tmp, house)

    END IF doproject

    rowzero: IF (otype == 8) THEN
       IF (screen > 0) &
            WRITE (*, '(a, 5x, a)') &
            'PDAF', '--- EnKF_Omega: Ensure that row sums are zero'

       DO i = 1, dim_ens

          rownorm = 0.0

          DO j = 1, r
             rownorm = rownorm + Omega(i,j)
          ENDDO
          rownorm = rownorm / REAL(r)

          DO j = 1, r
             Omega(i,j) = Omega(i,j) - rownorm
          ENDDO

       END DO

       IF (screen > 0) &
            WRITE (*, '(a, 5x, a)') &
            'PDAF', '--- EnKF_Omega: Ensure that variance in rows is one'

       DO i = 1, dim_ens

          rownorm = 0.0

          DO j = 1, r
             rownorm = rownorm + Omega(i,j)*Omega(i,j)
          ENDDO
          rownorm = rownorm / REAL(r-1)
          rownorm = SQRT(rownorm)

          DO j = 1, r
             Omega(i,j) = Omega(i,j) / rownorm
          ENDDO
        
       END DO
    END IF rowzero


! *** Initialize norm for ensemble transformation ***
    IF (otype == 1) THEN
       norm = 1.0
    ELSEIF (otype == 2) THEN
       norm = SQRT(REAL(dim_ens - 1))
    ELSEIF (otype == 3) THEN
       norm = 1.0
    ELSEIF (otype == 4) THEN
       norm = 1.0
    ELSEIF (otype == 6) THEN
       norm = SQRT(REAL(dim_ens - 1))
    ELSEIF (otype == 7) THEN
       norm = 1.0
    END IF


! ********************
! *** Finishing up ***
! ********************

    DEALLOCATE(rndvec)

  END SUBROUTINE PDAF_ens_Omega


!-------------------------------------------------------------------------------
!> Operate matrix Omega on some matrix
!!
!! Operate matrix Omega on another matrix as
!!                 B = Omega A\\
!! 
!! Omega is a dim_ens x (dim_ens-1) matrix with matximum
!! rank and zero column sums. It is computed by the 
!! Householder reflection associate with the vector
!! (N-1)^(-1) (1,...,1)^T.
!!
!! The values of Omega are
!!    1 - 1 / (N (1/sqrt(N) + 1)) for i=j
!!    - 1 / (N (1/sqrt(N) + 1)) for i/=j, i<N
!!    - 1 / sqrt(N) for i=N
!!
!! In this routine the product A Omega is implemented as
!! operations:
!! 1. Compute the column sums of A
!! 2. Normalize column sums by 1/(sqrt(N) + N)
!! 3. Subtract value of last row multiplied by 1/(1+sqrt(N))
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_OmegaA(rank, dim_col, A, B)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: rank               !< Rank of initial covariance matrix
    INTEGER, INTENT(in) :: dim_col            !< Number of columns in A and B
    REAL, INTENT(in)    :: A(rank, dim_col)   !< Input matrix
    REAL, INTENT(out)   :: B(rank+1, dim_col) !< Output matrix (TA)
  
! *** local variables ***
    INTEGER :: row, col             ! counters
    REAL :: normsum                 ! Normalization for row sum
    REAL :: normlast                ! Normalization for last column
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: colsums(:) ! Mean values of columns of A

!$OMP threadprivate(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(colsums(dim_col))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim_col)
       allocflag = 1
    END IF
    colsums = 0.0

    ! Initialize normalization values
    normsum = 1.0 / REAL(rank+1) / (1.0/SQRT(REAL(rank+1))+1.0)
    normlast = - 1.0 / SQRT(REAL(rank+1))


    ! *** Compute column sums of A ***
    DO col = 1, dim_col
       DO row = 1, rank
          colsums(col) = colsums(col) + A(row, col)
       END DO
    END DO

! ****************************************************
! ***  Operate Omega on A                          ***
! ****************************************************

    ! Initialize last row of B
    DO col = 1, dim_col
       B(rank+1, col) = colsums(col) * normlast
    END DO

    ! Scale by NORMSUM
    colsums = normsum * colsums

    ! first rank rows
    DO col = 1, dim_col
       DO row = 1, rank
          B(row, col) = A(row, col) - colsums(col)
       END DO
    END DO


! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(colsums)

  END SUBROUTINE PDAF_estkf_OmegaA


!-------------------------------------------------------------------------------
!> Operate matrix Omega on A as AOmega
!!
!! Operate matrix Omega on another matrix as
!!         $A = A Omega$.
!!
!! Omega is a dim_ens x (dim_ens-1) matrix with matximum
!! rank and zero column sums. It is computed by the 
!! Householder reflection associate with the vector
!! (N-1)^(-1) (1,...,1)^T.
!!
!! The values of Omega are
!!    1 - 1 / (N (1/sqrt(N) + 1)) for i=j
!!    - 1 / (N (1/sqrt(N) + 1)) for i/=j, i<N
!!    - 1 / sqrt(N) for i=N
!!
!! In this routine the product A Omega is implemented as
!! operations:
!! 1. Compute the row sums of A
!! 2. Normalize row sums by 1/(sqrt(N) + N)
!! 3. Subtract value of last column multiplied by 1/(1+sqrt(N))
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_AOmega(dim, dim_ens, A)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments
    INTEGER, INTENT(in) :: dim               !< dimension of states
    INTEGER, INTENT(in) :: dim_ens           !< Size of ensemble
    REAL, INTENT(inout) :: A(dim, dim_ens)   !< Input/output matrix

! *** local variables ***
    INTEGER :: row, col  ! Counters
    REAL :: normsum      ! Normalization for row sum
    REAL :: normlast     ! Normalization for last column
    REAL :: val          ! Temporary variable
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: rowsums(:) ! Row sums of A

!$OMP threadprivate(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(rowsums(dim))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim)
       allocflag = 1
    END IF
    rowsums   = 0.0

    ! Initialize normalization values
    normsum = 1.0 / REAL(dim_ens) / (1.0/SQRT(REAL(dim_ens))+1.0)
    normlast = 1.0 / (1.0 + SQRT(REAL(dim_ens)))


    ! *** Compute row sums of A ***
    DO col = 1, dim_ens
       DO row = 1, dim
          rowsums(row) = rowsums(row) + A(row, col)
       END DO
    END DO

    ! Scale by NORMSUM
    rowsums = normsum * rowsums

    ! Substract scale value for last column
    DO row = 1, dim
       val = A(row, dim_ens) * normlast
       rowsums(row) = rowsums(row) + val
    END DO


! **********************************************
! ***  Operate Omega on A                    ***
! **********************************************

    DO col = 1, dim_ens - 1
       DO row = 1, dim
          A(row, col) = A(row, col) - rowsums(row)
       END DO
    END DO

    DO row = 1, dim
       A(row, dim_ens) = 0.0
    END DO


! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(rowsums)

  END SUBROUTINE PDAF_estkf_AOmega

!-------------------------------------------------------------------------------
!> Subtract row-wise means from array, e.g. to generate ensemble perturbation matrix
!!
!! This routine subtracts the mean value of each row from
!! a matrix. A use case is to subtract the ensemble mean state
!! from an ensemble matrix.
!!
!! This is a copy of PDAF_etkf_Tright with a clearer name.
!!
!! ! This is a core routine of PDAF and
!!   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code of PDAF_etkf_Tright
!! * 2024-12 - Lars Nerger - Renaming to clearer name
!! * 2025-02 - Lars Nerger - moved into PDAF_analysis_utils in general revision
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_subtract_rowmean(dim, dim_ens, A)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim               ! dimension of states
    INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
    REAL, INTENT(inout) :: A(dim, dim_ens)   ! Input/output matrix
  
! *** Local variables ***
    INTEGER :: row, col             ! counters
    REAL :: invdimens               ! Inverse of ensemble size
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: rowmean(:) ! Mean values of rows of A


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(rowmean(dim))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim)
       allocflag = 1
    END IF
    rowmean   = 0.0
    invdimens = 1.0 / REAL(dim_ens)

    ! *** Compute row means of A ***
    DO col = 1, dim_ens
       DO row = 1, dim
          rowmean(row) = rowmean(row) + A(row, col)
       END DO
    END DO
    rowmean = invdimens * rowmean


! **********************************************
! ***  Operate T on A                        ***
! ***                                        ***
! *** v^TT = (v_1-mean(v), ... ,v_r-mean(v)) ***
! *** with v = (v_1,v_2, ... ,r_N)           ***
! **********************************************

    DO col = 1, dim_ens
       DO row = 1, dim
          A(row, col) = A(row, col) - rowmean(row)
       END DO
    END DO
  

! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(rowmean)

  END SUBROUTINE PDAF_subtract_rowmean

!-------------------------------------------------------------------------------
!> Subtract column-wise means from an array
!!
!! This routine subtracts the mean value of each row from
!! a matrix. Usually the array is the transpose of an ensemble
!! matrix. This is used e.g. in the ETKF using a formulation
!! with T-matrix. This is the multiplication from the left: B = T A
!! where T is a symmetric dim_ens x dim_ens matrix with zero column 
!! sums defined as:  
!!            diag(T)=1-1/dim_ens; nondiag(T)=-1/dim_ens
!!
!! This is a copy of PDAF_etkf_Tleft with a clearer name.
!!
!! ! This is a core routine of PDAF and
!!   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code of PDAF_etkf_Tleft
!! * 2024-11 - Lars Nerger - Renaming for clearer name
!! * 2025-02 - Lars Nerger - moved into PDAF_analysis_utils in general revision
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_subtract_colmean(dim_ens, dim, A)

    USE PDAF_memcounting, &
         ONLY: PDAF_memcount

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_ens           ! Rank of initial covariance matrix
    INTEGER, INTENT(in) :: dim               ! Number of columns in A and B
    REAL, INTENT(inout) :: A(dim_ens, dim)   ! Input/output matrix
  
! *** Local variables ***
    INTEGER :: row, col             ! counters
    REAL :: invdimens               ! Inverse of ensemble size
    INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
    REAL, ALLOCATABLE :: colmean(:) ! Mean values of columns of A


! **********************
! *** INITIALIZATION ***
! **********************

    ALLOCATE(colmean(dim))
    IF (allocflag == 0) THEN
       ! count allocated memory
       CALL PDAF_memcount(3, 'r', dim)
       allocflag = 1
    END IF
    colmean = 0.0
    invdimens = 1.0 / REAL(dim_ens)

! *** Compute column means of A ***
    DO col = 1, dim
       DO row = 1, dim_ens
          colmean(col) = colmean(col) + A(row, col)
       END DO
    END DO
    colmean = invdimens * colmean


! ****************************************************
! ***  Operate T on A                              ***
! ***                                              ***
! *** Tv_1 = (v_11-mean(v_1), ... ,v_r1-mean(v_1)) ***
! *** with v_1 = (v_11,v_21, ... ,v_N)             ***
! ****************************************************

    ! first DIM columns
    DO col = 1, dim
       DO row = 1, dim_ens
          A(row, col) = A(row, col) - colmean(col)
       END DO
    END DO


! ********************
! *** FINISHING UP ***
! ********************

    DEALLOCATE(colmean)

  END SUBROUTINE PDAF_subtract_colmean

!-------------------------------------------------------------------------------
!> Set adaptive forgetting factor
!!
!! Dynamically set the global forgetting factor.
!! This is a typical implementation that tries to ensure
!! statistical consistency by enforcing the condition\\
!! var\_resid = 1/forget var\_ens + var\_obs\\
!! where var\_res is the variance of the innovation residual,
!! var\_ens is the ensemble-estimated variance, and
!! var\_obs is the observation error variance.\\
!! This routine is used in SEIK. It can also be used in LSEIK. 
!! In this case a forgetting factor for the PE-local domain is 
!! computed. An alternative for LSEIK is PDAF\_set\_forget\_local, 
!! which computes a forgetting factor for each local analysis 
!! domain. The implementation used in both routines is 
!! experimental and not proven to improve the estimates.
!!
!! __Revision history:__
!! * 2006-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, mens_p, &
       mstate_p, obs_p, U_init_obsvar, forget_in, forget_out, &
       screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE mpi
    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_mod_filtermpi, &
         ONLY: mype, npes_filter, MPIerr, COMM_filter

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: step                      !< Current time step
    INTEGER, INTENT(in) :: dim_obs_p                 !< Dimension of observation vector
    INTEGER, INTENT(in) :: dim_ens                   !< Ensemble size
    REAL, INTENT(in) :: mens_p(dim_obs_p, dim_ens)   !< Observed PE-local ensemble
    REAL, INTENT(in) :: mstate_p(dim_obs_p)          !< Observed PE-local mean state
    REAL, INTENT(in) :: obs_p(dim_obs_p)             !< Observation vector
    REAL, INTENT(in) :: forget_in                    !< Prescribed forgetting factor
    REAL, INTENT(out) :: forget_out                  !< Adaptively estimated forgetting factor
    INTEGER, INTENT(in) :: localfilter               !< Whether filter is domain-local
    INTEGER, INTENT(in) :: screen                    !< Verbosity flag

! *** External subroutine ***
!  (PDAF-internal name, real name is defined in the call to PDAF)
    EXTERNAL :: U_init_obsvar                        !< Initialize mean obs. error variance
  
! *** local variables ***
    INTEGER :: i, j                            ! Counters
    INTEGER :: dim_obs                         ! Global observation dimension for non-local filters
                                               ! PE-local dimension for local filters
    INTEGER :: dim_obs_g                       ! Global observation dimension
    REAL :: var_ens_p, var_ens                 ! Variance of ensemble
    REAL :: var_resid_p, var_resid             ! Variance of residual
    REAL :: var_obs                            ! Variance of observation errors
    REAL :: forget_neg, forget_max, forget_min ! Limiting values of forgetting factor


! **********************
! *** INITIALIZATION ***
! **********************

    ! Define limiting values of forgetting factor
    ! These are set very arbitrarily for now
    forget_neg = forget_in
    forget_max = 100.0
    forget_min = 0.01

    IF (mype == 0) THEN
       WRITE (*, '(a, 5x, a)') &
            'PDAF', '--- Apply global adaptive forgetting factor'
       WRITE (*, '(a, 9x, a, es10.2)') &
            'PDAF', 'Maximum limit for forgetting factor', forget_max
       WRITE (*, '(a, 9x, a, es10.2)') &
            'PDAF', 'Minimum limit for forgetting factor', forget_min
       WRITE (*, '(a, 9x, a, es10.2)') &
            'PDAF', 'Forgetting factor if var(obs) > var(resid)', forget_neg
    ENDIF


! ******************************************
! *** Compute adaptive forgetting factor ***
! ******************************************

    ! *** Compute mean ensemble variance ***

    CALL PDAF_timeit(51, 'new')

    IF (npes_filter>1) THEN
       CALL MPI_allreduce(dim_obs_p, dim_obs_g, 1, &
            MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
    ELSE
       dim_obs_g = dim_obs_p
    END IF

    IF (localfilter==0) THEN
       ! global fitlers 
       dim_obs = dim_obs_g
    ELSE
       ! domain-local filters
       dim_obs = dim_obs_p
    ENDIF

    IF (dim_obs_g > 0) THEN

       ! local
       var_ens_p = 0.0
       IF (dim_obs > 0) THEN
          DO i = 1, dim_obs_p
             DO j = 1, dim_ens
                var_ens_p = var_ens_p + (mstate_p(i) - mens_p(i, j)) ** 2
             ENDDO
          ENDDO
          var_ens_p = var_ens_p / REAL(dim_ens - 1) / REAL(dim_obs)
       END IF

       IF (localfilter==0) THEN
          ! global 
          CALL MPI_allreduce(var_ens_p, var_ens, 1, MPI_REALTYPE, MPI_SUM, &
               COMM_filter, MPIerr)
       ELSE
          ! For domain-local filters use only PE-local variance
          var_ens = var_ens_p
       ENDIF

       ! *** Compute mean of innovation ***
   
       ! Compute variance
       var_resid_p = 0.0
       IF (dim_obs > 0) THEN
          DO i = 1, dim_obs_p
             var_resid_p = var_resid_p + (obs_p(i) - mstate_p(i)) ** 2
          ENDDO
          var_resid_p = var_resid_p / REAL(dim_obs)
       END IF

       IF (localfilter==0) THEN
          ! global 
          CALL MPI_allreduce(var_resid_p, var_resid, 1, &
               MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
       ELSE
          ! For domain-local filters use only PE-local variance
          var_resid = var_resid_p
       ENDIF

       CALL PDAF_timeit(51, 'old')

       ! *** Compute mean observation variance ***

       ! Get mean observation error variance
       CALL PDAF_timeit(49, 'new')
       IF (dim_obs_p>0) THEN
          CALL U_init_obsvar(step, dim_obs_p, obs_p, var_obs)
       ELSE
          var_obs=0.0
       END IF
       CALL PDAF_timeit(49, 'old')

       CALL PDAF_timeit(51, 'new')

       ! *** Compute optimal forgetting factor ***
       IF (var_resid>0.0 .AND. var_obs>0.0) THEN
          forget_out = var_ens / (var_resid - var_obs)
       ELSE
          forget_out = forget_in
       END IF

       ! Apply special condition if observation variance is larger than residual variance
       IF (forget_out < 0.0) forget_out = forget_neg

       ! Impose upper limit for forgetting factor
       IF (forget_out > forget_max) forget_out = forget_max

       ! Impose lower limit for forgetting factor
       IF (forget_out < forget_min) forget_out = forget_min

    ELSE
       ! No observations available in full model domain
       forget_out = forget_in
    END IF


! ********************
! *** FINISHING UP ***
! ********************

    IF (mype == 0 .AND. screen>0) THEN
       WRITE (*, '(a, 9x, a, es10.2)') &
            'PDAF', '--> Computed forgetting factor', forget_out
    ENDIF

    CALL PDAF_timeit(51, 'old')
   
  END SUBROUTINE PDAF_set_forget

!-------------------------------------------------------------------------------
!> Set local adaptive forgetting factor
!!
!! Dynamically set the global forgetting factor individually for
!! each local analysis domain in LSEIK. 
!! This is a typical implementation that tries to ensure
!! statistical consistency by enforcing the condition\\
!! var\_resid = 1/forget var\_ens + var\_obs\\
!! where var\_res is the variance of the innovation residual,
!! var\_ens is the ensemble-estimated variance, and
!! var\_obs is the observation error variance.\\
!! This variant is not proven to improve the estimates!
!!
!! __Revision history:__
!! * 2006-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_set_forget_local(domain, step, dim_obs_l, dim_ens, &
       HX_l, HXbar_l, obs_l, U_init_obsvar_l, forget, aforget)

    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_mod_filtermpi, &
         ONLY: mype
#if defined (_OPENMP)
    USE omp_lib, &
         ONLY: omp_get_num_threads, omp_get_thread_num
#endif

    IMPLICIT NONE

! *** Arguments
    INTEGER, INTENT(in) :: domain                !< Current local analysis domain
    INTEGER, INTENT(in) :: step                  !< Current time step
    INTEGER, INTENT(in) :: dim_obs_l             !< Dimension of local observation vector
    INTEGER, INTENT(in) :: dim_ens               !< Ensemble size
    REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens) !< Local observed ensemble
    REAL, INTENT(in) :: HXbar_l(dim_obs_l)       !< Local observed state estimate
    REAL, INTENT(in) :: obs_l(dim_obs_l)         !< Local observation vector
    REAL, INTENT(in) :: forget                   !< Prescribed forgetting factor
    REAL, INTENT(out) :: aforget                 !< Adaptive forgetting factor

! *** External subroutines ***
! ! (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_obsvar_l                  !< Initialize local mean obs. error variance
  
! *** local variables ***
    INTEGER :: i, j                     ! Counters
    REAL :: var_ens, var_resid, var_obs ! Variances
    INTEGER, SAVE :: first = 1          ! Flag for very first call to routine
    INTEGER, SAVE :: lastdomain = -1    ! store domain index
    LOGICAL, SAVE :: screenout = .true. ! Whether to print information to stdout
    INTEGER, SAVE :: mythread, nthreads ! Thread variables for OpenMP
    REAL :: forget_neg, forget_max, forget_min ! limiting values of forgetting factor


! **********************
! *** INITIALIZATION ***
! **********************

    ! Define limiting values of forgetting factor
    ! These are set very arbitrarily for now
    forget_neg = forget
    forget_max = 100.0
    forget_min = 0.01

#if defined (_OPENMP)
    nthreads = omp_get_num_threads()
    mythread = omp_get_thread_num()
#else
    nthreads = 1
    mythread = 0
#endif


    ! Control screen output
    IF (lastdomain<domain .AND. lastdomain>-1) THEN
       screenout = .false.
    ELSE
       screenout = .true.

       ! In case of OpenMP, let only thread 0 write output to the screen
       IF (mythread>0) screenout = .false.

       ! Output, only in case of OpenMP parallelization
    END IF


! ****************************************************
! *** Initialize adaptive local forgetting factors ***
! ****************************************************

    IF (screenout) THEN
       ! At first call during each forecast phase
       IF (mype == 0) THEN
          WRITE (*, '(a, 5x, a)') &
               'PDAF', '--- Apply dynamically estimated local forgetting factors'
          WRITE (*, '(a, 9x, a, es10.2)') &
               'PDAF', 'Maximum limit for forgetting factor', forget_max
          WRITE (*, '(a, 9x, a, es10.2)') &
               'PDAF', 'Minimum limit for forgetting factor', forget_min
          WRITE (*, '(a, 9x, a, es10.2)') &
               'PDAF', 'Forgetting factor if var(obs)>var(resid)', forget_neg
       ENDIF

       ! Set flag
       first = 0
    ENDIF
    lastdomain = domain


    ! ************************************************************
    ! *** Compute optimal forgetting factor for current domain ***
    ! ************************************************************

    ! *** Compute mean ensemble variance ***

    ! local
    var_ens = 0.0
    DO i = 1, dim_obs_l
       DO j = 1, dim_ens
          var_ens = var_ens + (HXbar_l(i) - HX_l(i, j)) ** 2
       ENDDO
    ENDDO
    var_ens = var_ens / REAL(dim_ens - 1) / REAL(dim_obs_l)


    ! *** Compute mean of observation-minus-forecast residual ***
   
    var_resid = 0.0
    DO i = 1, dim_obs_l
       var_resid = var_resid + (obs_l(i) - HXbar_l(i)) ** 2
    ENDDO
    var_resid = var_resid / REAL(dim_obs_l)


    ! *** Compute mean observation variance ***

    ! Get mean observation error variance
    CALL PDAF_timeit(52, 'new')
    CALL U_init_obsvar_l(domain, step, dim_obs_l, obs_l, var_obs)
    CALL PDAF_timeit(52, 'old')

    ! *** Compute optimal forgetting factor ***
    aforget = var_ens / (var_resid - var_obs)

    ! Apply special condition if observation variance is larger than residual variance
    IF (aforget < 0.0) aforget = forget_neg

    ! Impose upper limit for forgetting factor
    ! - the value for this is quite arbitary
    IF (aforget > forget_max) aforget = forget_max

    ! Impose lower limit for forgetting factor
    ! - the value for this is quite arbitary
    IF (aforget < forget_min) aforget = forget_min

  END SUBROUTINE PDAF_set_forget_local

!-------------------------------------------------------------------------------
!> Add noise to particles after resampling
!!
!! Adding noise to particles to avoid identical particles
!! after resampling.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-07 - Lars Nerger initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p, type_noise, noise_amp, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &
         ONLY: mype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< State dimension
    INTEGER, INTENT(in) :: dim_ens               !< Number of particles
    REAL, INTENT(inout) :: state_p(dim_p)        !< State vector (not filled)
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< Ensemble array
    INTEGER, INTENT(in) :: type_noise            !< Type of noise
    REAL, INTENT(in)    :: noise_amp             !< Noise amplitude
    INTEGER, INTENT(in) :: screen                !< Verbosity flag

! *** local variables ***
    INTEGER :: i, member                ! Loop counters
    INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
    INTEGER, SAVE :: first = 1          ! flag for init of random number seed
    INTEGER, SAVE :: iseed(4)           ! seed array for random number routine
    REAL, ALLOCATABLE :: ens_noise(:,:) ! Noise to be added for PF
    REAL :: noisenorm                   ! output argument of PDAF_ens_Omega (not used)
    REAL :: invdim_ens                  ! Inverse ensemble size
    REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
    REAL :: variance                    ! Ensmeble variance


! **********************
! *** INITIALIZATION ***
! **********************

    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 5x, a)') &
            'PDAF', 'Perturb particles:'
       IF (type_noise == 1) THEN
          WRITE (*, '(a, 5x, a, f10.4)') &
               'PDAF', '--- Gaussian noise with constant standard deviation', noise_amp
       ELSEIF (type_noise == 2) THEN
          WRITE (*, '(a, 5x, a, es10.3, a)') &
               'PDAF', '--- Gaussian noise with amplitude ', noise_amp,' of ensemble standard deviation'
       END IF
    END IF

    ! Initialized seed for random number routine
    IF (first == 1) THEN
       iseed(1) = 1000
       iseed(2) = 2045+mype
       iseed(3) = 10
       iseed(4) = 3
       first = 2
    END IF

    ! Initialize numbers
    invdim_ens    = 1.0 / REAL(dim_ens)  
    invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! ******************************
! *** Add noise to particles ***
! ******************************

    ALLOCATE(ens_noise(1, dim_ens))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)


    IF (type_noise == 1) THEN

       ! *** Noise with constant standard deviation ***

       DO i = 1, dim_p
          CALL PDAF_ens_Omega(iseed, dim_ens, 1, ens_noise(1,:), noisenorm, 8, 0)

          DO member = 1, dim_ens
             ens_p(i, member) = ens_p(i, member) + noise_amp * ens_noise(1,member)
          END DO
       END DO

    ELSEIF (type_noise == 2) THEN

       ! *** Noise with fraction of ensemble standard deviation ***

       ! Compute mean state
       state_p = 0.0
       DO member = 1, dim_ens
          DO i = 1, dim_p
             state_p(i) = state_p(i) + ens_p(i, member)
          END DO
       END DO
       state_p(:) = invdim_ens * state_p(:)

       ! Add noise
       DO i = 1, dim_p

          ! Initialize noise vector with zero mean and unit variance
          CALL PDAF_ens_Omega(iseed, dim_ens, 1, ens_noise(1,:), noisenorm, 8, 0)

          ! Compute sampled variance
          variance = 0.0
          DO member = 1, dim_ens
             variance = variance + (ens_p(i, member) - state_p(i))*(ens_p(i, member) - state_p(i))
          END DO
          variance = invdim_ensm1 * variance

          ! Add noise to particles
          DO member = 1, dim_ens
             ens_p(i, member) = ens_p(i, member) + noise_amp * sqrt(variance) * ens_noise(1, member)
          END DO
       END DO

    END IF


! ********************
! *** Finishing up ***
! ********************

    DEALLOCATE(ens_noise)

    IF (allocflag == 0) allocflag = 1

  END SUBROUTINE PDAF_add_particle_noise

!-------------------------------------------------------------------------------
!> Inflation for particle filters
!!
!! This routine compute an adaptive inflation using the effective sample
!! size N_eff according to
!!      N_eff / N >= alpha
!! whether N_eff is itertively computed with increasing inflation of the
!! observation error variance.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-08 - Lars Nerger
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_inflate_weights(screen, dim_ens, alpha, weights)

    USE PDAF_diag, ONLY: PDAF_diag_effsample

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: screen            !< verbosity flag
    INTEGER, INTENT(in) :: dim_ens           !< Ensemble size
    REAL, INTENT(inout) :: weights(dim_ens)  !< weights (before and after inflation)
    REAL, INTENT(in) :: alpha                !< Minimum limit of n_eff / N

! *** Local variables ***
    INTEGER :: i
    REAL :: alpha_iter
    REAL, ALLOCATABLE :: logw(:)      ! Logarith of particle weights
    REAL, ALLOCATABLE :: aweights(:)  ! temporary weights
    REAL :: a_step, tot_weight
    REAL :: alpha_lim, n_eff


! ******************
! *** Initialize ***
! ******************

    ALLOCATE(logw(dim_ens))
    ALLOCATE(aweights(dim_ens))

    ! Get logarithm of weights
    DO i=1, dim_ens
       logw(i) = LOG(weights(i))
    END DO

    ! Store initial weigts
    aweights = weights

    IF (screen>0) THEN
       WRITE (*,'(a, 5x, a, F10.3)') &
            'PDAF','--- Inflate weights according to N_eff/N > ', alpha
    END IF
  

! **********************************************
! *** Determine inflation according to alpha ***
! **********************************************

    alpha_iter = 0.0
    a_step = 0.05

    ! Set limit
    alpha_lim = alpha * REAL(dim_ens)

    aloop: DO

       IF (alpha_iter >= 1.0) THEN
          alpha_iter = 1.0
          EXIT aloop
       END IF

       ! scale 
       DO i = 1, dim_ens
          aweights(i) = EXP(logw(i) * (1.0-alpha_iter))
       END DO
  
       ! Normalize weights
       tot_weight = 0.0
       DO i = 1, dim_ens
          tot_weight = tot_weight + aweights(i)
       END DO
       IF (tot_weight /= 0.0) THEN
          aweights = aweights / tot_weight

          ! Compute effective ensemble size
          CALL PDAF_diag_effsample(dim_ens, aweights, n_eff)

          ! If limiting condition is fullfileed, exit loop
          IF (REAL(n_eff) > alpha_lim) EXIT aloop
       END IF

       alpha_iter = alpha_iter + a_step

    END DO aloop

    ! Store final inflated weigts
    weights = aweights


! ********************
! *** Finishing up ***
! ********************

    DEALLOCATE(aweights, logw)

  END SUBROUTINE PDAF_inflate_weights

!-------------------------------------------------------------------------------
!> Inflate ensemble spread according for forgetting factor
!!
!! This routine modifies an input ensemble such that its covariance 
!! is inflated by the factor 1/forget.  The ensemble perturbations
!! are inflated by 1/sqrt(forget). The ensemble mean is unchanged:
!!      Xnew = Xmean + 1/sqrt(forget) * (X - Xmean) 
!! The routine returns the inflated ensemble, replacing the input ensemble
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-11 - Julian Toedter - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_inflate_ens(dim, dim_ens, meanstate, ens, forget, do_ensmean)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim               !< dimension of states
    INTEGER, INTENT(in) :: dim_ens           !< Size of ensemble
    REAL, INTENT(inout) :: meanstate(dim)    !< state vector to hold ensemble mean
    REAL, INTENT(inout) :: ens(dim, dim_ens) !< Input/output ensemble matrix
    REAL, INTENT(in)    :: forget            !< Forgetting factor
    LOGICAL, INTENT(in) :: do_ensmean        !< Whether to compute the ensemble mean state

! *** local variables ***
    INTEGER :: row, col  ! counters
    REAL :: invdimens    ! Inverse of ensemble size
    REAL :: infl         ! inflation factor of perturbations = 1/sqrt(forget)


! **********************
! *** INITIALIZATION ***
! **********************

    IF (do_ensmean) THEN
       invdimens = 1.0 / REAL(dim_ens)

       ! Compute ensemble mean state 
       meanstate   = 0.0
       DO col = 1, dim_ens
          DO row = 1, dim
             meanstate(row) = meanstate(row) + ens(row, col)
          END DO
       END DO
       meanstate = invdimens * meanstate
    END IF


! **********************************************
! ***  Subtract mean and apply inflation     ***
! **********************************************

    ! Get perturbation matrix X'=X-xmean
    DO col = 1, dim_ens
       DO row = 1, dim
          ens(row, col) = ens(row, col) - meanstate(row)
       END DO
    END DO

    infl = 1.0/SQRT(forget)

    ! Inflation is done by X'new = infl * X 
    ! Add mean again to get Xnew = X'new + ensmean
    DO col = 1, dim_ens
       DO row = 1, dim
          ens(row, col) = infl * ens(row, col) + meanstate(row)
       END DO
    END DO

  END SUBROUTINE PDAF_inflate_ens


!-------------------------------------------------------------------------------
!> Compute weight for localization
!!
!! This routine initializates a single weight based on the given
!! distance and localization radii for the specified weighting
!! type.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2010-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
     nrows, ncols, A, var_obs, weight, verbose)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: wtype      !< Type of weight function
          !< * (0): unit weight (=1 up to distance=cradius)
          !< * (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
          !< * (2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
  INTEGER, INTENT(in) :: rtype      !< Type of regulated weighting
          !< * (0): no regulation
          !< * (1): regulation using mean variance
          !< * (2): regulation using variance of single observation point
  REAL, INTENT(in)    :: cradius    !< Cut-off radius
  REAL, INTENT(in)    :: sradius    !< Support radius 
  REAL, INTENT(in)    :: distance   !< Distance to observation
  INTEGER, INTENT(in) :: nrows      !< Number of rows in matrix A
  INTEGER, INTENT(in) :: ncols      !< Number of columns in matrix A
  REAL, INTENT(in) :: A(nrows, ncols) !< Input matrix
  REAL, INTENT(in)    :: var_obs    !< Observation variance
  REAL, INTENT(out)   :: weight     !< Weights
  INTEGER, INTENT(in) :: verbose    !< Verbosity flag

! *** Local variables ***
  INTEGER :: i, j                   ! Counters
  REAL    :: cfaci                  ! parameter for initialization of 5th-order polynomial
  REAL    :: meanvar                ! Mean variance in observation domain
  REAL    :: svarpovar              ! Mean state plus observation variance
  REAL    :: var                    ! variance for Gaussian
  REAL, PARAMETER :: pi=3.141592653589793   !Pi


! ********************************
! *** Print screen information ***
! ********************************

  ptype: IF (verbose == 1) THEN
     WRITE (*,'(a, 5x, a)') 'PDAF','Set localization weights'
     IF (wtype == 0) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize unit weights'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Support radius ', sradius
        IF (cradius < sradius) THEN
           WRITE (*, '(a, 5x, a, es12.4)') &
                'PDAF', '--- Use cut-off radius ', cradius
        END IF
     ELSE IF (wtype == 1) THEN
       WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize exponential weight function'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Distance for 1/e   ', sradius
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Cut-off radius ', cradius
     ELSE IF (wtype == 2) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize weights by 5th-order polynomial'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Support radius ', sradius
        IF (cradius < sradius) THEN
          WRITE (*, '(a, 5x, a, es12.4)') &
                'PDAF', '--- Use cut-off radius ', cradius
        END IF
     ELSE IF (wtype == 3) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize weights by scaled Gaussian function'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Standard deviation ', sradius
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Cut-off radius ', cradius
     END IF
  END IF ptype


! **************************
! *** Initialize weights ***
! **************************

  t_weight: IF (wtype == 0) THEN
     ! Unit weights

     IF (distance <= cradius) THEN
        weight = 1.0
     ELSE
        weight = 0.0
     END IF

  ELSE IF (wtype == 1) THEN t_weight
     ! Weighting by exponential decrease

     IF (cradius > 0.0 .AND. sradius > 0.0) THEN

        IF (distance <= cradius) THEN
           weight = EXP(-distance / sradius)
        ELSE
           weight = 0.0
        END IF

     ELSE

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF

  ELSE IF (wtype == 2) THEN t_weight
     ! Weighting by the square-root of a 5th-order function;
     ! equation (4.10) of Gaspari&Cohn, QJRMS125, 723 (1999)

     cfaci = REAL(sradius) / 2.0

     ! Compute weight
     cradnull: IF (cradius > 0.0 .and. sradius > 0.0) THEN

        cutoff: IF (distance <= cradius) THEN
           IF (distance <= sradius / 2.0) THEN
              weight = -0.25 * (distance / cfaci)**5 &
                   + 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   - 5.0 / 3.0 * (distance / cfaci)**2 + 1.0
           ELSEIF (distance > sradius / 2.0 .AND. distance < sradius * 0.9) THEN
              weight = 1.0 / 12.0 * (distance / cfaci)**5 &
                   - 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   + 5.0 / 3.0 * (distance / cfaci)**2 &
                   - 5.0 * (distance / cfaci) &
                   + 4.0 - 2.0 / 3.0 * cfaci / distance
           ELSEIF (distance >= sradius * 0.9 .AND. distance < sradius) THEN
              ! Ensure that weight is non-negative
              weight = MAX(1.0 / 12.0 * (distance / cfaci)**5 &
                   - 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   + 5.0 / 3.0 * (distance / cfaci)**2 &
                   - 5.0 * (distance / cfaci) &
                   + 4.0 - 2.0 / 3.0 * cfaci / distance, 0.0)
           ELSE
              weight = 0.0
           ENDIF
        ELSE cutoff
           weight = 0.0
        END IF cutoff

     ELSE cradnull

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF cradnull

  ELSE IF (wtype == 3) THEN t_weight
     ! Weighting by Gaussian function scaled for w(0)=1.0

     ! Compute weight
     IF (cradius > 0.0 .and. sradius > 0.0) THEN

        IF (distance <= cradius) THEN
           var = sradius*sradius
           weight = exp(-distance*distance/ (2.0*var))
        ELSE
           weight = 0.0
        END IF

     ELSE

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF



  END IF t_weight


! *********************************
! *** Compute weight regulation ***
! *********************************

  regweight: IF (rtype == 1) THEN
     ! Regulated weight with a function based on the fixed weight
     ! function, the variance of the observations and the local mean of
     ! the estimated state error variance for the local analysis domain.

     IF (verbose == 1) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Compute regulated weight'
     END IF

     ! Compute mean variance from ensemble perturbations
     meanvar = 0.0
     DO i = 1, nrows
        DO j = 1, ncols
           meanvar = meanvar + A(i,j)**2
        END DO
     END DO
     meanvar = meanvar / REAL(ncols) / REAL(nrows)

     svarpovar = meanvar + var_obs

     ! Compute regulated weight
     weight = weight * var_obs / svarpovar &
          / (1.0 - (weight * meanvar / svarpovar))

  END IF regweight


END SUBROUTINE PDAF_local_weight

!-------------------------------------------------------------------------------
!> Compute weight functions for localization
!!
!! This routine initializates a vector holding weight coefficients
!! for localization. The weights can be applied to localization
!! on the state covariance matrix of the observation error 
!! covariance matrix.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2010-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_local_weights(wtype, cradius, sradius, dim, distance, &
     weight, verbose)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: wtype          !< Type of weight function
          !< * (0): unit weight (=1 up to distance=cradius)
          !< * (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
          !< * (2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
  REAL, INTENT(in)    :: cradius        !< Parameter for cut-off
  REAL, INTENT(in)    :: sradius        !< Support radius 
  INTEGER, INTENT(in) :: dim            !< Size of distance and weight arrays
  REAL, INTENT(in)    :: distance(dim)  !< Array holding distances
  REAL, INTENT(out)   :: weight(dim)    !< Array for weights
  INTEGER, INTENT(in) :: verbose        !< Verbosity flag
  

! *** Local variables ***
  INTEGER :: i          ! Counter
  INTEGER :: verbose_s  ! verbosity flag
  REAL    :: distance_s ! Distance for single observation
  REAL    :: weight_s   ! Weight for single observation
  REAL    :: matA(1,1)  ! Tempoarary array required in call to PDAF_local_weight


! *******************************
! *** Initialize weight array ***
! *******************************

  DO i = 1, dim

     IF (verbose >= 1 .AND. i == 1) THEN
        verbose_s = 1
     ELSE
        verbose_s = 0
     END IF

     ! Set distance
     distance_s = distance(i)

     ! Compute weight
     CALL PDAF_local_weight(wtype, 0, cradius, sradius, distance_s, &
     1, 1, matA, 0.0, weight_s, verbose_s)

     ! Initialize element of weight array
     weight(i) = weight_s

  END DO

END SUBROUTINE PDAF_local_weights

END MODULE PDAF_analysis_utils
