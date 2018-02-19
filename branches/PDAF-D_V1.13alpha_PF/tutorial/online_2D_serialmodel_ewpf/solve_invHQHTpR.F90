!$Id$
!BOP
!
! !ROUTINE: solve_invHQHTpR --- invert HQH^T + R
!
! !INTERFACE:
subroutine solve_invHQHTpR(step, dim_obs, dim_p, dim_ens, input, output, model_state)

! !DESCRIPTION:
! User-supplied routine for PDAF.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner
! Later revisions - see svn log
!
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! Arguments
  integer, intent(in) :: dim_obs, dim_p, step, dim_ens
  real , intent(in) :: input(dim_obs)
  real, intent(out) :: output(dim_obs)
  real, intent(in) :: model_state(dim_p)


! Local variables
  INTEGER :: i
  INTEGER :: info
  REAL :: cp_input(dim_obs)

  REAL, ALLOCATABLE :: cov_r(:)
  REAL, ALLOCATABLE :: Htcov(:)
  REAL, ALLOCATABLE :: QHtcov(:)
  REAL, ALLOCATABLE, SAVE :: HQHtcov(:,:) 
  INTEGER, SAVE :: cnt_ens=0
  REAL, ALLOCATABLE :: HQHtcov_tmp(:,:)
  
! Increment ensemble counter
  cnt_ens =  cnt_ens + 1

  ALLOCATE(cov_r(dim_obs))
  ALLOCATE(HTcov(dim_p))
  ALLOCATE(QHTcov(dim_p))

  ALLOCATE(HQHtcov_tmp(dim_obs,dim_obs))

  ! initialize Identity in observation space

  IF (cnt_ens == 1) THEN
     
     ALLOCATE( HQHtcov(dim_obs,dim_obs) )
     HQHtcov = 0
     cov_r = 0

     do i =1,dim_obs
        Htcov = 0
        cov_r(i) = 1.0

        CALL adj_obs_op_pdaf(step, dim_obs, dim_p, cov_r, Htcov)
        CALL prodQx(dim_p,Htcov ,model_state,QHtcov)
        CALL obs_op_pdaf(step, dim_p, dim_obs, QHtcov, HQHtcov(:,i))

        cov_r(i) = 0
     enddo

     do i = 1, dim_obs
        HQHtcov(i,i) = HQHtcov(i,i) + rms_obs**2
     enddo

  ENDIF
 
  HQHtcov_tmp = HQHtcov
  cp_input = input

  !Inversion: solve x = Qy
  call DPOSV( 'U', dim_obs, 1, HQHtcov_tmp, dim_obs, &
       cp_input, dim_obs, info)

  IF (info .ne. 0 ) then
     Write(*,*) 'There is an error in the inversion!!!! info =',info
  end IF

  output = cp_input

  IF (cnt_ens == dim_ens) THEN
     cnt_ens = 0
     HQHtcov = 0 
     DEALLOCATE(HQHTcov)
  ENDIF

! Clean up

  DEALLOCATE(HQHtcov_tmp)
  DEALLOCATE(cov_r, HTcov, QHTcov)

END subroutine solve_invHQHTpR

