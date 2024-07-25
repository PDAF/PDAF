module analytical_correlation_functions

  implicit none

    ! Referenes Gaspari & cohn:
    !
    ! https://link.springer.com/chapter/10.1007/978-3-662-47100-5_15#Sec8
    !
    ! https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1256/qj.05.08 (C1 und C2)
    !
    ! https://link.springer.com/content/pdf/10.1007/s00190-007-0195-4.pdf?pdf=inline%20link (41) and (42)
   type :: gaspari_cohn_1
    real, dimension(4,6), private :: B
    real, dimension(4), private :: d
    real, private :: a = huge(1.)
    real, public :: R = huge(1.)
    contains
    procedure, pass(this) :: init => gaspari_cohn_1_init
    procedure, pass(this), private :: gaspari_cohn_1_evaluate_scalar
    procedure, pass(this), private :: gaspari_cohn_1_evaluate_array
    procedure, pass(this) :: gaspari_cohn_1_evaluate_norm_scalar
    procedure, pass(this) :: gaspari_cohn_1_evaluate_norm_array
    generic :: evaluate => gaspari_cohn_1_evaluate_norm_scalar, gaspari_cohn_1_evaluate_norm_array
  end type

  contains

  subroutine gaspari_cohn_1_init(this,a,R)
    implicit none

    class(gaspari_cohn_1) :: this
    real, intent(in) :: a
    real, intent(in) :: R

    if(this%R/=R) then
      this%R=R
    end if
    if(this%a/=a) then
      this%a=a

      this%B(1,6)=-16.*(3.-8.*a + 7.*a**2)/3.
      this%B(1,5)=+16.*(1.-2.*a + 2.*a**2)
      this%B(1,4)=+10.*(1.-4.*a + 8.*a**2)
      this%B(1,3)=-40.*(1.-2.*a + 8.*a**2)/3.
      this%B(1,2)=0.
      this%B(1,1)=2. + 6.*a + 44.*a**2

      this%B(2,6)=+16.*(1.-6.*a + 5.*a**2)/3.
      this%B(2,5)=-8.*(2. - 10.*a + 8.*a**2)
      this%B(2,4)=+10.
      this%B(2,3)=+20.*(2.-22.*a + 20.*a**2)/3.
      this%B(2,2)=-5.*(4. - 26.*a + 36.*a**2)
      this%B(2,1)= 8.-35.*a + 102.*a**2

      this%B(3,6)=+16.*a*(2.-3.*a)/3.
      this%B(3,5)=-16.*a*(3.-4.*a)
      this%B(3,4)=+40.*a*(1.-a)
      this%B(3,3)=+40.*a*(9.-10.*a)/3.
      this%B(3,2)=-10.*a*(27.-22.*a)
      this%B(3,1)=a*(189.-122.*a)

      this%B(4,6)=+16.*a**2/3.
      this%B(4,5)=-32.*a**2
      this%B(4,4)=+40.*a**2
      this%B(4,3)=+320.*a**2/3
      this%B(4,2)=-320.*a**2
      this%B(4,1)=+256.*a**2

      this%d(1) = 0.
      this%d(2) = (-4. + 29.*a - 42.*a**2)/6.
      this%d(3) = -a*(243. - 230.*a)/6.
      this%d(4) = -128.*a**2/3.
    end if

  end subroutine

  function gaspari_cohn_1_evaluate_norm_array(this,tau) result(c)

    implicit none

    ! arguments
    class(gaspari_cohn_1) :: this
    real, dimension(:), intent(in) :: tau

    ! result
    real, dimension(size(tau)) :: c

    ! local
    real :: c0

    call this%gaspari_cohn_1_evaluate_array(tau,c)
    call this%gaspari_cohn_1_evaluate_scalar(0.0,c0)

    c = c/c0

  end function

  function gaspari_cohn_1_evaluate_norm_scalar(this,tau) result(c)

    implicit none

    ! arguments
    class(gaspari_cohn_1) :: this
    real, intent(in) :: tau

    ! result
    real :: c

    ! local
    real :: c0

    call this%gaspari_cohn_1_evaluate_scalar(tau,c)
    call this%gaspari_cohn_1_evaluate_scalar(0.0,c0)

    c = c/c0

  end function

  subroutine gaspari_cohn_1_evaluate_array(this,tau,c)

    implicit none

    ! arguments
    class(gaspari_cohn_1) :: this
    real, dimension(:), intent(in) :: tau
    real, dimension(:), intent(out) :: c

    ! local
    real, dimension(:,:), allocatable :: F
    integer :: i,j
    integer, dimension(:), allocatable :: case_idx

    allocate(F(6,size(tau)))

    F(1,:) = 1.
    F(2,:) = tau/this%R
    do i=3,6
      F(i,:) = F(i-1,:)*F(2,:)
    end do

    allocate(case_idx(size(tau)))
    case_idx = 0
    where(  (0<=tau) .and. (tau<=this%R/2.) )
      case_idx = 1
    else where ( (this%R/2.<tau) .and. (tau<=this%R) )
      case_idx = 2
    else where ( (this%R<tau) .and. (tau<=3./2.*this%R) )
      case_idx = 3
    else where ( (3/2*this%R<tau) .and. (tau<=2*this%R) )
      case_idx = 4
    end where

    do i=1,size(tau)
      c(i) = 0.
      if(case_idx(i)>0) then
        ! scalar product: < B(case_idx(i),:),  F >
        do j=1,6
          c(i) = c(i) + this%B(case_idx(i),j) * F(j,i)
        end do
        if(case_idx(i)>1)then
          c(i) = c(i) + this%d(case_idx(i))*(this%R/tau(i))
        end if
      end if
    end do

    deallocate(F)
    deallocate(case_idx)

  end subroutine

  subroutine gaspari_cohn_1_evaluate_scalar(this,tau,c)

    implicit none

    ! arguments
    class(gaspari_cohn_1) :: this
    real, intent(in) :: tau
    real, intent(out) :: c

    ! local
    real, dimension(6) :: F
    integer :: i
    integer :: case_idx

    F(1) = 1.
    F(2) = tau/this%R
    do i=3,6
      F(i) = F(i-1)*F(2)
    end do

    case_idx = 0
    if(  (0<=tau) .and. (tau<=this%R/2.) ) then
      case_idx = 1
    else if ( (this%R/2.<tau) .and. (tau<=this%R) ) then
      case_idx = 2
    else if ( (this%R<tau) .and. (tau<=3./2.*this%R) ) then
      case_idx = 3
    else if ( (3/2*this%R<tau) .and. (tau<=2*this%R) ) then
      case_idx = 4
    end if

    c = 0.
    if(case_idx>0) then
      ! TODO use ddot?
      do i=1,6
        c = c + this%B(case_idx,i) * F(i)
      end do
      if(case_idx>1) then
        ! d(1) is always 0 and In case tau==0 this would result in invalid division
        c = c + this%d(case_idx)*(this%R/tau)
      end if
    end if

  end subroutine


end module
