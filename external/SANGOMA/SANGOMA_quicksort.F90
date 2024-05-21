module SANGOMA_quicksort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-11-19 17:57:23 vetra-carvalho>
!!!
!!!    Subsection of subroutines needed for the equivalent weights 
!!!    particle filter code supplied for SANGOMA project by 
!!!    Sanita Vetra-Carvalho. 
!!!
!!!    This code was taken from http://rosettacode.org/wiki/Quicksort#Fortran
!!!    and is distributed under GNU Free Documentation License 1.2.
!!!    see  http://www.gnu.org/licenses/fdl-1.2.html and was modified to also return
!!!    sorted index of the original array a.
!!!
!!!    Collection of subroutines to sort and return a one-dimensional array
!!!    as well as corresponding sorted index of the array a. 
!!!    Copyright (C) 2014  S. Vetra-Carvalho
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: s.vetra-carvalho @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!> subroutine to sort using the quicksort algorithm
!! @param[in,out] a, an array of doubles to be sorted
!! @param[in] na, dimension of the array a

  recursive subroutine quicksort_d(a,na)

    implicit none 

    ! DUMMY ARGUMENTS
    integer, intent(in) :: na ! nr or items to sort
    real, dimension(nA), intent(inout) :: a ! vector to be sorted
 
    ! LOCAL VARIABLES
    integer :: left, right, mid
    real :: pivot, temp
    integer :: marker

    if (nA > 1) then
       ! insertion sort limit of 47 seems best for sorting 10 million
       ! integers on Intel i7-980X CPU.  Derived data types that use
       ! more memory are optimized with smaller values - around 20 for a 16
       ! -byte type.
       if (nA > 47) then
          ! Do quicksort for large groups
          ! Get median of 1st, mid, & last points for pivot (helps reduce
          ! long execution time on some data sets, such as already
          ! sorted data, over simple 1st point pivot)
          mid = (nA+1)/2
          if (a(mid) >= a(1)) then
             if (a(mid) <= a(nA)) then
                pivot = a(mid)
             else if (a(nA) > a(1)) then
                pivot = a(nA)
             else
                pivot = a(1)
             end if
          else if (a(1) <= a(nA)) then
             pivot = a(1)
          else if (a(nA) > a(mid)) then
             pivot = a(nA)
          else
             pivot = a(mid)
          end if
 
          left = 0
          right = nA + 1
 
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
             end if
          end do
 
          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if
 
          call quicksort_d(A(:marker-1),marker-1)
          call quicksort_d(A(marker:),nA-marker+1)
 
       else
          call InsertionSort_d(A,nA)    ! Insertion sort for small groups is faster than Quicksort
       end if
    end if
 
  end subroutine quicksort_d


!> subroutine to sort using the insertionsort algorithm and return indecies
!! @param[in,out] a, an array of doubles to be sorted
!! @param[in] na, dimension of the array a 
  subroutine InsertionSort_d(a,na)
 
    ! DUMMY ARGUMENTS
    integer, intent(in) :: na                      !< nr or items to sort
    real, dimension(nA), intent(inout) :: a        !< vector to be sorted
 
    ! LOCAL VARIABLES
    real :: temp
    integer :: i, j
 
    do i = 2, nA
       j = i - 1
       temp = A(i)
       do
          if (j == 0) exit
          if (a(j) <= temp) exit
          A(j+1) = A(j)
          j = j - 1
       end do
       a(j+1) = temp
    end do

  end subroutine InsertionSort_d

!---------------------------------------------------------------------------
!> Quicksort for real and index array vectors     
!!
  recursive subroutine quicksort_idx_d(a,idx_a,na)

    implicit none 

    ! DUMMY ARGUMENTS
    integer, intent(in) :: na                      !< nr or items to sort
    real, dimension(nA), intent(inout) :: a        !< vector to be sorted
    integer, dimension(nA), intent(inout) :: idx_a !< sorted indecies of a
 
    ! LOCAL VARIABLES
    integer :: left, right, mid
    real :: pivot, temp
    integer :: marker, idx_temp
    integer :: i ! counter


    ! If this is the original call of the quicksort_d function 
    ! assign indecies to the array that we are sorting
    if (sum(idx_a) .eq. 0) then
       do i = 1,na
          idx_a(i) = 1
       end do
    end if

    if (nA > 1) then
       ! insertion sort limit of 47 seems best for sorting 10 million
       ! integers on Intel i7-980X CPU.  Derived data types that use
       ! more memory are optimized with smaller values - around 20 for a 16
       ! -byte type.
       if (nA > 47) then
          ! Do quicksort for large groups
          ! Get median of 1st, mid, & last points for pivot (helps reduce
          ! long execution time on some data sets, such as already
          ! sorted data, over simple 1st point pivot)
          mid = (nA+1)/2
          if (a(mid) >= a(1)) then
             if (a(mid) <= a(nA)) then
                pivot = a(mid)
             else if (a(nA) > a(1)) then
                pivot = a(nA)
             else
                pivot = a(1)
             end if
          else if (a(1) <= a(nA)) then
             pivot = a(1)
          else if (a(nA) > a(mid)) then
             pivot = a(nA)
          else
             pivot = a(mid)
          end if

          left = 0
          right = nA + 1
 
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp = A(left)
                idx_temp = idx_a(left)
                A(left) = A(right)
                idx_a(left) = idx_a(right)
                A(right) = temp
                idx_a(right) = idx_temp
             end if
          end do
 
          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if
 
          call quicksort_idx_d(A(:marker-1),idx_A(:marker-1),marker-1)
          call quicksort_idx_d(A(marker:),idx_A(marker:),nA-marker+1)
 
       else
          call InsertionSort_idx_d(A,idx_a,nA)    ! Insertion sort for small groups is faster than Quicksort
       end if
    end if
 
  end subroutine quicksort_idx_d


!> subroutine to sort using the insertionsort algorithm and return indecies
!! @param[in,out] a, an array of doubles to be sorted
!! @param[in,out] idx_a, an array of integers of sorted indecies
!! @param[in] na, dimension of the array a 
  subroutine InsertionSort_idx_d(a,idx_a,na)
 
     ! DUMMY ARGUMENTS
    integer, intent(in) :: na                      !< nr or items to sort
    real, dimension(nA), intent(inout) :: a        !< vector to be sorted
    integer, dimension(nA), intent(inout) :: idx_a !< sorted indecies of a
 
    ! LOCAL VARIABLES
    real :: temp
    integer :: i, j
 
    do i = 2, nA
       j = i - 1
       temp = A(i)
       do
          if (j == 0) exit
          if (a(j) <= temp) exit
          A(j+1) = A(j)
          idx_a(j+1) = idx_a(j)
          j = j - 1
       end do
       a(j+1) = temp
       idx_a(j+1) = i
    end do

  end subroutine InsertionSort_idx_d

end module SANGOMA_quicksort
