!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module fmm_sort_T_pairs

! Sorting module with a very simple O(N^2) insertion sort
! algorithm or using the average O(NlogN) Quicksort algorithm
!
! Sort array of interaction T-pairs using Quicksort wrt
!     (a) RHS index
!     (b) T-vector component
!     (c) T-vector modulus (ratio)

use fmm_global_paras, only: INTK, REALK, T_pair_single
implicit none
private
! Public procedures
public :: fmm_quicksort_wrt_RHS, &
          fmm_quicksort_wrt_vector, &
          fmm_quicksort_wrt_ratio

! criteria for switching from Quicksort to Insertion-Sort algorithms
integer(INTK), parameter :: quicksort_CUTOFF = 10

contains

!-------------------------------------------------------------------------------

subroutine insertion_sort_wrt_RHS(arr)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  type(T_pair_single) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%paras%RHS_id > tmp%paras%RHS_id) then
        arr(j+1) = arr(j)
      else
        arr(j+1) = tmp
        cycle iloop
      end if
    end do
    arr(1) = tmp
  end do iloop

end subroutine insertion_sort_wrt_RHS

!-------------------------------------------------------------------------------

subroutine RHS_swap_elements(arr,i,j)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: i, j
  type(T_pair_single) :: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

end subroutine RHS_swap_elements

!-------------------------------------------------------------------------------
! Subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

function RHS_median_of_three(arr,left,right)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: left, right
  integer(INTK) :: RHS_median_of_three
  integer(INTK) :: cntr

  cntr = (left+right)/2

  if (arr(left)%paras%RHS_id > arr(cntr)%paras%RHS_id) call RHS_swap_elements(arr,left,cntr)
  if (arr(left)%paras%RHS_id > arr(right)%paras%RHS_id) call RHS_swap_elements(arr,left,right)
  if (arr(cntr)%paras%RHS_id > arr(right)%paras%RHS_id) call RHS_swap_elements(arr,cntr,right)

  call RHS_swap_elements(arr,cntr,right-1)
  RHS_median_of_three = arr(right-1)%paras%RHS_id

end function RHS_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_RHS(arr)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)

  integer(INTK) :: left, right
  integer(INTK) :: i, j, pivot

  left = 1
  right = size(arr)

  if (left+quicksort_CUTOFF > right) then
    ! array too small for Quicksort to be efficient, so don't use it
    call insertion_sort_wrt_RHS(arr)
    return
  end if

  pivot = RHS_median_of_three(arr,left,right)
  ! partition array elements about chosen pivot
  i = left
  j = right-2
  do ! until left, right pointers cross
    ! increment left until large element found
    do while(arr(i)%paras%RHS_id < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%paras%RHS_id > pivot)
      j = j-1
    end do
    if (i < j) then
      ! swap items (i,j)
      call RHS_swap_elements(arr,i,j)
      ! now increment (i,j) to avoid infinite loop if both = pivot
      i = i+1
      j = j-1
    else
      exit
    end if
  end do
  ! i is now the position where the pivot should be
  ! swap pivot back to middle of array
  call RHS_swap_elements(arr,i,right-1)

  ! now sort sub-arrays either side of pivot
  call fmm_quicksort_wrt_RHS(arr(left:(i-1)))
  call fmm_quicksort_wrt_RHS(arr((i+1):right))

end subroutine fmm_quicksort_wrt_RHS

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt 1 component of T-vector
!-----------------------------------------------------------------------

subroutine insertion_sort_wrt_vector(arr,xyz)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: xyz

  type(T_pair_single) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%r_ab(xyz) > tmp%r_ab(xyz)) then
        arr(j+1) = arr(j)
      else
        arr(j+1) = tmp
        cycle iloop
      end if
    end do
    arr(1) = tmp
  end do iloop

end subroutine insertion_sort_wrt_vector

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt 1 component of T-vector
!------------------------------------------------------------------

subroutine vector_swap_elements(arr,i,j)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: i, j
  type(T_pair_single) :: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

end subroutine vector_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

function vector_median_of_three(arr,left,right,xyz)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: left, right
  integer(INTK), intent(in)          :: xyz

  real(REALK) :: vector_median_of_three
  integer(INTK) :: cntr

  cntr = (left+right)/2

  if (arr(left)%r_ab(xyz) > arr(cntr)%r_ab(xyz)) call vector_swap_elements(arr,left,cntr)
  if (arr(left)%r_ab(xyz) > arr(right)%r_ab(xyz)) call vector_swap_elements(arr,left,right)
  if (arr(cntr)%r_ab(xyz) > arr(right)%r_ab(xyz)) call vector_swap_elements(arr,cntr,right)

  call vector_swap_elements(arr,cntr,right-1)
  vector_median_of_three = arr(right-1)%r_ab(xyz)

end function vector_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_vector(arr,xyz)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: xyz

  integer(INTK) :: left, right
  integer(INTK) :: i, j
  real(REALK) :: pivot

  left = 1
  right = size(arr)

  if (left+quicksort_CUTOFF > right) then
    ! array too small for Quicksort to be efficient, so don't use it
    call insertion_sort_wrt_vector(arr,xyz)
    return
  end if

  pivot = vector_median_of_three(arr,left,right,xyz)
  ! partition array elements about chosen pivot
  i = left
  j = right-2
  do ! until left, right pointers cross
    ! increment left until large element found
    do while(arr(i)%r_ab(xyz) < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%r_ab(xyz) > pivot)
      j = j-1
    end do
    if (i < j) then
      ! swap items (i,j)
      call vector_swap_elements(arr,i,j)
      ! now increment (i,j) to avoid infinite loop if both = pivot
      i = i+1
      j = j-1
    else
      exit
    end if
  end do
  ! i is now the position where the pivot should be
  ! swap pivot back to middle of array
  call vector_swap_elements(arr,i,right-1)

  ! now sort sub-arrays either side of pivot
  call fmm_quicksort_wrt_vector(arr(left:(i-1)),xyz)
  call fmm_quicksort_wrt_vector(arr((i+1):right),xyz)

end subroutine fmm_quicksort_wrt_vector

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt T-vector modulus (ratio)
!------------------------------------------------------------------------

subroutine insertion_sort_wrt_ratio(arr)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)

  type(T_pair_single) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%paras%ratio > tmp%paras%ratio) then
        arr(j+1) = arr(j)
      else
        arr(j+1) = tmp
        cycle iloop
      end if
    end do
    arr(1) = tmp
  end do iloop

end subroutine insertion_sort_wrt_ratio

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt T-vector modulus (ratio)
!-------------------------------------------------------------------

subroutine ratio_swap_elements(arr,i,j)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: i, j
  type(T_pair_single) :: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

end subroutine ratio_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

function ratio_median_of_three(arr,left,right)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)
  integer(INTK), intent(in)          :: left, right

  real(REALK) :: ratio_median_of_three
  integer(INTK) :: cntr

  cntr = (left+right)/2

  if (arr(left)%paras%ratio > arr(cntr)%paras%ratio) call ratio_swap_elements(arr,left,cntr)
  if (arr(left)%paras%ratio > arr(right)%paras%ratio) call ratio_swap_elements(arr,left,right)
  if (arr(cntr)%paras%ratio > arr(right)%paras%ratio) call ratio_swap_elements(arr,cntr,right)

  call ratio_swap_elements(arr,cntr,right-1)
  ratio_median_of_three = arr(right-1)%paras%ratio

end function ratio_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_ratio(arr)

  implicit none
  type(T_pair_single), intent(inout) :: arr(:)

  integer(INTK) :: left, right
  integer(INTK) :: i, j
  real(REALK) :: pivot

  left = 1
  right = size(arr)

  if (left+quicksort_CUTOFF > right) then
    ! array too small for Quicksort to be efficient, so don't use it
    call insertion_sort_wrt_ratio(arr)
    return
  end if

  pivot = ratio_median_of_three(arr,left,right)
  ! partition array elements about chosen pivot
  i = left
  j = right-2
  do ! until left, right pointers cross
    ! increment left until large element found
    do while(arr(i)%paras%ratio < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%paras%ratio > pivot)
      j = j-1
    end do
    if (i < j) then
      ! swap items (i,j)
      call ratio_swap_elements(arr,i,j)
      ! now increment (i,j) to avoid infinite loop if both = pivot
      i = i+1
      j = j-1
    else
      exit
    end if
  end do
  ! i is now the position where the pivot should be
  ! swap pivot back to middle of array
  call ratio_swap_elements(arr,i,right-1)

  ! now sort sub-arrays either side of pivot
  call fmm_quicksort_wrt_ratio(arr(left:(i-1)))
  call fmm_quicksort_wrt_ratio(arr((i+1):right))

end subroutine fmm_quicksort_wrt_ratio

!-------------------------------------------------------------------------------

end module fmm_sort_T_pairs
