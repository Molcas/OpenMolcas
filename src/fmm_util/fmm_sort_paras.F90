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

module fmm_sort_paras

! Sorting module with a very simple O(N^2) insertion sort
! algorithm or using the average O(NlogN) Quicksort algorithm
!
! Sort multipole moment parameters wrt
!     (a) expansion centre component
!     (b) box component
!     (c) branch

use fmm_global_paras, only: INTK, REALK, box_mm_paras, raw_mm_paras
implicit none
private
! Public procedures
public :: fmm_quicksort_wrt_branches, &
          fmm_quicksort_wrt_vector, &
          fmm_quicksort_wrt_boxes

! criteria for switching from Quicksort to Insertion-Sort algorithms
integer(INTK), parameter :: quicksort_CUTOFF = 10

contains

!-------------------------------------------------------------------------------

subroutine insertion_sort_wrt_branches(arr)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  type(box_mm_paras) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%bra > tmp%bra) then
        arr(j+1) = arr(j)
      else
        arr(j+1) = tmp
        cycle iloop
      end if
    end do
    arr(1) = tmp
  end do iloop

end subroutine insertion_sort_wrt_branches

!-------------------------------------------------------------------------------

subroutine bra_swap_elements(arr,i,j)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: i, j
  type(box_mm_paras) :: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

end subroutine bra_swap_elements

!-------------------------------------------------------------------------------
! Subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

function bra_median_of_three(arr,left,right)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: left, right
  integer(INTK) :: bra_median_of_three
  integer(INTK) :: cntr

  cntr = (left+right)/2
  if (arr(left)%bra > arr(cntr)%bra) call bra_swap_elements(arr,left,cntr)
  if (arr(left)%bra > arr(right)%bra) call bra_swap_elements(arr,left,right)
  if (arr(cntr)%bra > arr(right)%bra) call bra_swap_elements(arr,cntr,right)
  call bra_swap_elements(arr,cntr,right-1)
  bra_median_of_three = arr(right-1)%bra

end function bra_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_branches(arr)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)

  integer(INTK) :: left, right
  integer(INTK) :: i, j, pivot

  left = 1
  right = size(arr)

  if (left+quicksort_CUTOFF > right) then
    ! array too small for Quicksort to be efficient, so don't use it
    call insertion_sort_wrt_branches(arr)
    return
  end if

  pivot = bra_median_of_three(arr,left,right)
  ! partition array elements about chosen pivot
  i = left
  j = right-2
  do ! until left, right pointers cross
    ! increment left until large element found
    do while(arr(i)%bra < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%bra > pivot)
      j = j-1
    end do
    if (i < j) then
      ! swap items (i,j)
      call bra_swap_elements(arr,i,j)
      ! now increment (i,j) to avoid infinite loop if both = pivot
      i = i+1
      j = j-1
    else
      exit
    end if
  end do
  ! i is now the position where the pivot should be
  ! swap pivot back to middle of array
  call bra_swap_elements(arr,i,right-1)

  ! now sort sub-arrays either side of pivot
  call fmm_quicksort_wrt_branches(arr(left:(i-1)))
  call fmm_quicksort_wrt_branches(arr((i+1):right))

end subroutine fmm_quicksort_wrt_branches

!-------------------------------------------------------------------------------
! Insertion sort for sorting of multipole moment centres
!-----------------------------------------------------------------------

subroutine insertion_sort_wrt_vector(arr,xyz)

  implicit none
  type(raw_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: xyz

  type(raw_mm_paras) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%cntr(xyz) > tmp%cntr(xyz)) then
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
  type(raw_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: i, j
  type(raw_mm_paras) :: tmp

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
  type(raw_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: left, right
  integer(INTK), intent(in)         :: xyz

  real(REALK) :: vector_median_of_three
  integer(INTK) :: cntr

  cntr = (left+right)/2

  if (arr(left)%cntr(xyz) > arr(cntr)%cntr(xyz)) call vector_swap_elements(arr,left,cntr)
  if (arr(left)%cntr(xyz) > arr(right)%cntr(xyz)) call vector_swap_elements(arr,left,right)
  if (arr(cntr)%cntr(xyz) > arr(right)%cntr(xyz)) call vector_swap_elements(arr,cntr,right)

  call vector_swap_elements(arr,cntr,right-1)
  vector_median_of_three = arr(right-1)%cntr(xyz)

end function vector_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_vector(arr,xyz)

  implicit none
  type(raw_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: xyz

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
    do while(arr(i)%cntr(xyz) < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%cntr(xyz) > pivot)
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
! Insertion sort for sorting wrt one box index
!-----------------------------------------------------------------------

subroutine insertion_sort_wrt_boxes(arr,xyz)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: xyz

  type(box_mm_paras) :: tmp
  integer(INTK) :: i, j

  iloop: do i=1,size(arr)
    tmp = arr(i)
    do j=(i-1),1,-1
      if (arr(j)%box(xyz) > tmp%box(xyz)) then
        arr(j+1) = arr(j)
      else
        arr(j+1) = tmp
        cycle iloop
      end if
    end do
    arr(1) = tmp
  end do iloop

end subroutine insertion_sort_wrt_boxes

!-------------------------------------------------------------------------------
! Quicksort for sorting wrt one box index
!------------------------------------------------------------------

subroutine boxes_swap_elements(arr,i,j)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: i, j
  type(box_mm_paras) :: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

end subroutine boxes_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

function boxes_median_of_three(arr,left,right,xyz)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: left, right
  integer(INTK), intent(in)         :: xyz

  integer(INTK) :: boxes_median_of_three
  integer(INTK) :: box

  box = (left+right)/2

  if (arr(left)%box(xyz) > arr(box)%box(xyz)) call boxes_swap_elements(arr,left,box)
  if (arr(left)%box(xyz) > arr(right)%box(xyz)) call boxes_swap_elements(arr,left,right)
  if (arr(box)%box(xyz) > arr(right)%box(xyz)) call boxes_swap_elements(arr,box,right)

  call boxes_swap_elements(arr,box,right-1)
  boxes_median_of_three = arr(right-1)%box(xyz)

end function boxes_median_of_three

!-------------------------------------------------------------------------------

recursive subroutine fmm_quicksort_wrt_boxes(arr,xyz)

  implicit none
  type(box_mm_paras), intent(inout) :: arr(:)
  integer(INTK), intent(in)         :: xyz

  integer(INTK) :: left, right
  integer(INTK) :: i, j
  integer(INTK) :: pivot

  left = 1
  right = size(arr)

  if (left+quicksort_CUTOFF > right) then
    ! array too small for Quicksort to be efficient, so don't use it
    call insertion_sort_wrt_boxes(arr,xyz)
    return
  end if

  pivot = boxes_median_of_three(arr,left,right,xyz)
  ! partition array elements about chosen pivot
  i = left
  j = right-2
  do ! until left, right pointers cross
    ! increment left until large element found
    do while(arr(i)%box(xyz) < pivot)
      i = i+1
    end do
    ! increment right until small element found
    do while(arr(j)%box(xyz) > pivot)
      j = j-1
    end do
    if (i < j) then
      ! swap items (i,j)
      call boxes_swap_elements(arr,i,j)
      ! now increment (i,j) to avoid infinite loop if both = pivot
      i = i+1
      j = j-1
    else
      exit
    end if
  end do
  ! i is now the position where the pivot should be
  ! swap pivot back to middle of array
  call boxes_swap_elements(arr,i,right-1)

  ! now sort sub-arrays either side of pivot
  call fmm_quicksort_wrt_boxes(arr(left:(i-1)),xyz)
  call fmm_quicksort_wrt_boxes(arr((i+1):right),xyz)

end subroutine fmm_quicksort_wrt_boxes

!-------------------------------------------------------------------------------

end module fmm_sort_paras

!==============================================================================
!
!module fmm_sort_sh_pairs
!
!! Sorting module with a very simple O(N^2) insertion sort
!! algorithm or using the average O(NlogN) Quicksort algorithm
!!
!! Sort shell pair array wrt
!!     (a) box component
!
!use fmm_global_paras, only: INTK, fmm_shell_pair_node
!implicit none
!private
!! Public procedures
!public :: fmm_quicksort_wrt_boxes
!
!! criteria for switching from Quicksort to Insertion-Sort algorithms
!integer(INTK), parameter :: quicksort_CUTOFF = 10
!
!contains
!
!!-----------------------------------------------------------------------
!! Insertion sort for sorting wrt one box index
!!-----------------------------------------------------------------------
!
!subroutine insertion_sort_wrt_boxes(arr,xyz)
!
!  implicit none
!  type(fmm_shell_pair_node), intent(inout) :: arr(:)
!  integer(INTK), intent(in)                :: xyz
!
!  type(fmm_shell_pair_node) :: tmp
!  integer(INTK) :: i, j
!
!  iloop: do i=1,size(arr)
!    tmp = arr(i)
!    do j=(i-1),1,-1
!      if (arr(j)%box(xyz) > tmp%box(xyz)) then
!        arr(j+1) = arr(j)
!      else
!        arr(j+1) = tmp
!        cycle iloop
!      end if
!    end do
!    arr(1) = tmp
!  end do iloop
!
!end subroutine insertion_sort_wrt_boxes
!
!!------------------------------------------------------------------
!! Quicksort for sorting wrt one box index
!!------------------------------------------------------------------
!
!subroutine boxes_swap_elements(arr,i,j)
!
!  implicit none
!  type(fmm_shell_pair_node), intent(inout) :: arr(:)
!  integer(INTK), intent(in)                :: i, j
!  type(fmm_shell_pair_node) :: tmp
!
!  tmp = arr(i)
!  arr(i) = arr(j)
!  arr(j) = tmp
!
!end subroutine boxes_swap_elements
!
!!-------------------------------------------------------------------------------
!! subroutine to return the median of first, last, and centre array elements.
!! order these and hide the pivot (median) as the second to last array element
!! (acts as a sentinel whilst the last is the largest and put in correct place.)
!
!function boxes_median_of_three(arr,left,right,xyz)
!
!  implicit none
!  type(fmm_shell_pair_node), intent(inout) :: arr(:)
!  integer(INTK), intent(in)                :: left, right
!  integer(INTK),intent(in)                 :: xyz
!
!  integer(INTK) :: boxes_median_of_three
!  integer(INTK) :: box
!
!  box = (left+right)/2
!
!  if (arr(left)%box(xyz) > arr(box)%box(xyz)) call boxes_swap_elements(arr,left,box)
!  if (arr(left)%box(xyz) > arr(right)%box(xyz)) call boxes_swap_elements(arr,left,right)
!  if (arr(box)%box(xyz) > arr(right)%box(xyz)) call boxes_swap_elements(arr,box,right)
!
!  call boxes_swap_elements(arr,box,right-1)
!  boxes_median_of_three = arr(right-1)%box(xyz)
!
!end function boxes_median_of_three
!
!!-------------------------------------------------------------------------------
!
!recursive subroutine fmm_quicksort_wrt_boxes(arr,xyz)
!
!   implicit none
!   type(fmm_shell_pair_node), intent(inout) :: arr(:)
!   integer(INTK), intent(in)                :: xyz
!
!   integer(INTK) :: left, right
!   integer(INTK) :: i, j
!   integer(INTK) :: pivot
!
!   left = 1
!   right = size(arr)
!
!   if (left+quicksort_CUTOFF > right) then
!      ! array too small for Quicksort to be efficient, so don't use it
!      call insertion_sort_wrt_boxes(arr,xyz)
!      return
!   end if
!
!   pivot = boxes_median_of_three(arr,left,right,xyz)
!   ! partition array elements about chosen pivot
!   i = left
!   j = right-2
!   do ! until left, right pointers cross
!      ! increment left until large element found
!      do while (arr(i)%box(xyz) < pivot)
!         i = i+1
!      end do
!      ! increment right until small element found
!      do while (arr(j)%box(xyz) > pivot)
!         j = j-1
!      end do
!      if (i < j) then
!         ! swap items (i,j)
!         call boxes_swap_elements(arr,i,j)
!         ! now increment (i,j) to avoid infinite loop if both = pivot
!         i = i+1
!         j = j-1
!      else
!         exit
!      end if
!   end do
!   ! i is now the position where the pivot should be
!   ! swap pivot back to middle of array
!   call boxes_swap_elements(arr,i,right-1)
!
!   ! now sort sub-arrays either side of pivot
!   call fmm_quicksort_wrt_boxes(arr(left:(i-1)),xyz)
!   call fmm_quicksort_wrt_boxes(arr((i+1):right),xyz)
!
!end subroutine fmm_quicksort_wrt_boxes
!
!!-------------------------------------------------------------------------------
!
!end module fmm_sort_sh_pairs
