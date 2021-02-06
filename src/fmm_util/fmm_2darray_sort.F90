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

module fmm_2darray_sort

! Sorting module with a very simple O(N^2) insertion sort algorithm
! for integer 2-d arrays like (3,:) according to input key 1-3

use fmm_global_paras, only: INTK
implicit none
private
! Public procedures
public :: fmm_insertion_sort

contains

!-------------------------------------------------------------------------------

subroutine fmm_insertion_sort(arr,xyz)

  implicit none
  integer(INTK), intent(inout) :: arr(:,:)
  integer(INTK), intent(in)    :: xyz

  integer(INTK) :: tmp(3)
  integer(INTK) :: i, j

  iloop: do i=1,size(arr,2)
    tmp(:) = arr(:,i)
    do j=(i-1),1,-1
      if (arr(xyz,j) > tmp(xyz)) then
        arr(:,j+1) = arr(:,j)
      else
        arr(:,j+1) = tmp(:)
        cycle iloop
      end if
    end do
    arr(:,1) = tmp(:)
  end do iloop

end subroutine fmm_insertion_sort

!-------------------------------------------------------------------------------

end module fmm_2darray_sort
