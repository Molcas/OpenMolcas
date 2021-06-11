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

module fmm_scale_T_buffer

use fmm_global_paras, only: INTK, REALK, T_pair_batch, T_pair_single, Zero, One
use fmm_stats, only: stat_tpack_chunks, stat_tpack_unique, stat_tpack_total
use fmm_utils, only: fmm_quit

implicit none
private
! public procedures
public :: fmm_init_scale_T_buffer, &
          fmm_free_scale_T_buffer, &
          fmm_scale_T_buffer_add

integer(INTK), parameter :: BUFFER_SIZE = 500000
! module wide variables
type(T_pair_batch), save :: T_pair_buffer

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_scale_T_buffer()

  implicit none

  nullify(T_pair_buffer%items)
  allocate(T_pair_buffer%items(BUFFER_SIZE))
  T_pair_buffer%ndim = 0

end subroutine fmm_init_scale_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_free_scale_T_buffer(T_contractor)

  implicit none
  external :: T_contractor

  if (.not. associated(T_pair_buffer%items)) call fmm_quit('T_pair_buffer not alloc.')
  if (T_pair_buffer%ndim /= 0) then
    call expunge_scale_buffer(T_contractor)
    T_pair_buffer%ndim = 0
  end if
  deallocate(T_pair_buffer%items)
  nullify(T_pair_buffer%items)

end subroutine fmm_free_scale_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_scale_T_buffer_add(T_contractor,T_pair)

  implicit none
  type(T_pair_single), intent(in) :: T_pair
  external                        :: T_contractor
  real(REALK) :: ratio

  stat_tpack_total = stat_tpack_total+one
  T_pair_buffer%ndim = T_pair_buffer%ndim+1
  T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

  ! normalise T-vectors held in buffer for sorting purposes
  ratio = sqrt(sum(T_pair%r_ab*T_pair%r_ab))
  ! we sort wrt x-axis first, and we want opposite vectors to be together
  if (T_pair%r_ab(1) < zero) ratio = -ratio
  T_pair_buffer%items(T_pair_buffer%ndim)%paras%ratio = ratio
  T_pair_buffer%items(T_pair_buffer%ndim)%r_ab = T_pair%r_ab/ratio

  if (T_pair_buffer%ndim == BUFFER_SIZE) then
    ! sort the buffer and pass all T-pairs to contractor
    call expunge_scale_buffer(T_contractor)
  end if

end subroutine fmm_scale_T_buffer_add

!-------------------------------------------------------------------------------

subroutine expunge_scale_buffer(T_contractor)

  use fmm_sort_T_pairs, only: fmm_quicksort_wrt_vector, &
                              fmm_quicksort_wrt_ratio

  implicit none
  external :: T_contractor

  type(T_pair_batch) :: ptr, ptr2
  integer(INTK) :: i, lo, hi
  real(REALK) :: q1, q2

  ptr%ndim = min(BUFFER_SIZE,T_pair_buffer%ndim)
  ptr%items => T_pair_buffer%items(1:ptr%ndim)

  ! recursively sort wrt T-vectors, starting with x-component
  call sort_wrt_axis(1,ptr%items)

  ! expunge in batches of the same T-vector direction
  lo = 1
  do i=2,ptr%ndim
    q1 = ptr%items(i)%r_ab(1)
    q2 = ptr%items(i-1)%r_ab(1)
    if (q1 == q2) then
      q1 = ptr%items(i)%r_ab(2)
      q2 = ptr%items(i-1)%r_ab(2)
      if (q1 == q2) then
        q1 = ptr%items(i)%r_ab(3)
        q2 = ptr%items(i-1)%r_ab(3)
        if (q1 == q2) cycle
      end if
    end if
    hi = i-1
    ptr2%ndim = hi-lo+1
    ptr2%items => ptr%items(lo:hi)
    stat_tpack_unique = stat_tpack_unique+one
    call T_contractor(ptr2)
    lo = i
  end do

  ! finally do last batch
  hi = ptr%ndim
  ptr2%ndim = hi-lo+1
  ptr2%items => ptr%items(lo:hi)
  stat_tpack_unique = stat_tpack_unique+one
  call T_contractor(ptr2)

  T_pair_buffer%ndim = 0
  stat_tpack_chunks = stat_tpack_chunks+one

contains

  !-------------------------------------------------------------------------------

recursive subroutine sort_wrt_axis(xyz,items)

  implicit none
  integer(INTK), intent(in)          :: xyz
  type(T_pair_single), intent(inout) :: items(:)

  integer(INTK) :: i, lo, hi
  real(REALK) :: q1, q2

  if (size(items) == 1) return

  ! sort only if needed
  q1 = items(1)%r_ab(xyz)
  do i=2,size(items)
    q2 = items(i)%r_ab(xyz)
    if (q2 < q1) then
      call fmm_quicksort_wrt_vector(items,xyz)
      exit
    end if
    q1 = q2
  end do

  ! sub-sort next T-vector component
  lo = 1
  do i=2,size(items)
    q1 = items(i-1)%r_ab(xyz)
    q2 = items(i)%r_ab(xyz)
    if (q2 /= q1) then
      hi = i-1
      if (xyz == 3) then
        call fmm_quicksort_wrt_ratio(items(lo:hi))
        return
      else
        call sort_wrt_axis(xyz+1,items(lo:hi))
      end if
      lo = i
    end if
  end do

  ! do last batch
  hi = size(items)
  if (xyz == 3) then
    call fmm_quicksort_wrt_ratio(items(lo:hi))
    return
  else
    call sort_wrt_axis(xyz+1,items(lo:hi))
  end if

end subroutine sort_wrt_axis

!-------------------------------------------------------------------------------

end subroutine expunge_scale_buffer

!-------------------------------------------------------------------------------

end module fmm_scale_T_buffer
