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

module fmm_multi_T_buffer

use fmm_global_paras, only: INTK, T_pair_batch, T_pair_single, One
use fmm_stats, only: stat_tpack_chunks, stat_tpack_total
use fmm_utils, only: fmm_quit

implicit none
private
! public procedures
public :: fmm_init_multi_T_buffer, &
          fmm_free_multi_T_buffer, &
          fmm_multi_T_buffer_add

integer(INTK), parameter :: BUFFER_SIZE = 1000
! module wide variables
integer(INTK), save :: ndim_max
type(T_pair_batch), save :: T_pair_buffer

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_multi_T_buffer(ndim_max_in)

  implicit none
  integer(INTK), intent(in) :: ndim_max_in

  !ndim_max = ndim_max_in*2  ! we multiply by two for "paired" algorithm
  ndim_max = ndim_max_in
  if (ndim_max < 1) call fmm_quit('invalid multiple T-matrix dimension!')
  nullify(T_pair_buffer%items)
  allocate(T_pair_buffer%items(BUFFER_SIZE))
  T_pair_buffer%ndim = 0

end subroutine fmm_init_multi_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_free_multi_T_buffer(T_contractor)

  implicit none
  external T_contractor

  if (.not. associated(T_pair_buffer%items)) call fmm_quit('T_pair_buffer not alloc.')
  if (T_pair_buffer%ndim /= 0) then
    call expunge_multi_buffer(T_contractor)
    T_pair_buffer%ndim = 0
  end if
  deallocate(T_pair_buffer%items)
  nullify(T_pair_buffer%items)

end subroutine fmm_free_multi_T_buffer

!-------------------------------------------------------------------------------
!
!subroutine fmm_multi_T_buffer_add(T_contractor,T_pair)
!
!  implicit none
!  type(T_pair_single), intent(in) :: T_pair
!  external :: T_contractor
!
!  integer(INTK), save :: iRHS_last = 0
!  integer(INTK) :: iRHS
!
!  iRHS = T_pair%paras%RHS_id
!
!  !!!!if (BTEST(T_pair_buffer%ndim+1,0)) then
!    ! number of buffer entries is even; try to expunge
!    if ((T_pair_buffer%ndim == ndim_max) .or. ((iRHS /= iRHS_last) .and. (iRHS_last /= 0)) ) then
!      ! expunge buffer and build all the T-matrices at once
!      call T_contractor(T_pair_buffer)
!      T_pair_buffer%ndim = 0
!    end if
!    iRHS_last = T_pair%paras%RHS_id
!  !!!!end if
!
!  T_pair_buffer%ndim = T_pair_buffer%ndim +1
!  T_pair_buffer%items(T_pair_buffer%ndim) = T_pair
!
!end subroutine fmm_multi_T_buffer_add
!
!-------------------------------------------------------------------------------

subroutine fmm_multi_T_buffer_add(T_contractor,T_pair)

  implicit none
  type(T_pair_single), intent(in) :: T_pair
  external :: T_contractor

  if (T_pair_buffer%ndim == BUFFER_SIZE) then
    ! expunge buffer and build all the T-matrices at once
    call expunge_multi_buffer(T_contractor)
  end if

  stat_tpack_total = stat_tpack_total+one
  T_pair_buffer%ndim = T_pair_buffer%ndim+1
  T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

end subroutine fmm_multi_T_buffer_add

!-------------------------------------------------------------------------------

subroutine expunge_multi_buffer(T_contractor)

  use fmm_sort_T_pairs, only: fmm_quicksort_wrt_RHS

  implicit none
  external T_contractor

  integer(INTK) :: i, lo
  integer(INTK) :: iRHS, iRHS_next, item_max
  type(T_pair_batch) :: ptr

  lo = 1
  item_max = min((BUFFER_SIZE-1),(T_pair_buffer%ndim-1))

  ! sort only if needed
  iRHS = T_pair_buffer%items(1)%paras%RHS_id
  do i=2,item_max
    iRHS_next = T_pair_buffer%items(i)%paras%RHS_id
    if (iRHS_next < iRHS) then
      call fmm_quicksort_wrt_RHS(T_pair_buffer%items(1:item_max))
      exit
    end if
    iRHS = iRHS_next
  end do

  do i=1,item_max
    iRHS = T_pair_buffer%items(i)%paras%RHS_id
    iRHS_next = T_pair_buffer%items(i+1)%paras%RHS_id
    ptr%ndim=i-lo+1
    if ((iRHS /= iRHS_next) .or. (ptr%ndim == ndim_max)) then
      ptr%items => T_pair_buffer%items(lo:i)
      call T_contractor(ptr)
      lo = i+1
    end if
  end do

  ptr%ndim = (item_max+1)-lo+1
  ptr%items => T_pair_buffer%items(lo:(item_max+1))
  call T_contractor(ptr)

  T_pair_buffer%ndim = 0
  stat_tpack_chunks = stat_tpack_chunks+one

end subroutine expunge_multi_buffer

!-------------------------------------------------------------------------------

end module fmm_multi_T_buffer
