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

module fmm_W_buffer

use fmm_global_paras, only: T_pair_single, scheme_paras, SKIP_W_BUFFER, NULL_W_BUFFER, TREE_W_BUFFER, TREE_LENGTH
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_add_to_W_buffer, &
          fmm_open_W_buffer, &
          fmm_close_W_buffer

! diagnostic flag
character(len=4), save :: W_buffer_stat

contains

!-------------------------------------------------------------------------------

subroutine fmm_add_to_W_buffer(W_pair)

  implicit none
  type(T_pair_single), intent(in) :: W_pair
  external fmm_selected_w_buffer
  external fmm_selected_w_contractor

  call fmm_selected_w_buffer(fmm_selected_w_contractor,W_pair)

end subroutine fmm_add_to_W_buffer

!-------------------------------------------------------------------------------

subroutine fmm_close_W_buffer(scheme)

  use fmm_W_contractors, only: fmm_lock_W_con
  use fmm_tree_buffer, only: fmm_tree_buffer_finish

  implicit none
  type(scheme_paras), intent(in) :: scheme
  external fmm_selected_w_contractor

  if (W_buffer_stat /= 'OPEN') call fmm_quit('W_buffer already closed!')
  select case (scheme%W_con%W_buffer)
    case (SKIP_W_BUFFER)
      ! do nothing
    case (NULL_W_BUFFER)
      ! do nothing
    case (TREE_W_BUFFER)
      call fmm_tree_buffer_finish(fmm_selected_w_contractor)
    case default
      call fmm_quit('cannot reconcile list type in fmm_close_W_buffer')
  end select
  W_buffer_stat = 'FREE'
  fmm_lock_W_con = .false.

end subroutine fmm_close_W_buffer

!-------------------------------------------------------------------------------

subroutine fmm_null_W_buffer(W_contractor,W_pair)

  implicit none
  type(T_pair_single), intent(in) :: W_pair
  external W_contractor

  call W_contractor(W_pair)

end subroutine fmm_null_W_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

subroutine fmm_skip_W_buffer(W_contractor,W_pair)

# include "macros.fh"

  implicit none
  type(T_pair_single), intent(in) :: W_pair
  external W_contractor

  unused_var(W_contractor)
  unused_var(W_pair)

  return

end subroutine fmm_skip_W_buffer

!-------------------------------------------------------------------------------

subroutine fmm_open_W_buffer(scheme)

  use fmm_W_contractors, only: fmm_lock_W_con
  use fmm_tree_buffer, only: fmm_tree_buffer_init, &
                             fmm_tree_buffer_add

  implicit none
  type(scheme_paras), intent(in) :: scheme
  external fmm_store_w_buffer

  if (W_buffer_stat == 'OPEN') call fmm_quit('cannot reopen W_buffer')

  select case (scheme%W_con%W_buffer)
    case (SKIP_W_BUFFER)
      ! all W-contractions will be skipped by this choice of buffer
      call fmm_store_w_buffer(fmm_skip_W_buffer)
    case (NULL_W_BUFFER)
      call fmm_store_w_buffer(fmm_null_W_buffer)
    case (TREE_W_BUFFER)
      ! use tree-based sorting/evaluating module
      call fmm_store_w_buffer(fmm_tree_buffer_add)
      call fmm_tree_buffer_init(TREE_LENGTH,scheme%W_con%sort_para)
    case default
      call fmm_quit('cannot reconcile list type in fmm_open_W_buffer')
  end select

  W_buffer_stat = 'OPEN'
  fmm_lock_W_con = .true.

end subroutine fmm_open_W_buffer

!-------------------------------------------------------------------------------

end module fmm_W_buffer
