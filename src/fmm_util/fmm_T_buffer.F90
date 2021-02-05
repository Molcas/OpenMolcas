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

module fmm_T_buffer

use fmm_global_paras, only: INTK, T_pair_single, scheme_paras, SKIP_T_BUFFER, NULL_T_BUFFER, TREE_T_BUFFER, MULTI_T_BUFFER, &
                            SCALE_T_BUFFER, NEAR_FIELD, TMATM_DF, TREE_LENGTH, One
use fmm_stats, only: stat_T_mat_builds, stat_tpack_total, fmm_init_buffer_stats
use fmm_utils, only: fmm_quit

implicit none
private
! public procedures
public :: fmm_add_to_T_buffer, &
          fmm_open_T_buffer, &
          fmm_close_T_buffer

! only one T_buffer can be open at once
integer(INTK), save :: buffer = -1

! diagnostic flag
character(len=4), save :: T_buffer_stat

contains

!-------------------------------------------------------------------------------

subroutine fmm_add_to_T_buffer(T_pair)

  implicit none
  type(T_pair_single), intent(in) :: T_pair
  external fmm_selected_t_buffer
  external fmm_selected_t_contractor

  call fmm_selected_t_buffer(fmm_selected_t_contractor,T_pair)

end subroutine fmm_add_to_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_close_T_buffer()

  use fmm_T_contractors, only: fmm_lock_T_con
  use fmm_tree_buffer, only: fmm_tree_buffer_finish
  use fmm_multi_T_buffer, only: fmm_free_multi_T_buffer
  use fmm_scale_T_buffer, only: fmm_free_scale_T_buffer

  implicit none
  external fmm_selected_t_contractor

  if (T_buffer_stat /= 'OPEN') call fmm_quit('T_buffer already closed!')

  select case (buffer)
    case (SKIP_T_BUFFER)
      ! do nothing
    case (NULL_T_BUFFER)
      ! do nothing
    case (TREE_T_BUFFER)
      call fmm_tree_buffer_finish(fmm_selected_t_contractor)
    case (MULTI_T_BUFFER)
      call fmm_free_multi_T_buffer(fmm_selected_t_contractor)
    case (SCALE_T_BUFFER)
      call fmm_free_scale_T_buffer(fmm_selected_t_contractor)
    case default
      call fmm_quit('cannot reconcile list type in fmm_close_T_buffer')
  end select

  T_buffer_stat = 'FREE'
  fmm_lock_T_con = .false.

end subroutine fmm_close_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_null_T_buffer(T_contractor,T_pair)

  implicit none
  type(T_pair_single), intent(in) :: T_pair
  external T_contractor

  stat_tpack_total = stat_tpack_total+one
  call T_contractor(T_pair)

end subroutine fmm_null_T_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

subroutine fmm_skip_T_buffer(T_contractor,T_pair)

# include "macros.fh"

  implicit none
  type(T_pair_single), intent(in) :: T_pair
  external T_contractor

  unused_var(T_contractor)
  unused_var(T_pair)

  stat_T_mat_builds = stat_T_mat_builds+1
  return

end subroutine fmm_skip_T_buffer

!-------------------------------------------------------------------------------

subroutine fmm_open_T_buffer(scheme)

  use fmm_T_contractors, only: fmm_lock_T_con
  use fmm_tree_buffer, only: fmm_tree_buffer_init, &
                             fmm_tree_buffer_add
  use fmm_multi_T_buffer, only: fmm_init_multi_T_buffer, &
                                fmm_multi_T_buffer_add
  use fmm_scale_T_buffer, only: fmm_init_scale_T_buffer, &
                                fmm_scale_T_buffer_add

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK) :: sort_para
  external fmm_store_t_buffer

  call fmm_init_buffer_stats('T')
  if (T_buffer_stat == 'OPEN') call fmm_quit('cannot reopen T_buffer')

  if (scheme%phase == NEAR_FIELD) then
    buffer = scheme%T_con%NF_T_buffer
    sort_para = scheme%T_con%NF_sort_para
  else
    buffer = scheme%T_con%FF_T_buffer
    sort_para = scheme%T_con%FF_sort_para
  end if

  select case (buffer)
    case (SKIP_T_BUFFER)
      ! all T-contractions will be skipped by this choice of buffer
      call fmm_store_t_buffer(fmm_skip_T_buffer)
    case (NULL_T_BUFFER)
      call fmm_store_t_buffer(fmm_null_T_buffer)
    case (TREE_T_BUFFER)
      ! use tree-based sorting/evaluating module
      call fmm_store_t_buffer(fmm_tree_buffer_add)
      call fmm_tree_buffer_init(TREE_LENGTH,sort_para)
    case (SCALE_T_BUFFER)
      ! use tree-based sorting/evaluating module
      call fmm_store_t_buffer(fmm_scale_T_buffer_add)
      call fmm_init_scale_T_buffer()
    case (MULTI_T_BUFFER)
      ! use buffer to drive multiple T matrix simultaneous build
      call fmm_store_t_buffer(fmm_multi_T_buffer_add)
      call fmm_init_multi_T_buffer(TMATM_DF)
    case default
      call fmm_quit('cannot reconcile list type in fmm_open_T_buffer')
  end select

  T_buffer_stat = 'OPEN'
  fmm_lock_T_con = .true.

end subroutine fmm_open_T_buffer

!-------------------------------------------------------------------------------

end module fmm_T_buffer
