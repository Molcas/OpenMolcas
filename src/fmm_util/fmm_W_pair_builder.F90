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

module fmm_W_pair_builder

use fmm_global_paras, only: INTK,REALK, scheme_paras, old_new, T_pair_single, raw_mm_data, raw_mm_paras, box_mm_data, &
                            box_mm_paras,One
use fmm_stats, only: fmm_init_matrix_stats, fmm_init_buffer_stats
use fmm_W_contractors, only: fmm_select_W_con
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_translate_raw_moments, &
          fmm_translate_boxed_moments, &
          fmm_translate_parents_Vff, &
          fmm_get_raw_Vff_from_boxed_Vff

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_W_pair_builder(scheme)

  use fmm_W_buffer, only: fmm_open_W_buffer
  implicit none
  type(scheme_paras), intent(in) :: scheme
  call fmm_open_W_buffer(scheme)

end subroutine fmm_init_W_pair_builder

!-------------------------------------------------------------------------------

subroutine fmm_free_W_pair_builder(scheme)

  use fmm_W_buffer, only: fmm_close_W_buffer
  implicit none
  type(scheme_paras), intent(in) :: scheme
  call fmm_close_W_buffer(scheme)

end subroutine fmm_free_W_pair_builder

!-------------------------------------------------------------------------------

subroutine fmm_get_W_pair(addr,r_ab,new_LMAX,old_LMAX,object,W_pair)

  implicit none
  type(old_new), intent(in) :: addr
  real(REALK), intent(in) :: r_ab(3)
  integer(INTK), intent(in) :: new_LMAX, old_LMAX
  character(len=3), intent(in) :: object
  type(T_pair_single), intent(out) :: W_pair

  W_pair%paras%ratio = one ! i.e. W_pair%r_ab is the actual vector

  ! indices to map back to actual moments in separate array
  W_pair%paras%LHS_id = addr%new
  W_pair%paras%RHS_id = addr%old
  ! orders of contraction with T-matrix
  W_pair%paras%LHS_LMAX = new_LMAX
  W_pair%paras%RHS_LMAX = old_LMAX

  ! (see fmm_W_worker)
  select case(object)
    case('qlm')
      W_pair%r_ab(:) = r_ab(:)
      W_pair%N_or_T = 'N'
    case('Vff')
      ! For Vff translations, we require W'(-r_ab)
      W_pair%r_ab(:) = -r_ab(:)
      ! Use 'T' (transpose) for W_matrix contraction with DTRMV
      W_pair%N_or_T = 'T'
    case default
      call fmm_quit('cannot resolve translation object in fmm_get_W_pair!')
  end select

  W_pair%LMAX = max(W_pair%paras%LHS_LMAX,W_pair%paras%RHS_LMAX)
  W_pair%lm_max = (1+W_pair%LMAX)**2

end subroutine fmm_get_W_pair

!-------------------------------------------------------------------------------

subroutine fmm_translate_raw_moments(scheme,mms_in,mms_out)

  use fmm_W_contractors, only: fmm_set_W_con_ptrs
  use fmm_W_buffer, only: fmm_add_to_W_buffer
  use fmm_qlm_utils, only: fmm_get_T_sym_qlm

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  !FIXME: do not need to pass all this in! (just the moments)
  type(raw_mm_data), intent(in)    :: mms_in
  type(box_mm_data), intent(inout) :: mms_out

  type(T_pair_single) :: W_pair
  type(old_new) :: addr
  integer(INTK) :: i, new_LMAX, old_LMAX
  real(REALK) :: new_centre(3), old_centre(3), trans_vec(3)

  call fmm_select_W_con(scheme%W_con%ID)
  old_LMAX = scheme%raw_LMAX
  new_LMAX = scheme%trans_LMAX
  ! set pointers in W_contractor
  call fmm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

  ! generate W_pairs and pass them to the W_buffer
  call fmm_init_buffer_stats('W','RAW_BOX')
  call fmm_init_matrix_stats('W','RAW_BOX')
  call fmm_init_W_pair_builder(scheme)
  do i=1,size(mms_in%paras)
    ! indices to actual moments to translate from and to
    addr%old = mms_in%paras(i)%id
    addr%new = mms_in%paras(i)%map_up
    if (addr%new == 0) call fmm_quit('parameter mappings incomplete! 1')
    old_centre = mms_in%paras(i)%cntr
    new_centre = mms_in%paras(i)%box_cntr
    ! make W_pair and throw into W-buffer
    trans_vec = new_centre-old_centre
    call fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
    call fmm_add_to_W_buffer(W_pair)
  end do

  ! scale new translated moments for T_matrix symmetry
  ! but first ensure W_buffer is empty
  call fmm_free_W_pair_builder(scheme)
  call fmm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

end subroutine fmm_translate_raw_moments

!-------------------------------------------------------------------------------

subroutine fmm_translate_boxed_moments(scheme,mms_in,mms_out)

  use fmm_W_contractors, only: fmm_set_W_con_ptrs
  use fmm_W_buffer, only: fmm_add_to_W_buffer
  use fmm_qlm_utils, only: fmm_get_T_sym_qlm

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  type(box_mm_data), intent(in)    :: mms_in
  type(box_mm_data), intent(inout) :: mms_out

  type(T_pair_single) :: W_pair
  type(old_new) :: addr
  integer(INTK) :: i, new_LMAX, old_LMAX
  real(REALK) :: new_centre(3), old_centre(3), trans_vec(3)

  call fmm_select_W_con(scheme%W_con%ID)
  old_LMAX = scheme%trans_LMAX
  new_LMAX = scheme%trans_LMAX
  ! set pointers in W_contractor
  call fmm_set_W_con_ptrs(mms_in%qlm_W,mms_out%qlm_W)

  ! generate W_pairs and pass them to the W_buffer
  call fmm_init_buffer_stats('W','BOX_BOX')
  call fmm_init_matrix_stats('W','BOX_BOX')
  call fmm_init_W_pair_builder(scheme)
  do i=1,size(mms_in%RHS_paras)
    ! indices to actual moments to translate from and to
    addr%old = mms_in%RHS_paras(i)%id
    addr%new = mms_in%RHS_paras(i)%map_up
    if (addr%new == 0) call fmm_quit('parameter mappings incomplete! 2')
    old_centre = mms_in%RHS_paras(i)%cntr
    new_centre = mms_in%RHS_paras(i)%cntr_up
    ! make W_pair and throw into W-buffer
    trans_vec = new_centre-old_centre
    call fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'qlm',W_pair)
    call fmm_add_to_W_buffer(W_pair)
  end do

  ! scale new translated moments for T_matrix symmetry
  ! but first ensure W_buffer is empty
  call fmm_free_W_pair_builder(scheme)
  call fmm_get_T_sym_qlm(new_LMAX,mms_out%qlm_W,mms_out%qlm_T)

end subroutine fmm_translate_boxed_moments

!-------------------------------------------------------------------------------

subroutine fmm_translate_parents_Vff(level,scheme,Vff_p,Vff_c,c_box_paras)

  use fmm_W_contractors, only: fmm_set_W_con_ptrs
  use fmm_W_buffer, only: fmm_add_to_W_buffer

  implicit none
  integer(INTK), intent(in)      :: level
  type(scheme_paras), intent(in) :: scheme
  real(REALK), intent(in)        :: Vff_p(:,:)
  real(REALK), intent(inout)     :: Vff_c(:,:)
  type(box_mm_paras), intent(in) :: c_box_paras(:)   ! child

  type(T_pair_single) :: W_pair
  type(old_new) :: addr
  integer(INTK) :: i, new_LMAX, old_LMAX
  real(REALK) :: new_centre(3), old_centre(3), trans_vec(3)

  if (level <= 2) return

  new_LMAX = scheme%trans_LMAX
  old_LMAX = scheme%trans_LMAX

  call fmm_select_W_con(scheme%W_con%ID)
  call fmm_set_W_con_ptrs(Vff_p,Vff_c)
  call fmm_init_buffer_stats('W','BOX_BOX')
  call fmm_init_matrix_stats('W','BOX_BOX')
  call fmm_init_W_pair_builder(scheme)
  do i=1,size(Vff_c,2)
    addr%new = c_box_paras(i)%id
    addr%old = c_box_paras(i)%map_up
    if (addr%old == 0) call fmm_quit('parameter mappings incomplete! 3')
    new_centre = c_box_paras(i)%cntr
    old_centre = c_box_paras(i)%cntr_up
    trans_vec = new_centre-old_centre
    call fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
    call fmm_add_to_W_buffer(W_pair)
  end do
  call fmm_free_W_pair_builder(scheme)

end subroutine fmm_translate_parents_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_raw_Vff_from_boxed_Vff(raw_paras,scheme,boxed_Vff,Vff)

  use fmm_W_contractors, only: fmm_set_W_con_ptrs
  use fmm_W_buffer, only: fmm_add_to_W_buffer

  implicit none
  type(raw_mm_paras), intent(in) :: raw_paras(:)
  type(scheme_paras), intent(in) :: scheme
  real(REALK), intent(in)        :: boxed_Vff(:,:)
  real(REALK), intent(inout)     :: Vff(:,:)

  type(T_pair_single) :: W_pair
  type(old_new) :: addr
  integer(INTK) :: i, new_LMAX, old_LMAX
  real(REALK) :: new_centre(3), old_centre(3), trans_vec(3)

  new_LMAX = scheme%raw_LMAX
  old_LMAX = scheme%trans_LMAX

  call fmm_select_W_con(scheme%W_con%BR_ID)
  call fmm_set_W_con_ptrs(boxed_Vff,Vff)
  call fmm_init_buffer_stats('W','BOX_RAW')
  call fmm_init_matrix_stats('W','BOX_RAW')
  call fmm_init_W_pair_builder(scheme)
  do i=1,size(raw_paras)
    addr%new = raw_paras(i)%id
    addr%old = raw_paras(i)%map_up
    if (addr%old == 0) call fmm_quit('parameter mappings incomplete! 4')
    new_centre = raw_paras(i)%cntr
    old_centre = raw_paras(i)%box_cntr
    trans_vec = new_centre-old_centre
    call fmm_get_W_pair(addr,trans_vec,new_LMAX,old_LMAX,'Vff',W_pair)
    call fmm_add_to_W_buffer(W_pair)
  end do
  call fmm_free_W_pair_builder(scheme)

end subroutine fmm_get_raw_Vff_from_boxed_Vff

!-------------------------------------------------------------------------------

end module fmm_W_pair_builder
