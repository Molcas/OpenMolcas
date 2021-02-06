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

module fmm_box_builder

use fmm_global_paras, only: INTK, REALK, raw_mm_data, raw_mm_paras, scheme_paras, box_mm_data, box_mm_paras, TOP_LEVEL, Zero, One
use fmm_stats, only: stat_deepest_level, stat_RHS_boxes, stat_LHS_boxes
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_box_builder, &
          fmm_free_box_builder, &
          fmm_get_box_paras_at_level, &
          fmm_get_box_qlm_at_level

! Pointers to unboxed LHS parameters and RHS moments & parameters
type(raw_mm_data), pointer, save :: RHS_raw_mms
type(raw_mm_paras), pointer, save :: LHS_raw_paras(:)
! Packed LHS & RHS parameters and RHS moments at all levels
type(box_mm_data), pointer, save :: mms_at_lev(:)
! Deepest level of boxes in the hierarchy
integer(INTK), save :: deepest_level

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_box_builder(LHS_paras,RHS_mms,scheme)

  use fmm_box_utils, only: fmm_deepest_level

  implicit none
  type(raw_mm_paras), target, intent(inout) :: LHS_paras(:)
  type(raw_mm_data), target, intent(inout)  :: RHS_mms
  type(scheme_paras), intent(in)            :: scheme

  deepest_level = fmm_deepest_level(scheme)
  stat_deepest_level = deepest_level

  ! initialise RHS pointers to raw moments and parameters
  RHS_raw_mms => RHS_mms
  ! initialise LHS pointers to raw parameters
  LHS_raw_paras => LHS_paras(:)
  ! allocate foundation of hierarchical moment data structure
  call preallocate_levels()
  ! update raw_mm_paras with box and branch data
  call init_box_paras(LHS_paras,RHS_mms%paras,scheme)

end subroutine fmm_init_box_builder

!-------------------------------------------------------------------------------

subroutine init_box_paras(LHS,RHS,scheme)

  use fmm_box_utils, only: fmm_box, fmm_branch, fmm_grain, fmm_box_centre

  implicit none
  type(raw_mm_paras), intent(inout) :: LHS(:)
  type(raw_mm_paras), intent(inout) :: RHS(:)
  type(scheme_paras), intent(in)    :: scheme

  integer(INTK) :: i
  real(REALK) :: grain, grain_inv

  grain = fmm_grain(scheme,deepest_level)
  grain_inv = one/grain

  do i=1,size(RHS)
    RHS(i)%box = fmm_box(RHS(i)%cntr,grain_inv)
    RHS(i)%bra = fmm_branch(RHS(i)%ext,grain_inv)
    RHS(i)%box_cntr = fmm_box_centre(RHS(i)%box,grain)
    RHS(i)%map_up = 0     ! not defined yet
  end do

  do i=1,size(LHS)
    LHS(i)%box = fmm_box(LHS(i)%cntr,grain_inv)
    LHS(i)%bra = fmm_branch(LHS(i)%ext,grain_inv)
    LHS(i)%box_cntr = fmm_box_centre(LHS(i)%box,grain)
    LHS(i)%map_up = 0     ! not defined yet
  end do

end subroutine init_box_paras

!-------------------------------------------------------------------------------

subroutine build_paras_at_level(level,scheme)

  use fmm_box_packer, only: fmm_init_pkd_paras

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK), intent(in)      :: level

  type(raw_mm_paras), pointer :: ptr(:)

  if ((level < TOP_LEVEL) .or. (level > deepest_level)) then
    call fmm_quit('cannot iterate paras to this level!')
  end if

  ! build packed parameters at deepest level if not already done
  if (.not. associated(mms_at_lev(deepest_level)%RHS_paras)) then
    ptr => RHS_raw_mms%paras(:)
    call fmm_init_pkd_paras(deepest_level,scheme,ptr,mms_at_lev(deepest_level)%RHS_paras)
  end if
  if (.not. associated(mms_at_lev(deepest_level)%LHS_paras)) then
    ptr => LHS_raw_paras(:)
    call fmm_init_pkd_paras(deepest_level,scheme,ptr,mms_at_lev(deepest_level)%LHS_paras)
  end if

  if (level < deepest_level) then
    ! iterate paras through the box hierarchy to the reqd level
    call iterate_paras_to_level(level,scheme,'RHS')
    call iterate_paras_to_level(level,scheme,'LHS')
  end if

end subroutine build_paras_at_level

!-------------------------------------------------------------------------------

recursive subroutine iterate_paras_to_level(l,scheme,side)

  use fmm_box_packer, only: fmm_shift_and_pack_paras

  implicit none
  integer(INTK), intent(in)      :: l         ! level in box hierarchy
  type(scheme_paras), intent(in) :: scheme
  character(len=3), intent(in)   :: side

  type(box_mm_paras), pointer :: ptr(:)
  integer(INTK) :: l_down

  l_down = l+1
  select case (side)
    case ('RHS')
      if (.not. associated(mms_at_lev(l_down)%RHS_paras)) then
        ! must build paras at deeper levels before higher levels
        call iterate_paras_to_level(l_down,scheme,'RHS')
      end if
      ptr => mms_at_lev(l_down)%RHS_paras(:)
      ! build new paras from paras at previous level (including packing)
      call fmm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%RHS_paras)
    case ('LHS')
      if (.not. associated(mms_at_lev(l_down)%LHS_paras)) then
        ! must build paras at deeper levels before higher levels
        call iterate_paras_to_level(l_down,scheme,'LHS')
      end if
      ptr => mms_at_lev(l_down)%LHS_paras(:)
      ! build new paras from paras at previous level (including packing)
      call fmm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%LHS_paras)
    case default
      call fmm_quit('must build LHS or RHS paras!')
  end select

end subroutine iterate_paras_to_level

!-------------------------------------------------------------------------------

subroutine build_qlm_at_level(level,scheme,memory)

  use fmm_W_pair_builder, only: fmm_translate_raw_moments

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK), intent(in)      :: level
  character(len=4), intent(in)   :: memory

  type(raw_mm_data), pointer :: ptr1
  type(box_mm_data), pointer :: ptr2
  integer(INTK) :: mms_dim

  if ((level < TOP_LEVEL) .or. (level > deepest_level)) then
    call fmm_quit('cannot iterate boxed moments to this level!')
  end if

  ! build boxed moments at the deepest level if not already done
  if (.not. associated(mms_at_lev(deepest_level)%qlm_W)) then
    mms_dim = size(mms_at_lev(deepest_level)%RHS_paras)
    call allocate_lm_at_level(deepest_level,mms_dim,scheme%trans_LMAX)
    if (.not. associated(RHS_raw_mms)) call fmm_quit('mm_box_builder not correctly initialised!')
    ptr1 => RHS_raw_mms
    ptr2 => mms_at_lev(deepest_level)
    call fmm_translate_raw_moments(scheme,ptr1,ptr2)
  end if

  if (level < deepest_level) then
    ! iterate RHS moments through the box hierarchy to the reqd level
    call iterate_qlm_to_level(level,scheme,memory)
  end if

end subroutine build_qlm_at_level

!-------------------------------------------------------------------------------

recursive subroutine iterate_qlm_to_level(level,scheme,memory)

  use fmm_W_pair_builder, only: fmm_translate_boxed_moments

  implicit none
  integer(INTK), intent(in)      :: level
  type(scheme_paras), intent(in) :: scheme
  character(len=4), intent(in)   :: memory

  integer(INTK) :: mms_dim, l_down
  type(box_mm_data), pointer :: ptr1, ptr2

  l_down = level+1
  ! note we always build qlm_W first, then scale to give qlm_T
  if (.not. associated(mms_at_lev(l_down)%qlm_W)) then
    ! must have qlm at deeper levels before higher levels
    call iterate_qlm_to_level(l_down,scheme,memory)
  end if

  ! must build boxed paras before boxed moments
  if (.not. associated(mms_at_lev(l_down)%RHS_paras)) then
    call build_paras_at_level(l_down,scheme)
  end if

  mms_dim = size(mms_at_lev(level)%RHS_paras)
  call allocate_lm_at_level(level,mms_dim,scheme%trans_LMAX)
  ptr1 => mms_at_lev(l_down)
  ptr2 => mms_at_lev(level)
  call fmm_translate_boxed_moments(scheme,ptr1,ptr2)

end subroutine iterate_qlm_to_level

!-------------------------------------------------------------------------------

subroutine fmm_get_box_paras_at_level(l,scheme,box_paras,side)

  implicit none
  integer(INTK), intent(in)      :: l
  type(scheme_paras), intent(in) :: scheme
  type(box_mm_paras), pointer    :: box_paras(:)
  character(len=3), intent(in)   :: side

  if (.not. associated(mms_at_lev)) call fmm_quit('mms_at_lev should be allocated!')

  select case (side)
    case ('RHS')
      if (.not. associated(mms_at_lev(l)%RHS_paras)) then
        ! RHS paras not available at this level yet. So make them!
        call build_paras_at_level(l,scheme)
      end if
      box_paras => mms_at_lev(l)%RHS_paras(:)
      stat_RHS_boxes(l) = size(box_paras)
    case ('LHS')
      if (.not. associated(mms_at_lev(l)%LHS_paras)) then
        ! LHS paras not available at this level yet. So make them!
        call build_paras_at_level(l,scheme)
      end if
      box_paras => mms_at_lev(l)%LHS_paras(:)
      stat_LHS_boxes(l) = size(box_paras)
    case default
      call fmm_quit('must select just LHS or RHS paras to use')
  end select

end subroutine fmm_get_box_paras_at_level

!-------------------------------------------------------------------------------

subroutine fmm_get_box_qlm_at_level(l,scheme,qlm_T,side,memory)

  implicit none
  integer(INTK), intent(in)      :: l
  type(scheme_paras), intent(in) :: scheme
  real(REALK), pointer           :: qlm_T(:,:)
  character(len=3), intent(in)   :: side
  character(len=4), intent(in)   :: memory

  if (.not. associated(mms_at_lev)) call fmm_quit('mms_at_lev should be allocated!')
  if (.not. associated(mms_at_lev(l)%qlm_T)) then
    ! qlm data not available at this level yet. So make them!
    call build_qlm_at_level(l,scheme,memory)
  end if

  if (side == 'LHS') call fmm_quit('currently no LHS boxed mms built!')
  if (side == 'RHS') then
    qlm_T => mms_at_lev(l)%qlm_T(:,:)
  else
    call fmm_quit('must select LHS or RHS boxed moments!')
  end if

end subroutine fmm_get_box_qlm_at_level

!-------------------------------------------------------------------------------

subroutine preallocate_levels()

  implicit none
  integer(INTK) :: i

  if (deepest_level == 0) return
  if (associated(mms_at_lev)) call fmm_quit('mms_at_lev should not be allocated!')
  if (deepest_level < TOP_LEVEL) then
    call fmm_quit('error allocating levels in box hierarchy')
  end if
  ! We allow for possibility of all levels being used, but sub-variables
  ! are only allocated if required, so this is not a big deal
  allocate(mms_at_lev(deepest_level))
  do i=lbound(mms_at_lev,1),ubound(mms_at_lev,1)
    nullify(mms_at_lev(i)%LHS_paras)
    nullify(mms_at_lev(i)%RHS_paras)
    nullify(mms_at_lev(i)%qlm_W)
    nullify(mms_at_lev(i)%qlm_T)
  end do

end subroutine preallocate_levels

!-------------------------------------------------------------------------------

subroutine allocate_lm_at_level(l,mms_dim,LMAX)

  implicit none
  integer(INTK), intent(in) :: l, mms_dim, LMAX

  integer(INTK) :: lm_dim
  logical :: fail

  lm_dim = (1+LMAX)**2
  if (l > deepest_level) call fmm_quit('invalid level to allocate!')
  if (l < TOP_LEVEL) call fmm_quit('invalid level to allocate!')

  ! must test if pointers already allocated (compiler may not notice)
  fail = .false.
  if (associated(mms_at_lev(l)%qlm_T)) fail = .true.
  if (associated(mms_at_lev(l)%qlm_W)) fail = .true.
  if (fail) call fmm_quit('box lm data already allocated!')
  allocate(mms_at_lev(l)%qlm_T(lm_dim,mms_dim))
  allocate(mms_at_lev(l)%qlm_W(lm_dim,mms_dim))
  ! must zero as they are built additively
  mms_at_lev(l)%qlm_T = zero
  mms_at_lev(l)%qlm_W = zero

end subroutine allocate_lm_at_level

!-------------------------------------------------------------------------------

subroutine fmm_free_box_builder()

  implicit none
  integer(INTK) :: l

  nullify(RHS_raw_mms)
!FIXME
  nullify(LHS_raw_paras)

  if (associated(mms_at_lev)) then
    do l=lbound(mms_at_lev,1),ubound(mms_at_lev,1)
      if (associated(mms_at_lev(l)%LHS_paras,mms_at_lev(l)%RHS_paras)) then
        ! LHS and RHS paras the same and only point the same space
        deallocate(mms_at_lev(l)%RHS_paras)
      else
        if (associated(mms_at_lev(l)%RHS_paras)) deallocate(mms_at_lev(l)%RHS_paras)
        if (associated(mms_at_lev(l)%LHS_paras)) deallocate(mms_at_lev(l)%LHS_paras)
      end if
      if (associated(mms_at_lev(l)%qlm_W)) deallocate(mms_at_lev(l)%qlm_W)
      if (associated(mms_at_lev(l)%qlm_T)) deallocate(mms_at_lev(l)%qlm_T)
      nullify(mms_at_lev(l)%LHS_paras)
      nullify(mms_at_lev(l)%RHS_paras)
      nullify(mms_at_lev(l)%qlm_W)
      nullify(mms_at_lev(l)%qlm_T)
    end do
    deallocate(mms_at_lev)
  end if
  deepest_level = 0

end subroutine fmm_free_box_builder

!-------------------------------------------------------------------------------

end module fmm_box_builder
