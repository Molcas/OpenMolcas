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

module fmm_box_utils

use fmm_global_paras, only: INTK, REALK, scheme_paras, gen_mm_paras, LHS_RHS_TYPE, box_mm_paras, MAX_LEVEL, TOP_LEVEL, WS_MIN, &
                            Two, Half
use fmm_stats, only: stat_level_saturation, stat_max_branch, stat_min_branch
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_box, &
          fmm_box_centre, &
          fmm_branch, &
          fmm_grain, &
          fmm_parent_box, &
          fmm_parent_bra, &
          fmm_deepest_level, &
          fmm_def_WS_NF, &
          fmm_def_WS_RFF, &
          fmm_same_box, &
          fmm_RFF_boxes, &
          fmm_NF_boxes, &
          fmm_translate_to_common_grid

contains

!-------------------------------------------------------------------------------

function fmm_parent_box(box)

  implicit none
  integer(INTK), intent(in) :: box(3)
  integer(INTK) :: fmm_parent_box(3)

  fmm_parent_box = 1+((box-1)/2)    ! get largest integer after /2

end function fmm_parent_box

!-------------------------------------------------------------------------------

function fmm_parent_bra(branch)

# include "macros.fh"

  implicit none
  integer(INTK), intent(in) :: branch
  integer(INTK) :: fmm_parent_bra

  unused_var(branch)

  !fmm_parent_bra = (((branch-1)/2)/2)*2 +2
  !fmm_parent_bra = max(WS_MIN,fmm_parent_bra)
  fmm_parent_bra = WS_MIN

end function fmm_parent_bra

!-------------------------------------------------------------------------------
! We use the coordinate shift to ensure all the box indices are positive

function fmm_box(centre,grain_inv)

  use fmm_qlm_builder, only: fmm_coord_shift

  implicit none
  real(REALK), intent(in) :: centre(3)
  real(REALK), intent(in) :: grain_inv
  integer(INTK) :: fmm_box(3)

  !fmm_box(:) = 1+int(centre(:)*grain_inv)
  fmm_box(:) = 1+int((centre(:)-fmm_coord_shift(:))*grain_inv)

end function fmm_box

!-------------------------------------------------------------------------------
! We take account of the coordinate shift used in fmm_box

function fmm_box_centre(box,grain)

  use fmm_qlm_builder, only: fmm_coord_shift

  implicit none
  integer(INTK), intent(in) :: box(3)
  real(REALK), intent(in)   :: grain
  real(REALK) :: fmm_box_centre(3)

  !fmm_box_centre(:) = grain*(box(:)-half)
  fmm_box_centre(:) = grain*(box(:)-half)+fmm_coord_shift(:)

end function fmm_box_centre

!-------------------------------------------------------------------------------

function fmm_branch(extent,grain_inv)

# include "macros.fh"

  implicit none
  !type(scheme_paras), intent(in) :: scheme
  real(REALK), intent(in) :: extent, grain_inv
  integer(INTK) :: fmm_branch

  unused_var(extent)
  unused_var(grain_inv)

  ! Note branches are always even in CFMM algorithm
  !if (scheme%branch_free) then
  !  fmm_branch = WS_MIN
  !else
  !  fmm_branch = 2*ceiling(extent*grain_inv)
  !end IF

  !fmm_branch = 2*ceiling(extent*grain_inv)
  !fmm_branch = max(WS_MIN, fmm_branch)

  fmm_branch = WS_MIN

  stat_max_branch = max(stat_max_branch,fmm_branch)
  stat_min_branch = min(stat_min_branch,fmm_branch)

end function fmm_branch

!-------------------------------------------------------------------------------

function fmm_deepest_level(scheme)

  use fmm_qlm_builder, only: fmm_system_size
  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK) :: fmm_deepest_level
  real(REALK) :: x

  x = fmm_system_size/scheme%grain
  ! Note we round UP the integer
  fmm_deepest_level = max(TOP_LEVEL,(1+int(log(x)/log(two))))

  if (fmm_deepest_level > MAX_LEVEL) then
    stat_level_saturation = fmm_deepest_level
    fmm_deepest_level = MAX_LEVEL
  end if

end function fmm_deepest_level

!-------------------------------------------------------------------------------
! Return the box size at a given level, based on the smallest box size.

function fmm_grain(scheme,level)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK), intent(in)      :: level
  real(REALK) :: fmm_grain
  integer(INTK) :: deepest_level

  deepest_level = fmm_deepest_level(scheme)
  fmm_grain = (scheme%grain)*(2**(deepest_level-level))

end function fmm_grain

!-------------------------------------------------------------------------------
! Routine to define the separation parameter for NF interactions
! i.e. boxes closer than this separation are treated as near field.
! only designed for cases where boxes are at the same level (same size)

subroutine fmm_def_WS_NF(LHS,RHS,id,WS_para)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  integer(INTK), intent(out)     :: WS_para

  WS_para = (LHS%box_paras(id%LHS)%bra+RHS%box_paras(id%RHS)%bra)/2

end subroutine fmm_def_WS_NF

!-------------------------------------------------------------------------------
! Routine to define local space that contains both LFF and NF space;
! i.e. boxes further apart than this separation are RFF boxes;
! note that this routine is generalised for cases where the boxes are
! at different levels (of different sizes) as is needed for the NlogN
! contraction scheme; in this case, we must choose which box size to express
! WS_para in terms of; by default we take the box size at highest level

subroutine fmm_def_WS_RFF(LHS_in,RHS_in,id,WS_para)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS_in, RHS_in
  type(LHS_RHS_type), intent(in) :: id
  integer(INTK), intent(out)     :: WS_para
  type(box_mm_paras) :: LHS, RHS

  LHS = LHS_in%box_paras(id%LHS)
  RHS = RHS_in%box_paras(id%RHS)

  ! we define WS_para in terms of the highest common grid
  call fmm_translate_to_common_grid(LHS,RHS)

  ! LFF is related to the parent's NF space
  LHS%bra = fmm_parent_bra(LHS%bra)
  RHS%bra = fmm_parent_bra(RHS%bra)

  ! hence NF-space for _parents_ of common grid
  WS_para = (LHS%bra+RHS%bra)/2
  ! hence LFF-space in common grid itself is
  WS_para = 1+2*WS_para
  ! +1 because of asymmetry in LFF space, and we are conservative.

end subroutine fmm_def_WS_RFF

!-------------------------------------------------------------------------------
! Are 2 boxes within their near-field space (hence too close to be
! treated by boxed multipole expansions exactly)?
! If number of boxes between A and B < ((Br(A)+Br(B))/2) they are NF

function fmm_NF_boxes(LHS,RHS)

  implicit none
  type(box_mm_paras), intent(in) :: LHS, RHS
  logical :: fmm_NF_boxes
  integer(INTK) :: WS_NF
  integer(INTK) :: space

  ! test only valid if box grids are at the same level depth
  if (LHS%level /= RHS%level) call fmm_quit('levels not equal in NF_boxes')

  WS_NF = (LHS%bra+RHS%bra)/2

  fmm_NF_boxes = .false.
  space = abs(LHS%box(3)-RHS%box(3))
  if (space <= WS_NF) then
    space = abs(LHS%box(2)-RHS%box(2))
    if (space <= WS_NF) then
      space = abs(LHS%box(1)-RHS%box(1))
      if (space <= WS_NF) then
        ! pair are in near field
        fmm_NF_boxes = .true.
      end if
    end if
  end if

end function fmm_NF_boxes

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Remote Far Field"
! based on separation.  Used by FMM and NlogN algorithms.

function fmm_RFF_boxes(LHS,RHS)

  implicit none
  type(box_mm_paras), intent(in) :: LHS, RHS
  logical :: fmm_RFF_boxes

  type(box_mm_paras) :: LHS_up, RHS_up

  ! test only valid if box grids are at the same level depth
  if (LHS%level /= RHS%level) call fmm_quit('levels in fmm_RFF_boxes')

  LHS_up = LHS
  RHS_up = RHS
  LHS_up%box = fmm_parent_box(LHS%box)
  LHS_up%bra = fmm_parent_bra(LHS%bra)
  RHS_up%box = fmm_parent_box(RHS%box)
  RHS_up%bra = fmm_parent_bra(RHS%bra)

  fmm_RFF_boxes = .not.(fmm_NF_boxes(LHS_up,RHS_up))

end function fmm_RFF_boxes

!-------------------------------------------------------------------------------

function fmm_same_box(LHS,RHS)

  implicit none
  type(box_mm_paras), intent(in) :: LHS, RHS
  logical :: fmm_same_box

  if (LHS%level /= RHS%level) call fmm_quit('levels not equal in same_box')
  fmm_same_box = .false.
  if (LHS%box(1) == RHS%box(1) .and. &
      LHS%box(2) == RHS%box(2) .and. &
      LHS%box(3) == RHS%box(3)) then
    fmm_same_box = .true.
  end if

end function fmm_same_box

!-------------------------------------------------------------------------------
! Routine to generate a temporary set of box parameters at a common grid
! level.  Parameters at the deeper level are translated to their parent's
! grid until the levels match.

subroutine fmm_translate_to_common_grid(LHS,RHS)

  implicit none
  type(box_mm_paras), intent(inout) :: LHS, RHS

  if (LHS%level == RHS%level) return
  if (LHS%level > RHS%level) then
    do while (LHS%level > RHS%level)
      LHS%box = fmm_parent_box(LHS%box)
      LHS%bra = fmm_parent_bra(LHS%bra)
      LHS%level = LHS%level-1
    end do
  else
    do while (RHS%level > LHS%level)
      RHS%box = fmm_parent_box(RHS%box)
      RHS%bra = fmm_parent_bra(RHS%bra)
      RHS%level = RHS%level-1
    end do
  end if

end subroutine fmm_translate_to_common_grid

!-------------------------------------------------------------------------------

end module fmm_box_utils
