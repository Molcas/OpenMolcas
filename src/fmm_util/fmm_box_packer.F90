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

module fmm_box_packer

use fmm_global_paras, only: INTK, REALK, scheme_paras, raw_mm_paras, box_mm_paras
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_pkd_paras, &
          fmm_shift_and_pack_paras

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_pkd_paras(deepest_level,scheme,raw_paras,pkd_paras)

   use fmm_box_utils, only: fmm_box_centre, fmm_parent_box, fmm_grain

   implicit none
   integer(INTK), intent(in)         :: deepest_level
   type(scheme_paras), intent(in)    :: scheme
   type(raw_mm_paras), intent(inout) :: raw_paras(:)
   type(box_mm_paras), pointer       :: pkd_paras(:)

   type(box_mm_paras) :: tmp_paras(size(raw_paras))
   integer(INTK)      :: i, l_up, box(3), foo
   integer(INTK)      :: tmp_map(size(raw_paras))
   real(REALK)        :: grain

  l_up = deepest_level-1
  grain = fmm_grain(scheme,l_up)
  do i=1,size(raw_paras)
    tmp_paras(i)%box(:) = raw_paras(i)%box(:)
    tmp_paras(i)%cntr(:) = raw_paras(i)%box_cntr(:)
    tmp_paras(i)%bra = raw_paras(i)%bra
    tmp_paras(i)%level = deepest_level
    tmp_paras(i)%id = i
    box(:) = fmm_parent_box(tmp_paras(i)%box(:))
    tmp_paras(i)%cntr_up(:) = fmm_box_centre(box,grain)
    tmp_paras(i)%map_up = 0   ! not defined yet
  end do

  if (associated(pkd_paras)) call fmm_quit('paras should be nullified!')
  if (.false.) then
    ! just build unpacked paras
    foo = size(raw_paras)
    allocate(pkd_paras(foo))
    pkd_paras(:) = tmp_paras(:)
    tmp_map = [(i,i=1,size(raw_paras))]
  else
    ! combine paras in same box and branch (i.e. packed form)
    call pack_boxed_paras(tmp_paras,pkd_paras,tmp_map)
  end if
  ! store map of packing indices (raw:boxed)
  raw_paras(:)%map_up = tmp_map(:)

end subroutine fmm_init_pkd_paras

!-------------------------------------------------------------------------------

subroutine fmm_shift_and_pack_paras(level,scheme,paras_in,paras_out)

  use fmm_box_utils, only: fmm_box_centre,fmm_grain, &
                           fmm_parent_box,fmm_parent_bra

  implicit none
  integer(INTK), intent(in)         :: level
  type(scheme_paras), intent(in)    :: scheme
  type(box_mm_paras), intent(inout) :: paras_in(:)
  type(box_mm_paras), pointer       :: paras_out(:)

  type(box_mm_paras) :: tmp_paras(size(paras_in))
  integer(INTK)      :: i, l_up, box(3), foo
  integer(INTK)      :: tmp_map(size(paras_in))
  real(REALK)        :: grain, grain_up

  ! build tmp array for unpacked paras at next level up
  l_up = level-1
  grain = fmm_grain(scheme,level)
  grain_up = fmm_grain(scheme,l_up)
  do i=1,size(paras_in)
    tmp_paras(i)%box = fmm_parent_box(paras_in(i)%box)
    tmp_paras(i)%cntr(:) = fmm_box_centre(tmp_paras(i)%box,grain)
    tmp_paras(i)%bra = fmm_parent_bra(paras_in(i)%bra)
    tmp_paras(i)%level = level
    tmp_paras(i)%id = i
    box(:) = fmm_parent_box(tmp_paras(i)%box(:))
    tmp_paras(i)%cntr_up(:) = fmm_box_centre(box,grain_up)
    tmp_paras(i)%map_up = 0   ! not defined yet
  end do

  if (.false.) then
    ! just build unpacked paras
    foo = size(paras_in)
    allocate(paras_out(foo))
    paras_out(:) = tmp_paras(:)
    tmp_map = [(i,i=1,size(paras_in))]
  else
    ! combine paras in same box and branch (i.e. packed form)
    call pack_boxed_paras(tmp_paras,paras_out,tmp_map)
  end if
  ! store map of packing indices (raw:boxed)
  paras_in(:)%map_up = tmp_map(:)

end subroutine fmm_shift_and_pack_paras

!-------------------------------------------------------------------------------

subroutine pack_boxed_paras(paras_in,paras_out,map)

  implicit none
  type(box_mm_paras), intent(inout) :: paras_in(:)
  type(box_mm_paras), pointer       :: paras_out(:)
  integer(INTK), intent(out)        :: map(:)

  type(box_mm_paras) :: tmp_paras(size(paras_in))
  integer(INTK) :: i, k

  call fmm_sort_wrt_boxes_and_branches(1_INTK,paras_in)

  map(:) = 0
  tmp_paras(1) = paras_in(1)
  tmp_paras(1)%id = 1
  map(paras_in(1)%id) = 1

  k = 1
  do i=2,size(paras_in)

    if (paras_in(i)%box(3) /= paras_in(i-1)%box(3)) then
      k=k+1
    else if (paras_in(i)%box(2) /= paras_in(i-1)%box(2)) then
      k=k+1
    else if (paras_in(i)%box(1) /= paras_in(i-1)%box(1)) then
      k=k+1
    else if (paras_in(i)%bra /= paras_in(i-1)%bra) then
      k=k+1
    end if

    tmp_paras(k) = paras_in(i)
    tmp_paras(k)%id = k
    map(paras_in(i)%id) = k

  end do

  ! store packed parameters in exactly allocated array
  allocate(paras_out(k))
  paras_out = tmp_paras(1:k)

end subroutine pack_boxed_paras

!-------------------------------------------------------------------------------

recursive subroutine fmm_sort_wrt_boxes_and_branches(xyz,paras)

  use fmm_sort_paras, only: fmm_quicksort_wrt_boxes, &
                            fmm_quicksort_wrt_branches

  implicit none
  integer(INTK), intent(in)         :: xyz
  type(box_mm_paras), intent(inout) :: paras(:)

  integer(INTK) :: i, lo, hi
  real(REALK)   :: q1, q2

  if (size(paras)==1) return

  ! sort only if needed
  q1 = paras(1)%box(xyz)
  do i=2,size(paras)
    q2 = paras(i)%box(xyz)
    if (q2 < q1) then
      call fmm_quicksort_wrt_boxes(paras,xyz)
      exit
    end if
    q1 = q2
  end do

  ! sub-sort next box component
  lo = 1
  do i=2,size(paras)
    q1 = paras(i-1)%box(xyz)
    q2 = paras(i)%box(xyz)
    if (q2 /= q1) then
      hi = i-1
      if (xyz == 3) then
        call fmm_quicksort_wrt_branches(paras(lo:hi))
      else
        call fmm_sort_wrt_boxes_and_branches(xyz+1,paras(lo:hi))
      endif
      lo = i
    end if
  end do

  ! do last batch
  hi = size(paras)
  if (xyz == 3) then
    call fmm_quicksort_wrt_branches(paras(lo:hi))
  else
    call fmm_sort_wrt_boxes_and_branches(xyz+1,paras(lo:hi))
  end if

end subroutine fmm_sort_wrt_boxes_and_branches

!-------------------------------------------------------------------------------

endmodule fmm_box_packer
