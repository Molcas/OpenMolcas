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

module fmm_boundary

use fmm_global_paras, only: INTK, REALK, LUPRI, raw_mm_paras, fmm_planes, scheme_paras, Zero, One
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_opt_near_field

contains

!-------------------------------------------------------------------------------

subroutine verify_planes(grid,planes)

  implicit none
  type(raw_mm_paras), intent(in)  :: grid(:)
  type(fmm_planes), intent(inout) :: planes

  real(REALK), parameter :: THR=1.0e-15_REALK
  integer(INTK) :: i

  do i=1,size(grid)
    if (abs(grid(i)%cntr(1)-planes%xmin) > THR .and. &
        abs(grid(i)%cntr(1)-planes%xmax) > THR .and. &
        abs(grid(i)%cntr(2)-planes%ymin) > THR .and. &
        abs(grid(i)%cntr(2)-planes%ymax) > THR .and. &
        abs(grid(i)%cntr(3)-planes%zmin) > THR .and. &
        abs(grid(i)%cntr(3)-planes%zmax) > THR) then
      call fmm_quit('boundary planes not // to coordinate axes!')
    end if
  end do

end subroutine verify_planes

!-------------------------------------------------------------------------------

subroutine get_boundary_planes(grid,planes)

  implicit none
  type(raw_mm_paras), intent(in)  :: grid(:)
  type(fmm_planes), intent(inout) :: planes

  integer(INTK) :: i

  planes%xmin = zero
  planes%xmax = zero
  planes%ymin = zero
  planes%ymax = zero
  planes%zmin = zero
  planes%zmax = zero
  do i=1,size(grid)
    planes%xmin = min(planes%xmin,grid(i)%cntr(1))
    planes%xmax = max(planes%xmax,grid(i)%cntr(1))
    planes%ymin = min(planes%ymin,grid(i)%cntr(2))
    planes%ymax = max(planes%ymax,grid(i)%cntr(2))
    planes%zmin = min(planes%zmin,grid(i)%cntr(3))
    planes%zmax = max(planes%zmax,grid(i)%cntr(3))
  end do

end subroutine get_boundary_planes

!-------------------------------------------------------------------------------

subroutine get_closest_approach(RHS,planes,gap)

  implicit none
  type(raw_mm_paras), intent(in)  :: RHS(:)
  type(fmm_planes), intent(in)  :: planes
  real(REALK), intent(out) :: gap

  integer(INTK) :: i

  gap = 1.0e10_REALK
  do i=1,size(RHS)
    gap = min(gap,abs(RHS(i)%cntr(1)-planes%xmin))
    gap = min(gap,abs(RHS(i)%cntr(1)-planes%xmax))
    gap = min(gap,abs(RHS(i)%cntr(2)-planes%ymin))
    gap = min(gap,abs(RHS(i)%cntr(2)-planes%ymax))
    gap = min(gap,abs(RHS(i)%cntr(3)-planes%zmin))
    gap = min(gap,abs(RHS(i)%cntr(3)-planes%zmax))
  end do

end subroutine get_closest_approach

!-------------------------------------------------------------------------------

subroutine fmm_opt_near_field(scheme,LHS,RHS)

   use fmm_box_utils, only: fmm_deepest_level, fmm_branch, fmm_grain

  implicit none
  type(raw_mm_paras), intent(in)    :: LHS(:)
  type(raw_mm_paras), intent(in)    :: RHS(:)
  type(scheme_paras), intent(inout) :: scheme

  type(fmm_planes) :: planes
  real(REALK)      :: gap, dummy
  real(REALK)      :: grain
  integer(INTK)    :: branch

  ! We have only considered branch-free algorithm
  if (.not.scheme%branch_free) return

  call get_boundary_planes(LHS,planes)
  call verify_planes(LHS,planes)
  call get_closest_approach(RHS,planes,gap)
  write(LUPRI,'(a,es15.7)') ' Minimum distance to boundary =',gap

  if (gap < scheme%extent_min) then
    call fmm_quit('conflict between branch-free radius and boundary gap!')
  end if

  grain = fmm_grain(scheme,fmm_deepest_level(scheme))
  branch = fmm_branch(dummy,one/grain)

  ! we choose this condition to skip the NF conservatively
  if (gap > (branch+2)*grain) then
    write(LUPRI,*) 'There are no near-field interactions!'
    scheme%include_near_field = .false.
  end if

end subroutine fmm_opt_near_field

!-------------------------------------------------------------------------------

endmodule fmm_boundary
