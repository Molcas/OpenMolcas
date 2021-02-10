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

module fmm_local_search

use fmm_global_paras, only: INTK, id_list, id_node, scheme_paras, box_mm_data, box_mm_paras, gen_mm_paras, LHS_raw_RHS_raw, &
                            LHS_box_RHS_box, TOP_LEVEL, DO_FMM
use fmm_box_utils, only: fmm_parent_box, fmm_same_box, fmm_RFF_boxes
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_init_local_search, &
          fmm_free_local_search, &
          fmm_get_local_paras

type fmm_local_map
  type(id_list), pointer :: box_list(:)
end type fmm_local_map

type(fmm_local_map), allocatable, save :: map_at_level(:)
integer(INTK), save :: deepest_level

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_local_search(scheme)

  use fmm_box_builder, only: fmm_get_box_paras_at_level
  use fmm_box_utils, only: fmm_deepest_level

  implicit none
  type(scheme_paras), intent(in) :: scheme

  type(fmm_local_map), allocatable :: RHS_map(:)
  type(box_mm_data), allocatable :: mm_paras(:)
  type(box_mm_paras), pointer :: ptr(:)
  integer(INTK) :: level, j

  if (.not. scheme%branch_free) call fmm_quit('local search module may not work with CFMM!')
  if (.not. scheme%algorithm == DO_FMM) call fmm_quit('local search module is only for FMM algorithm')

  deepest_level = fmm_deepest_level(scheme)
  nullify(ptr)
  allocate(mm_paras(deepest_level))
  allocate(RHS_map(deepest_level))

  ! Get multipole moment box information over whole hierarchy
  do level=TOP_LEVEL,deepest_level
    nullify(RHS_map(level)%box_list)
    nullify(mm_paras(level)%LHS_paras)
    nullify(mm_paras(level)%RHS_paras)
    nullify(mm_paras(level)%qlm_T)
    nullify(mm_paras(level)%qlm_W)
    call fmm_get_box_paras_at_level(level,scheme,ptr,'LHS')
    mm_paras(level)%LHS_paras => ptr
    call fmm_get_box_paras_at_level(level,scheme,ptr,'RHS')
    mm_paras(level)%RHS_paras => ptr
  end do

  ! Build parent:child map of RHS parameters and use it to make
  ! local interaction map of LHS-RHS pairs
  call fmm_init_local_map(mm_paras)
  call fmm_build_RHS_map(mm_paras,RHS_map)
  call fmm_build_local_map(mm_paras,RHS_map)

  ! Deallocations
  do level=TOP_LEVEL,deepest_level
    if (associated(RHS_map(level)%box_list)) then
      do j=1,size(RHS_map(level)%box_list)
        call free_linked_list(RHS_map(level)%box_list(j)%head)
      end do
      deallocate(RHS_map(level)%box_list)
    end if
    nullify(RHS_map(level)%box_list)
    nullify(mm_paras(level)%LHS_paras)
    nullify(mm_paras(level)%RHS_paras)
  end do
  deallocate(mm_paras,RHS_map)

end subroutine fmm_init_local_search

!-------------------------------------------------------------------------------

subroutine fmm_init_local_map(mm_paras)

  implicit none
  type(box_mm_data), intent(in) :: mm_paras(:)

  integer(INTK) :: i, j, ndim

  allocate(map_at_level(deepest_level))

  do i=TOP_LEVEL,deepest_level
    ndim = size(mm_paras(i)%LHS_paras)
    nullify(map_at_level(i)%box_list)
    allocate(map_at_level(i)%box_list(ndim))
    do j=1,ndim
      nullify(map_at_level(i)%box_list(j)%head)
      map_at_level(i)%box_list(j)%occ = 0
    end do
  end do

end subroutine fmm_init_local_map

!-------------------------------------------------------------------------------

subroutine fmm_free_local_search()

  implicit none
  integer(INTK) :: i, j, ndim

  if (deepest_level == TOP_LEVEL) return

  do i=TOP_LEVEL,deepest_level
    if (associated(map_at_level(i)%box_list)) then
      ndim = size(map_at_level(i)%box_list)
      do j=1,ndim
        call free_linked_list(map_at_level(i)%box_list(j)%head)
      end do
      deallocate(map_at_level(i)%box_list)
    end if
    nullify(map_at_level(i)%box_list)
  end do

  if (allocated(map_at_level)) deallocate(map_at_level)

end subroutine fmm_free_local_search

!-------------------------------------------------------------------------------

recursive subroutine free_linked_list(occ_node)

  implicit none
  type(id_node), pointer :: occ_node

  if (.not. associated(occ_node)) return

  ! traverse linked-list freeing from bottom up recursively
  if (associated(occ_node%next)) then
    call free_linked_list(occ_node%next)
    if (associated(occ_node)) deallocate(occ_node)
    nullify(occ_node)
  end if
  if (associated(occ_node)) deallocate(occ_node)
  nullify(occ_node)

end subroutine free_linked_list

!-------------------------------------------------------------------------------
! For each box, a list of the children boxes is generated.
! Box occupancy is only computed for levels _above_ the deepest
! level since deepest level occupancy corresponds to raw moments,
! not boxed moments, which is done separately.

subroutine fmm_build_RHS_map(mm_paras,RHS_map)

  implicit none
  type(box_mm_data), intent(in)    :: mm_paras(:)
  type(fmm_local_map), intent(out) :: RHS_map(:)

  type(box_mm_paras) :: box
  integer(INTK) :: i, id, level, ndim

  do level=deepest_level-1,TOP_LEVEL,-1

    ndim = size(mm_paras(level)%RHS_paras)
    allocate(RHS_map(level)%box_list(ndim))

    ! Initialise
    do i=1,ndim
      nullify(RHS_map(level)%box_list(i)%head)
      RHS_map(level)%box_list(i)%occ = 0
    end do

    ! We drive the loop over the child box level
    do i=1,size(mm_paras(level+1)%RHS_paras)
      ! 'map_up' indexes moments but, if in sync, then also the para set
      id = mm_paras(level+1)%RHS_paras(i)%map_up
      box%box = fmm_parent_box(mm_paras(level+1)%RHS_paras(i)%box)
      box%level = level
      ! 'id' should point to RHS parent box parameters
      if (.not. fmm_same_box(mm_paras(level)%RHS_paras(id),box)) call fmm_quit('RHS paras map-up to parent broken')
      ! Append item to linked list
      call fmm_add_item(RHS_map(level)%box_list(id),i)
    end do

  end do

end subroutine fmm_build_RHS_map

!-------------------------------------------------------------------------------

subroutine fmm_add_item(box_map,raw_id)

  implicit none
  type(id_list), intent(inout) :: box_map
  integer(INTK), intent(in)    :: raw_id
  type(id_node), pointer       :: new_node

  if (box_map%occ == 0) then
    ! Start new linked list for these boxes
    box_map%occ = 1
    allocate(box_map%head)
    box_map%head%id = raw_id
    nullify(box_map%head%next)
    return
  end if

  box_map%occ = box_map%occ+1

  allocate(new_node)
  new_node%id = raw_id
  if (associated(box_map%head%next)) then
    ! More than one entry in list (including head)
    ! so point new_node to old second entry
    new_node%next => box_map%head%next
    ! Point head to new_node
    nullify(box_map%head%next)
    box_map%head%next => new_node
  else
    ! Only head so far; make new_node our second entry
    box_map%head%next => new_node
    nullify(new_node%next)   ! end of list
  end if

end subroutine fmm_add_item

!-------------------------------------------------------------------------------
! Here we generate a list for each LHS box of all
! the RHS boxes spatially local to the LHS box.
! For now 'local' includes all the parent's near-field boxes.

subroutine fmm_build_local_map(mm_paras,RHS_map)

  implicit none
  type(box_mm_data), intent(in)   :: mm_paras(:)
  type(fmm_local_map), intent(in) :: RHS_map(:)

  type(box_mm_paras) :: LHS_box, RHS_box
  integer(INTK) :: i, j, id, lev
  type(id_node), pointer :: LHS_ptr, RHS_ptr

  lev = TOP_LEVEL
  do i=1,size(mm_paras(lev)%LHS_paras)
    do j=1,size(mm_paras(lev)%RHS_paras)
      LHS_box = mm_paras(lev)%LHS_paras(i)
      RHS_box = mm_paras(lev)%RHS_paras(j)
      if (.not. fmm_RFF_boxes(LHS_box,RHS_box)) then
        call fmm_add_item(map_at_level(lev)%box_list(i),j)
      end if
    end do
  end do

  do lev=TOP_LEVEL+1,deepest_level
    do i=1,size(mm_paras(lev)%LHS_paras)

      ! 'map_up' indexes moments but, if in sync, then also the para set
      id = mm_paras(lev)%LHS_paras(i)%map_up
      LHS_box%box = fmm_parent_box(mm_paras(lev)%LHS_paras(i)%box)
      LHS_box%level = lev-1
      ! 'id' should point to LHS parent box parameters
      if (.not. fmm_same_box(mm_paras(lev-1)%LHS_paras(id),LHS_box)) call fmm_quit('LHS paras map-up to parent broken')

      if (map_at_level(lev-1)%box_list(id)%occ == 0) cycle
      LHS_ptr => map_at_level(lev-1)%box_list(id)%head
      ! Loop over parent's near-field boxes
      parent: do
        if (RHS_map(lev-1)%box_list(LHS_ptr%id)%occ == 0) then
          if (.not. associated(LHS_ptr%next)) exit parent
          LHS_ptr => LHS_ptr%next
          cycle parent
        end if
        RHS_ptr => RHS_map(lev-1)%box_list(LHS_ptr%id)%head
        ! Get children and add to list of local boxes
        LHS_box = mm_paras(lev)%LHS_paras(i)
        child: do
          RHS_box = mm_paras(lev)%RHS_paras(RHS_ptr%id)
          if (.not. fmm_RFF_boxes(LHS_box,RHS_box)) then
            call fmm_add_item(map_at_level(lev)%box_list(i),RHS_ptr%id)
          end if
          if (.not. associated(RHS_ptr%next)) exit child
          RHS_ptr => RHS_ptr%next
        end do child

        if (.not. associated(LHS_ptr%next)) exit parent
        LHS_ptr => LHS_ptr%next
      end do parent

    end do
  end do

end subroutine fmm_build_local_map

!-------------------------------------------------------------------------------

subroutine fmm_get_local_paras(id,RHS_all,pair_type,RHS_local,ndim)

  implicit none
  type(gen_mm_paras), intent(out) :: RHS_local
  type(gen_mm_paras), intent(in)  :: RHS_all
  integer(INTK), intent(out)      :: ndim
  integer(INTK), intent(in)       :: pair_type
  integer(INTK), intent(in)       :: id

  type(id_node), pointer :: ptr
  integer(INTK) :: level, i

  select case (pair_type)
    case (LHS_raw_RHS_raw)
      call fmm_quit('local_paras: raw_raw NYI')

    case (LHS_box_RHS_box)

      if (associated(RHS_local%box_paras)) call fmm_quit('RHS_local')
      level = RHS_all%box_paras(1)%level
      ndim = map_at_level(level)%box_list(id)%occ
      if (ndim == 0) return
      allocate(RHS_local%box_paras(ndim))

      ptr => map_at_level(level)%box_list(id)%head
      i = 0
      do ! over linked-list of local boxes
        i = i+1
        RHS_local%box_paras(i) = RHS_all%box_paras(ptr%id)
        if (.not. associated(ptr%next)) exit
        ptr => ptr%next
      end do

    case default
      call fmm_quit('local_paras: requested T_pair type!')
  end select

end subroutine fmm_get_local_paras

!-------------------------------------------------------------------------------

end module fmm_local_search
