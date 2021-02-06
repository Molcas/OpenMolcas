!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2002, Pawel Salek                                      *
!               2002, Mark A. Watson                                   *
!***********************************************************************

module fmm_tree_buffer
! Pawel Salek
! approx. 2002.05.13
! Mark Watson : Sep. 2002.
! Generalised to allow for translation pairs too.
! i.e. can store and evaluate on the fly W-pairs
! Also included field for dynamic LMAX
!
! fmm_tree_buffer packs together interaction pairs having same R vector
! to reduce the cost of generating interactions matrices T - the matrix
! is constructed only once for entire set.
! This module needs a routine that will evaluate entire set corresponding
! to same R vector.
!
! One limitation cast on the implementation is to use possibly
! large but limited memory resources. This requires an ability to evaluate
! and expunge occasionally the interaction sets as they are created.
!
! The pack_inter_tree is used to find out unique translation vectors.
! This is done by sorting all translation vectors first into
! surfaces, then into lines within surfaces, and finally into
! unique points.

use fmm_global_paras, only: INTK, REALK, T_paras, T_pair_single, T_pair_list, DISTINCT_T_TOL, MAX_AVG_PER_NODE, &
                            SORT_BY_RHS_MMS, START_LEN, One
use fmm_stats, only: stat_tpack_chunks, stat_tpack_total, stat_tpack_unique
use fmm_utils, only: fmm_quit

! Public functions:
public :: fmm_tree_buffer_init, &
          fmm_tree_buffer_finish, &
          fmm_tree_buffer_add

! Public types:
public :: PointNode

private
integer(INTK), save :: pack_sort_order

! PointNode: private structure used
! for on-the-fly packing of the interactions.
!
type PointNode
  integer(INTK)            :: level !for debugging:1-plane, 2-line, etc
  real(REALK)              :: coord ! normalized, the vector scale==1
  type(PointNode), pointer :: left, right
  ! the one below used by non-leaves only
  type(PointNode), pointer :: this ! for level==1: this plane,
  !     level==2: this line
  !     level==3: self (for consistency)
  ! the three below used by leaves only (i.e. ones with level==3)
  type(T_paras), pointer   :: abl1l2(:)
  ! largest LMAX in each abl1l2 set stored here
  integer(INTK)            :: entries_used
  integer(INTK)            :: LHS_LMAX, RHS_LMAX, LMAX
  ! For translations must know if qlm or Vff translation
  character                :: N_or_T
end type PointNode

! PRIVATE DATA
! pack_inter_tree is a pack-on-the-fly, evaluate-when-needed structure for
! efficient, memory-adapted, R-vector packed evaluation of the
! interaction pairs.
type(PointNode), allocatable, target, save :: pack_inter_tree(:)
integer(INTK), save :: pack_inter_tree_used
integer(INTK), save :: pack_total_in_current_chunk

contains
!---------------------------------------------------------------------------
! fmm_tree_buffer_init:
! initialize data structures and counters when starting packing the tree.
! interaction sorting depends on whether the contraction scheme is
! T-matrix driven or U-matrix driven.  Here we initialise a
! module-wide logical to determine this for all schemes
!
subroutine fmm_tree_buffer_init(t_size,sort_order)
  implicit none
  integer(INTK), intent(in) :: t_size
  integer(INTK), intent(in) :: sort_order

  allocate(pack_inter_tree(t_size))
  pack_inter_tree_used = 0
  pack_sort_order = sort_order
  pack_total_in_current_chunk = 0
end subroutine fmm_tree_buffer_init

!---------------------------------------------------------------------------

subroutine fmm_tree_buffer_finish(pack_ev)
  implicit none
  external :: pack_ev
  type(PointNode), pointer :: root

  if (pack_inter_tree_used > 0) then
    root => pack_inter_tree(1)
    call fmm_tpack_process(root,pack_ev)
    pack_inter_tree_used = 0
    pack_total_in_current_chunk = 0
  end if
  if (allocated(pack_inter_tree)) deallocate(pack_inter_tree)
end subroutine fmm_tree_buffer_finish

!---------------------------------------------------------------------------

subroutine fmm_tpack_process(root,pack_ev)
  implicit none
  type(PointNode), pointer :: root
  external                 :: pack_ev

  integer(INTK) :: i

  call node_evaluator(root,1,pack_ev)

  do i=1,pack_inter_tree_used
    if (associated(pack_inter_tree(i)%abl1l2)) deallocate(pack_inter_tree(i)%abl1l2)
    nullify(pack_inter_tree(i)%abl1l2) ! not strictly neccesary
  end do
  pack_inter_tree_used = 0
  pack_total_in_current_chunk = 0
  nullify(root)
end subroutine fmm_tpack_process

!====================================================================
! point_node_new:
! NOTE: we choose for conveninence that:
! a). non-leaves have abl1l2 not associated.
! b). leaves have leaf%this => leaf
recursive subroutine point_node_new(node,level,r,ll,lr,lm,NT,T_pair)
  implicit none
  type(PointNode), pointer  :: node
  type(T_paras), intent(in) :: T_pair
  real(REALK), intent(in)   :: r(3)
  integer(INTK), intent(in) :: level
  ! put these variables in a TYPE ??
  integer(INTK), intent(in) :: ll, lr, lm
  ! For translations must know if qlm or Vff translation
  character, intent(in)     :: NT

  pack_inter_tree_used = pack_inter_tree_used+1
  node => pack_inter_tree(pack_inter_tree_used)
  node%level = level
  node%coord = r(level)
  nullify(node%left); nullify(node%right)
  if (level < 3) then
    nullify(node%abl1l2)
    call point_node_new(node%this,level+1,r,ll,lr,lm,NT,T_pair)
  else
    stat_tpack_unique = stat_tpack_unique+1
    node%this => node
    node%entries_used = 1
    node%LMAX = lm
    node%LHS_LMAX = ll
    node%RHS_LMAX = lr
    node%N_or_T = NT
    allocate(node%abl1l2(START_LEN))
    node%abl1l2(1) = T_pair
  end if
end subroutine point_node_new

!---------------------------------------------------------------------------
! node_evaluator:
! walks over the tree and evaluates the interaction pairs.
!
recursive subroutine node_evaluator(node,level,pack_ev)
  implicit none
  type(PointNode), pointer  :: node
  integer(INTK), intent(in) :: level
  external                  ::  pack_ev   ! pair evaluator

  real(REALK), save :: r(3)

  if (.not. associated(node)) return
  if (associated(node%left)) call node_evaluator(node%left,level,pack_ev)
  if (associated(node%right)) call node_evaluator(node%right,level,pack_ev)
  r(level) = node%coord
  if (level < 3) then
    call node_evaluator(node%this,level+1,pack_ev)
  else
    ! evaluate current point using provided evaluator.
    ! We sort node%ratio node%abl1l2 wrt to abl1l2(:)%RHS_id first.
    if (node%entries_used > 1) then
      call momentsort(node%abl1l2,node%entries_used)
    end if
    call fmm_interface_T_pair_out(pack_ev,r,node)
  end if

end subroutine node_evaluator

!---------------------------------------------------------------------------

subroutine momentsort(abl1l2,N)
  implicit none
  type(T_paras), intent(inout) :: abl1l2(:)
  integer(INTK), intent(in)    :: N
  type(T_paras) :: idxs
  integer(INTK) :: i

  do i=N/2,1,-1
    call downheap(abl1l2,i,N)
  end do

  ! abl1l2[1..N] is a heap now

  do i=N,1,-1
    idxs = abl1l2(i)
    abl1l2(i) = abl1l2(1)
    abl1l2(1) = idxs
    call downheap(abl1l2,1,i-1) ! restore a[1..i-1] heap
  end do
end subroutine momentsort

!---------------------------------------------------------------------------

subroutine downheap(abl1l2,lo,hi)
  implicit none
  type(T_paras), intent(inout) :: abl1l2(:)
  integer(INTK), intent(in)    :: lo, hi

  !  PRE: a[lo+1..hi] is a heap
  ! POST:  a[lo..hi]  is a heap
  type(T_paras) :: idxs
  integer(INTK) :: idx, child

  idxs = abl1l2(lo)

  idx = lo
  makeheap: do while (idx <= hi/2) ! while k has child(s)
    child = 2*idx
    ! pick larger child...
    if (child < hi) then ! exists right branch
      if (pack_sort_order == SORT_BY_RHS_MMS) then
        if ((abl1l2(child)%RHS_id < abl1l2(child+1)%RHS_id) .or. &
            ((abl1l2(child)%RHS_id == abl1l2(child+1)%RHS_id) .and. &
             (abl1l2(child)%ratio < abl1l2(child+1)%ratio))) child = child+1
      else
        if ((abl1l2(child)%ratio < abl1l2(child+1)%ratio) .or. &
            ((abs(abl1l2(child)%ratio-abl1l2(child+1)%ratio) < DISTINCT_T_TOL) .and. &
             (abl1l2(child)%RHS_id < abl1l2(child+1)%RHS_id))) child = child+1
      end if
    end if
    if (pack_sort_order == SORT_BY_RHS_MMS) then
      if (idxs%RHS_id >= abl1l2(child)%RHS_id) exit makeheap
    else
      if (idxs%ratio >= abl1l2(child)%ratio) exit makeheap
    end if
    abl1l2(idx) = abl1l2(child)
    idx = child
  end do makeheap

  abl1l2(idx) = idxs
end subroutine downheap

!---------------------------------------------------------------------------
! fmm_tree_buffer_add:
! creates packed list of interactions corresponding to the same R.
! packing is done as-you-go, without creation of the intermediate
! unpacked list. Tree structure is used for sorting.
! Appends given data to specified tree.
! direction is a normalized r_pq with positive x component.
subroutine fmm_tree_buffer_add(pack_ev,T_pair_in)
  implicit none
  type(T_pair_single), intent(in) :: T_pair_in
  external                        ::  pack_ev ! packed pair evaluator

  ! Local variables.
  type(PointNode), pointer :: node
  type(T_paras), pointer :: newblock(:)
  type(T_paras) :: new_T_pair
  integer(INTK) :: level
  real(REALK) :: r_pq(3)
  real(REALK) :: direction(3), r_pq_scale
  ! orders of multipole expanions per node
  integer(INTK) :: ll, lr, lm
  ! for translations must know if qlm or Vff translation per node
  character :: NT

  r_pq = T_pair_in%r_ab
  new_T_pair = T_pair_in%paras
  NT = T_pair_in%N_or_T
  ll = T_pair_in%paras%LHS_LMAX
  lr = T_pair_in%paras%RHS_LMAX
  lm = T_pair_in%LMAX

  node => pack_inter_tree(1) ! root is the first element of pack_inter_tree
  if (pack_inter_tree_used+3 > size(pack_inter_tree) .or. pack_total_in_current_chunk > MAX_AVG_PER_NODE*size(pack_inter_tree)) then
    stat_tpack_chunks = stat_tpack_chunks+one
    call fmm_tpack_process(node,pack_ev)
  end if

  ! do the job
  stat_tpack_total = stat_tpack_total+1
  pack_total_in_current_chunk = pack_total_in_current_chunk+1
  r_pq_scale = sqrt(sum(r_pq*r_pq))
  if (r_pq(1) < 0.0) r_pq_scale = -r_pq_scale
  direction = r_pq/r_pq_scale
  ! update ratio (r_ab is now the normalised vector)
  new_T_pair%ratio = r_pq_scale

  if (pack_inter_tree_used == 0) then
    call point_node_new(node,1_INTK,direction,ll,lr,lm,NT,new_T_pair)
    return
  end if

  ! level = 1 - searching for plane, 2 - searching for line, 3 - for point
  do level=1,3
    nodesearch: do while (abs(direction(level)-node%coord) > DISTINCT_T_TOL)
      if (direction(level) < node%coord) then
        if (.not. associated(node%left)) then
          call point_node_new(node%left,level,direction,ll,lr,lm,NT,new_T_pair)
          return
        else
          node => node%left
        end if
      else
        if (.not. associated(node%right)) then
          call point_node_new(node%right,level,direction,ll,lr,lm,NT,new_T_pair)
          return
        else
          node => node%right
        end if
      end if
    end do nodesearch
    node => node%this ! node%this points to node for leaves
  end do

  ! append to the current point. It is guaranteed to be the final point
  ! expand first data block if needed.
  if (node%entries_used == size(node%abl1l2)) then
    allocate(newblock(node%entries_used*2))
    !node%abl1l2(1:node%entries_used): index removed to avoid crash on IBM
    newblock(1:node%entries_used) = node%abl1l2
    deallocate(node%abl1l2)
    node%abl1l2 => newblock
  end if
  node%entries_used = node%entries_used+1
  if (node%N_or_T /= NT) call fmm_quit('inconsistent data in buffer node!')
  node%N_or_T = NT
  node%LHS_LMAX = max(node%LHS_LMAX,ll)
  node%RHS_LMAX = max(node%RHS_LMAX,lr)
  node%LMAX = max(node%LMAX,lm)
  node%abl1l2(node%entries_used) = new_T_pair

end subroutine fmm_tree_buffer_add

!---------------------------------------------------------------------------
! Here we transform the data of this packer (old) into the new format

subroutine fmm_interface_T_pair_out(pack_ev,r,node)

  implicit none
  type(PointNode), intent(in) :: node
  real(REALK), intent(in)     :: r(3)
  external                    :: pack_ev

  type(T_pair_list) :: T_pairs_out

  T_pairs_out%N_or_T = node%N_or_T
  T_pairs_out%r_ab = r  ! normalised vector
  T_pairs_out%paras => node%abl1l2(1:node%entries_used)

  T_pairs_out%LMAX = node%LMAX
  T_pairs_out%lm_max = (1+node%LMAX)**2
  T_pairs_out%LHS_LMAX = node%LHS_LMAX
  T_pairs_out%RHS_LMAX = node%RHS_LMAX

  call pack_ev(T_pairs_out)

end subroutine fmm_interface_T_pair_out

!---------------------------------------------------------------------------

end module fmm_tree_buffer
