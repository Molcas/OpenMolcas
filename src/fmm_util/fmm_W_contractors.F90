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
! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs.

!FIXME:  ************************************************************
!FIXME:  *                                                          *
!FIXME:  *  Array bounds when W-vector is ZERO !!!!!!!!!!           *
!FIXME:  *                                                          *
!FIXME:  ************************************************************

module fmm_W_contractors

use fmm_global_paras, only: INTK, REALK, T_pair_single, T_pair_list, W_CONTRACTOR_DIRECT, W_CONTRACTOR_X, W_CONTRACTOR_FAST, &
                            W_CONTRACTOR_BOUNDARY, ZERO_VECT_TOL, DISTINCT_T_TOL, Zero, One
use fmm_stats, only: stat_W_mat_builds
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_W_contractors, &
          fmm_free_W_contractors, &
          fmm_select_W_con, &
          fmm_set_W_con_ptrs

! Public variable to stop the resetting of W_con pointers with open T-buffer
public :: fmm_lock_W_con

real(REALK), allocatable, save :: W_matrix(:,:)
integer(INTK), save :: WLDA
! Pointers to actual moments and potentials described elsewhere
real(REALK), pointer, save :: old_ptr(:,:)
real(REALK), pointer, save :: new_ptr(:,:)
! Diagnostic variables
character(len=11), save :: W_con_stat
logical, save :: fmm_lock_W_con

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_W_contractors(LMAX)

  implicit none
  integer(INTK), intent(in) :: LMAX

  if (allocated(W_matrix)) call fmm_quit('W_matrix not deallocated!')
  allocate(W_matrix((1+LMAX)**2,(1+LMAX)**2))
  WLDA = (1+LMAX)**2
  W_matrix = zero

end subroutine fmm_init_W_contractors

!-------------------------------------------------------------------------------

subroutine fmm_free_W_contractors()

  implicit none
  deallocate(W_matrix)

end subroutine fmm_free_W_contractors

!-------------------------------------------------------------------------------

subroutine fmm_set_W_con_ptrs(old,new)

  implicit none
  real(REALK), target, intent(in) :: old(:,:), new(:,:)

  if (W_con_stat /= 'initialised') call fmm_quit('no W_contractor preselected!')
  if (fmm_lock_W_con) call fmm_quit('W_buffer not empty! Cannot reset W_con!')
  nullify(old_ptr,new_ptr)
  old_ptr => old(:,:)
  new_ptr => new(:,:)

end subroutine fmm_set_W_con_ptrs

!-------------------------------------------------------------------------------

subroutine fmm_check_W_status()

  implicit none
  if ((.not. associated(old_ptr)) .or. (.not. associated(new_ptr))) then
    call fmm_quit('W_contractor pointers not associated as reqd.')
  end if

end subroutine fmm_check_W_status

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one W_pair at a time.

subroutine fmm_W_con_direct(W_pair)

  use fmm_W_worker, only: fmm_get_ltsqr_W_matrix, &
                          fmm_contract_Wq

  implicit none
  type(T_pair_single), intent(in) :: W_pair

  real(REALK) :: arr_tmp(W_pair%lm_max)
  real(REALK) :: r_mod2
  integer(INTK) :: n, m, lm_dim, hi, p, q, LMAX, JMAX
  character :: NT

  call fmm_check_W_status()

  NT = W_pair%N_or_T
  p = W_pair%paras%LHS_id
  q = W_pair%paras%RHS_id
  lm_dim = W_pair%lm_max

  r_mod2 = dot_product(W_pair%r_ab,W_pair%r_ab)
  if (r_mod2 > ZERO_VECT_TOL*ZERO_VECT_TOL) then
    LMAX = W_pair%paras%LHS_LMAX
    JMAX = W_pair%paras%RHS_LMAX
    stat_W_mat_builds = stat_W_mat_builds+one
    call fmm_get_ltsqr_W_matrix(LMAX,JMAX,W_pair%r_ab,W_matrix)
    if (LMAX /= JMAX) then
      n = (1+JMAX)**2
      m = (1+LMAX)**2
      call fmm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
    else
      arr_tmp(:lm_dim) = old_ptr(:lm_dim,q)
      call DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
      new_ptr(:lm_dim,p) = new_ptr(:lm_dim,p)+arr_tmp(:lm_dim)
    end if
  else
    hi = (1+W_pair%paras%LHS_LMAX)**2
    new_ptr(:hi,p) = new_ptr(:hi,p)+old_ptr(:hi,q)
  end if

end subroutine fmm_W_con_direct

!-------------------------------------------------------------------------------

subroutine fmm_W_con_X(W_pairs)

  use fmm_W_worker, only: fmm_get_ltsqr_W_matrix

  implicit none
  type(T_pair_list), intent(in) :: W_pairs

  real(REALK) :: arr_tmp(W_pairs%lm_max)
  real(REALK) :: lastlen, r_pq(3), r_pq_mod(3)
  integer(INTK) :: LMAX, i, p, q, lm_dim, hi
  character :: NT

  call fmm_check_W_status()

  LMAX = W_pairs%LMAX
  NT = W_pairs%N_or_T
  r_pq_mod(:) = W_pairs%r_ab(:)
  lm_dim = W_pairs%lm_max
  lastlen = zero

  do i=1,size(W_pairs%paras)

    p = W_pairs%paras(i)%LHS_id
    q = W_pairs%paras(i)%RHS_id

    hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
    if (hi < size(arr_tmp)) arr_tmp = zero
    arr_tmp(:hi) = old_ptr(:hi,q)

    if (abs(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) then
      ! zero translation, so just copy over old to new
      hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
      new_ptr(:hi,p) = new_ptr(:hi,p)+arr_tmp(:hi)
      lastlen = W_pairs%paras(i)%ratio
      cycle
    end if
    ! get matrix if different
    if (abs(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) then
      r_pq = r_pq_mod*W_pairs%paras(i)%ratio
      lastlen = W_pairs%paras(i)%ratio
      stat_W_mat_builds = stat_W_mat_builds+one
      call fmm_get_ltsqr_W_matrix(LMAX,LMAX,r_pq,W_matrix)
    end if
    ! now translate
    hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
    if (hi < lm_dim) arr_tmp = zero
    arr_tmp(:hi) = old_ptr(:hi,q)
    call DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
    hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
    new_ptr(:hi,p) = new_ptr(:hi,p)+arr_tmp(:hi)

  end do

end subroutine fmm_W_con_X

!-------------------------------------------------------------------------------

subroutine fmm_W_con_FAST(W_pairs)

  use fmm_W_worker, only: fmm_get_ltsqr_W_matrix, fmm_contract_Wq

  implicit none
  type(T_pair_list), intent(in) :: W_pairs

  real(REALK) :: arr_tmp(W_pairs%lm_max)
  real(REALK) :: lastlen, r_pq(3), r_pq_mod(3)
  integer(INTK) :: LMAX, JMAX, n, m, lm_dim, i, p, q, hi
  character :: NT

  call fmm_check_W_status()

  NT = W_pairs%N_or_T
  r_pq_mod(:) = W_pairs%r_ab(:)
  lastlen = zero
  LMAX = W_pairs%LHS_LMAX
  JMAX = W_pairs%RHS_LMAX
  n = (1+JMAX)**2
  m = (1+LMAX)**2
  lm_dim = W_pairs%lm_max

  do i=1,size(W_pairs%paras)

    p = W_pairs%paras(i)%LHS_id
    q = W_pairs%paras(i)%RHS_id
    if (abs(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) then
      ! zero translation, so just copy over old to new
      hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
      new_ptr(:hi,p) = new_ptr(:hi,p)+old_ptr(:hi,q)
      lastlen = W_pairs%paras(i)%ratio
      cycle
    end if
    if (abs(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) then
      r_pq = r_pq_mod*W_pairs%paras(i)%ratio
      lastlen = W_pairs%paras(i)%ratio
      stat_W_mat_builds = stat_W_mat_builds+one
      call fmm_get_ltsqr_W_matrix(LMAX,JMAX,r_pq,W_matrix)
    end if
    if (LMAX /= JMAX) then
      call fmm_contract_Wq(NT,W_matrix,WLDA,old_ptr(:,q),n,new_ptr(:,p),m)
    else
      hi = (1+W_pairs%paras(i)%RHS_LMAX)**2
      if (hi < lm_dim) arr_tmp = zero
      arr_tmp(:hi) = old_ptr(:hi,q)
      call DTRMV('L',NT,'U',lm_dim,W_matrix,WLDA,arr_tmp,1)
      hi = (1+W_pairs%paras(i)%LHS_LMAX)**2
      new_ptr(:hi,p) = new_ptr(:hi,p)+arr_tmp(:hi)
    end if

  end do

end subroutine fmm_W_con_FAST

!-------------------------------------------------------------------------------

subroutine fmm_W_con_BOUNDARY(W_pairs)

  use fmm_W_worker, only: fmm_get_boundary_W_matrix

  implicit none
  type(T_pair_list), intent(in) :: W_pairs

  real(REALK) :: lastlen, r_pq(3), r_pq_mod(3)
  integer(INTK) :: LMAX, n, i, p, q

  call fmm_check_W_status()

  r_pq_mod(:) = W_pairs%r_ab(:)
  lastlen = zero
  LMAX = W_pairs%RHS_LMAX
  n = (1+LMAX)**2

  do i=1,size(W_pairs%paras)

    p = W_pairs%paras(i)%LHS_id
    q = W_pairs%paras(i)%RHS_id
    if (abs(W_pairs%paras(i)%ratio) < ZERO_VECT_TOL) then
      ! zero translation, so just copy over old to new
      new_ptr(1,p) = new_ptr(1,p)+old_ptr(1,q)
      lastlen = W_pairs%paras(i)%ratio
      cycle
    end if
    if (abs(W_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) then
      r_pq = r_pq_mod*W_pairs%paras(i)%ratio
      lastlen = W_pairs%paras(i)%ratio
      stat_W_mat_builds = stat_W_mat_builds+one
      call fmm_get_boundary_W_matrix(LMAX,r_pq,W_matrix)
    end if

    ! Always perform transpose contraction to get boundary potential
    new_ptr(1,p) = new_ptr(1,p)+dot_product(W_matrix(1:n,1),old_ptr(:,q))

  end do

end subroutine fmm_W_con_BOUNDARY

!-------------------------------------------------------------------------------

subroutine fmm_select_W_con(W_con_ID)

  implicit none
  integer(INTK), intent(in) :: W_con_ID
  external :: fmm_store_w_contractor

  if (.not. allocated(W_matrix)) call fmm_quit('W_matrix not allocated!')

  select case (W_con_ID)
    case (W_CONTRACTOR_DIRECT)
      call fmm_store_w_contractor(fmm_W_con_DIRECT)
    case (W_CONTRACTOR_X)
      call fmm_store_w_contractor(fmm_W_con_X)
    case (W_CONTRACTOR_FAST)
      call fmm_store_w_contractor(fmm_W_con_FAST)
    case (W_CONTRACTOR_BOUNDARY)
      call fmm_store_w_contractor(fmm_W_con_BOUNDARY)
    case default
      call fmm_quit('invalid W_contractor requested!')
  end select
  ! initialise diagnostics
  W_con_stat = 'initialised'
  W_con_stat = 'initialised'
  fmm_lock_W_con = .false.

end subroutine fmm_select_W_con

!-------------------------------------------------------------------------------

end module fmm_W_contractors
