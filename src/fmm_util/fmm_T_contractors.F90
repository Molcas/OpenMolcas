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

module fmm_T_contractors

! Module containing routines to generate contributions to a
! multipole potential on the fly.
! Note that another module allows contraction of T-pairs to generate
! an energy or J-matrix directly (e.g.to exploit symmetry).

use fmm_global_paras, only: INTK, REALK, scheme_paras, T_pair_single, T_pair_list, T_pair_batch, T_CONTRACTOR_MULTI,  &
                            T_CONTRACTOR_BOUNDARY, T_CONTRACTOR_DIRECT, T_CONTRACTOR_TREE, T_CONTRACTOR_SCALE_TREE,  &
                            T_CONTRACTOR_SCALE, T_CONTRACTOR_FULL, NEAR_FIELD, DISTINCT_T_TOL, TMATM_DF, Zero, One, Half
use fmm_stats, only: stat_T_mat_builds, fmm_init_matrix_stats
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_T_contractors, &
          fmm_free_T_contractors, &
          fmm_select_T_con, &
          fmm_set_T_con_ptrs

! Public variable to stop the resetting of T_con pointers with open T-buffer
public :: fmm_lock_T_con

real(REALK), allocatable, save :: T_matrix(:,:)
integer(INTK), save :: TLDA
! when building multiple T matrices together
real(REALK), allocatable, save :: T_mats(:,:,:)

! Pointers to actual moments and potentials described elsewhere
real(REALK), pointer, save :: Vff_ptr(:,:)
real(REALK), pointer, save :: qlm_ptr(:,:)

! Diagnostic variables
character(len=11), save :: T_con_stat
logical, save :: fmm_lock_T_con

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_T_contractors(scheme)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK) :: LMAX, T_con
  integer(INTK) :: lm_max

  LMAX = scheme%trans_lmax
  lm_max = (LMAX+1)**2
  if (scheme%phase == NEAR_FIELD) then
    T_con = scheme%T_con%NF_id
  else
    T_con = scheme%T_con%FF_id
  end if

  select case(T_con)
    case(T_CONTRACTOR_MULTI)
      if (allocated(T_mats)) call fmm_quit('T_mats not deallocated!')
      !FIXME: will need to change TMATM_DF if changed elsewhere
      allocate(T_mats(TMATM_DF,lm_max,lm_max))
      T_mats = zero
    case(T_CONTRACTOR_BOUNDARY)
      if (allocated(T_matrix)) call fmm_quit('T_matrix not deallocated!')
      allocate(T_matrix(lm_max,1))
      T_matrix = zero
    case default
      if (allocated(T_matrix)) call fmm_quit('T_matrix not deallocated!')
      allocate(T_matrix(lm_max,lm_max))
      T_matrix = zero
  end select
  TLDA = lm_max

  ! code for statistics only
  call fmm_init_matrix_stats('T')

end subroutine fmm_init_T_contractors

!-------------------------------------------------------------------------------

subroutine fmm_free_T_contractors()

  implicit none

  if (allocated(T_matrix)) deallocate(T_matrix)
  if (allocated(T_mats)) deallocate(T_mats)

end subroutine fmm_free_T_contractors

!-------------------------------------------------------------------------------

subroutine fmm_set_T_con_ptrs(Vff,qlm)

  implicit none
  real(REALK), target, intent(in) :: Vff(:,:), qlm(:,:)

  if (T_con_stat /= 'initialised') call fmm_quit('no T_contractor preselected!')
  if (fmm_lock_T_con) call fmm_quit('T_buffer not empty! Cannot reset T_con!')
  nullify(Vff_ptr,qlm_ptr)
  Vff_ptr => Vff(:,:)
  qlm_ptr => qlm(:,:)

end subroutine fmm_set_T_con_ptrs

!-------------------------------------------------------------------------------
! Built for DIRECT scheme.  Can only take one T_pair at a time.
! Only generates a square T-matrix with dim = MAX(LHS_LMAX,RHS_LMAX)
! Generates local potential to order MAX(LHS_LMAX,RHS_LMAX)
! Adds to Vff only up to LHS_LMAX.
! FIXME: because T_matrix is only built for (l+j)<=LMAX then with dynamic
! contraction we need LEXTRA>0 even for primitive basis sets and FQUAD.

subroutine fmm_T_con_DIRECT(T_pair)

  use fmm_T_worker, only: fmm_get_SPLTSQ_T_matrix, &
                          fmm_contract_Tq

  implicit none
  type(T_pair_single), intent(in) :: T_pair

  real(REALK) :: Vff_tmp(T_pair%lm_max)
  integer(INTK) :: weight, iLHS, iRHS, hi

  stat_T_mat_builds = stat_T_mat_builds+one
  call fmm_get_SPLTSQ_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

  iRHS = T_pair%paras%RHS_id
  hi = T_pair%lm_max

  Vff_tmp = fmm_contract_Tq(T_pair%LMAX,qlm_ptr(:hi,iRHS),T_matrix)
  iLHS = T_pair%paras%LHS_id
  hi = (1+T_pair%paras%LHS_LMAX)**2
  weight = T_pair%paras%weight
  !write(LUPRI,*) iRHS,qlm_ptr(1,iRHS),weight,Vff_tmp(:hi)
  Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)+weight*Vff_tmp(:hi)

end subroutine fmm_T_con_DIRECT

!-------------------------------------------------------------------------------
! Modification of T_con_DIRECT for evaluation of boundary potential
! when only the first element of Vff needs to be built

subroutine fmm_T_con_BOUNDARY(T_pair)

  use fmm_T_worker, only: fmm_get_boundary_T_matrix

  implicit none
  type(T_pair_single), intent(in) :: T_pair

  real(REALK) :: Vff_tmp
  integer(INTK) :: weight, iLHS, iRHS, hi

  stat_T_mat_builds = stat_T_mat_builds+one
  call fmm_get_boundary_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

  iRHS = T_pair%paras%RHS_id
  hi = T_pair%lm_max

  Vff_tmp = half*dot_product(qlm_ptr(:hi,iRHS),T_matrix(:hi,1))

  iLHS = T_pair%paras%LHS_id
  hi = (1+T_pair%paras%LHS_LMAX)**2
  weight = T_pair%paras%weight
  !write(LUPRI,*) iRHS,qlm_ptr(1,iRHS),weight,Vff_tmp

  Vff_ptr(1,iLHS) = Vff_ptr(1,iLHS)+weight*Vff_tmp

end subroutine fmm_T_con_BOUNDARY

!-------------------------------------------------------------------------------

subroutine fmm_T_con_TREE(T_pairs)

  use fmm_T_worker, only: fmm_get_SPLTSQ_T_matrix, &
                          fmm_contract_Tq

  implicit none
  type(T_pair_list), intent(in) :: T_pairs

  real(REALK) :: Vff_tmp(T_pairs%lm_max)
  real(REALK) :: lastlen, r_pq(3), r_pq_mod(3)
  integer(INTK) :: LMAX, weight, i, p, q, hi

  r_pq_mod(:) = T_pairs%r_ab(:)
  LMAX = T_pairs%LMAX
  lastlen = zero

  do i=1,size(T_pairs%paras)

    p = T_pairs%paras(i)%LHS_id
    q = T_pairs%paras(i)%RHS_id

    if (abs(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) then
      stat_T_mat_builds = stat_T_mat_builds+one
      r_pq = r_pq_mod*T_pairs%paras(i)%ratio
      call fmm_get_SPLTSQ_T_matrix(LMAX,r_pq,T_matrix)
      lastlen = T_pairs%paras(i)%ratio
    end if

    hi = T_pairs%lm_max
    Vff_tmp = fmm_contract_Tq(LMAX,qlm_ptr(:hi,q),T_matrix)
    hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
    weight = T_pairs%paras(i)%weight
    Vff_ptr(:hi,p) = Vff_ptr(:hi,p)+weight*Vff_tmp(:hi)

  end do

end subroutine fmm_T_con_TREE

!-------------------------------------------------------------------------------

subroutine fmm_T_con_SCALE(T_pairs)

  use fmm_T_worker, only: fmm_get_SPLTSQ_T_matrix, &
                          fmm_contract_Tq, &
                          fmm_scale_vec

  implicit none
  type(T_pair_batch), intent(in) :: T_pairs

  integer(INTK) :: LMAX, i, p, q, hi, lastq
  real(REALK) :: ratio, lastlen, pref
  real(REALK) :: Vff_tmp(T_pairs%items(1)%lm_max)
  real(REALK) :: scaled_qlm(T_pairs%items(1)%lm_max)
  real(REALK) :: scale_vec(T_pairs%items(1)%lm_max)
  logical :: new_vec

  pref = one
  lastq = -1
  lastlen = one
  scale_vec(:) = one
  LMAX = T_pairs%items(1)%LMAX
  stat_T_mat_builds = stat_T_mat_builds+one
  call fmm_get_SPLTSQ_T_matrix(LMAX,T_pairs%items(1)%r_ab,T_matrix)

  iloop: do i=1,T_pairs%ndim

    p = T_pairs%items(i)%paras%LHS_id
    q = T_pairs%items(i)%paras%RHS_id

    new_vec = .false.
    ratio = T_pairs%items(i)%paras%ratio
    if (abs(ratio-lastlen) > DISTINCT_T_TOL) then
      call fmm_scale_vec(LMAX,ratio,scale_vec,pref)
      lastlen = ratio
      new_vec = .true.
    end if

    hi = T_pairs%items(i)%lm_max
    if (new_vec .or. (q /= lastq)) then
      scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
      lastq = q
    end if

    Vff_tmp = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

    hi = (1+T_pairs%items(i)%paras%LHS_LMAX)**2
    Vff_ptr(:hi,p) = Vff_ptr(:hi,p)+pref*scale_vec(:hi)*Vff_tmp(:hi)

  end do iloop

end subroutine fmm_T_con_SCALE

!-------------------------------------------------------------------------------

subroutine fmm_T_con_SCALE_TREE(T_pairs)

  use fmm_T_worker, only: fmm_get_SPLTSQ_T_matrix, &
                          fmm_contract_Tq, &
                          fmm_scale_vec

  implicit none
  type(T_pair_list), intent(in) :: T_pairs

  integer(INTK) :: LMAX, i, p, q, hi, lastq
  real(REALK) :: weight, lastlen, pref
  real(REALK) :: Vff_tmp(T_pairs%lm_max)
  real(REALK) :: scaled_qlm(T_pairs%lm_max)
  real(REALK) :: scale_vec(T_pairs%lm_max)
  logical :: new_vec

  pref = one
  lastq = -1
  lastlen = one
  scale_vec(:) = one
  LMAX = T_pairs%LMAX
  stat_T_mat_builds = stat_T_mat_builds+one
  call fmm_get_SPLTSQ_T_matrix(LMAX,T_pairs%r_ab,T_matrix)

  iloop: do i=1,size(T_pairs%paras)

    p = T_pairs%paras(i)%LHS_id
    q = T_pairs%paras(i)%RHS_id

    new_vec = .false.
    if (abs(T_pairs%paras(i)%ratio-lastlen) > DISTINCT_T_TOL) then
      call fmm_scale_vec(LMAX,T_pairs%paras(i)%ratio,scale_vec,pref)
      lastlen = T_pairs%paras(i)%ratio
      new_vec = .true.
    end if

    hi = T_pairs%lm_max
    if (new_vec .or. (q /= lastq)) then
      scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,q)
      lastq = q
    end if

    Vff_tmp = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_matrix(:hi,:hi))

    hi = (1+T_pairs%paras(i)%LHS_LMAX)**2
    weight = pref*T_pairs%paras(i)%weight
    Vff_ptr(:hi,p) = Vff_ptr(:hi,p)+weight*scale_vec(:hi)*Vff_tmp(:hi)

  end do iloop

end subroutine fmm_T_con_SCALE_TREE

!-------------------------------------------------------------------------------
! Special contractor designed to take batch of T-pairs ordered in pairs
! such that ( b T1 1; a -T1 2; c T2 1; a -T2 3 ....)
! to halve the number of T matrix builds (using fmm_scale...)
! But only the common RHS qlm can be contracted simultaneously which slows
! it down a lot.
!
!subroutine fmm_T_con_MULTI(T_pairs)
!
!  use fmm_T_worker_multi, only: fmm_get_SPLTSQ_T_matrices,    &
!                                fmm_contract_multi_Tq
!  use fmm_T_worker, only: fmm_contract_Tq, fmm_scale_vec
!
!  implicit none
!  type(T_pair_batch), intent(in) :: T_pairs
!
!  real(REALK), allocatable :: Vff_tmp(:,:)
!
!  real(REALK), allocatable :: scaled_qlm(:)
!  real(REALK), allocatable :: scale_vec(:)
!
!  real(REALK) :: T_vectors((T_pairs%ndim/2),3)
!  real(REALK) :: pref, weight
!  integer(INTK) :: i, j, iLHS, iRHS, iRHS_last, hi, LMAX, nT
!
!  ! FIRST DO MULTIPLE T contracted with SAME RHS
!  !---------------------------------------------
!
!  if (btest(T_pairs%ndim,0)) call fmm_quit('ndim not EVEN!')
!  nT = T_pairs%ndim/2
!
!  LMAX = 0
!  iRHS_last = T_pairs%items(1)%paras%RHS_id
!  do i=1,nT
!    j = 2*i-1
!    LMAX = max(LMAX,T_pairs%items(i)%LMAX)
!    iRHS = T_pairs%items(j)%paras%RHS_id
!    ! get *distinct* T-vectors as every second item in batch
!    T_vectors(i,:) = T_pairs%items(j)%r_ab(:)
!    if (iRHS /= iRHS_last) call fmm_quit('must have same qlm on RHS')
!    iRHS_last = iRHS
!  end do
!
!  stat_T_mat_builds = stat_T_mat_builds + nT
!  call fmm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats(:nT,:,:))
!
!  hi = (1+LMAX)**2
!  allocate(Vff_tmp(nT,hi))
!  allocate(scaled_qlm(hi))
!  allocate(scale_vec(hi))
!  Vff_tmp(:,:) = zero
!  scale_vec(:) = one
!
!  iRHS = T_pairs%items(1)%paras%RHS_id
!
!  Vff_tmp(:,:hi) = fmm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),T_mats(:nT,:,:),nT)
!
!  do i=1,nT
!    j = 2*i-1
!    iLHS = T_pairs%items(j)%paras%LHS_id
!    hi = (1+T_pairs%items(j)%paras%LHS_LMAX)**2
!    weight = T_pairs%items(j)%paras%weight
!    Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)+weight*Vff_tmp(i,:hi)
!  end do
!
!  ! NOW DO remaining half of batched T-pairs corresponding to same LHS
!  !-------------------------------------------------------------------
!
!  iLHS = T_pairs%items(2)%paras%LHS_id
!  do i=1,nT
!    ! get T_matrix correpxonding to the minus T-vector
!    call fmm_scale_vec(LMAX,-one,scale_vec,pref)
!    iRHS = T_pairs%items(2*i)%paras%RHS_id
!    scaled_qlm = scale_vec(:hi)*qlm_ptr(:hi,iRHS)
!    Vff_tmp(1,:hi) = fmm_contract_Tq(LMAX,scaled_qlm(:hi),T_mats(i,:,:))
!    weight = pref
!    Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)+weight*scale_vec(:hi)*Vff_tmp(1,:hi)
!  end do
!
!  deallocate(Vff_tmp)
!  deallocate(scaled_qlm)
!  deallocate(scale_vec)
!
!end subroutine fmm_T_con_MULTI

!-------------------------------------------------------------------------------
! Contractor designed to build multiple T matrices simultaneously.
! Exact performance will depend on architecture and choice of NDIM.
! Assumes RHS moments are the same, but T-vectors are different.

subroutine fmm_T_con_MULTI(T_pairs)

  use fmm_multiple_T_worker, only: fmm_get_SPLTSQ_T_matrices, &
                                   fmm_contract_multi_Tq

  implicit none
  type(T_pair_batch), intent(in) :: T_pairs

  real(REALK) :: Vff_tmp(T_pairs%ndim,T_pairs%items(1)%lm_max)
  real(REALK) :: T_vectors(T_pairs%ndim,3)
  integer(INTK) :: i, iLHS, iRHS, hi, LMAX, nT

  nT = T_pairs%ndim
  LMAX = T_pairs%items(1)%LMAX

  do i=1,nT
    T_vectors(i,:) = T_pairs%items(i)%r_ab
    iRHS = T_pairs%items(i)%paras%RHS_id
    if (iRHS /= T_pairs%items(max(i-1,1))%paras%RHS_id) then
      call fmm_quit('RHS moments not sorted in fmm_T_con_MULTI')
    end if
  end do

  stat_T_mat_builds = stat_T_mat_builds+nT
  call fmm_get_SPLTSQ_T_matrices(nT,LMAX,T_vectors,T_mats(:nt,:,:))

  iRHS = T_pairs%items(1)%paras%RHS_id
  hi = (1+LMAX)**2

  Vff_tmp(:,:hi) = fmm_contract_multi_Tq(LMAX,qlm_ptr(:hi,iRHS),T_mats(:nT,:,:),nT)

  do i=1,nT
    iLHS = T_pairs%items(i)%paras%LHS_id
    Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)+Vff_tmp(i,:hi)
  end do

end subroutine fmm_T_con_MULTI

!-------------------------------------------------------------------------------
! Builds full T-matrix for exact contraction at low orders

subroutine fmm_T_con_FULL(T_pair)

  use fmm_T_worker, only: fmm_get_FLTSQ_T_matrix, fmm_postfac_Vff

  implicit none
  type(T_pair_single), intent(in) :: T_pair

  real(REALK) :: Vff_tmp(T_pair%lm_max)
  integer(INTK) :: iLHS, iRHS, hi

  stat_T_mat_builds = stat_T_mat_builds+one
  call fmm_get_FLTSQ_T_matrix(T_pair%LMAX,T_pair%r_ab,T_matrix)

  iRHS = T_pair%paras%RHS_id
  hi = T_pair%lm_max
  call DSYMV('L',hi,one,T_matrix,TLDA,qlm_ptr(:,iRHS),1,zero,Vff_tmp,1)

  iLHS = T_pair%paras%LHS_id
  hi = (1+T_pair%paras%LHS_LMAX)**2
  call fmm_postfac_Vff(T_pair%paras%LHS_LMAX,Vff_tmp)
  Vff_ptr(:hi,iLHS) = Vff_ptr(:hi,iLHS)+Vff_tmp(:hi)

end subroutine fmm_T_con_FULL

!-------------------------------------------------------------------------------

subroutine fmm_select_T_con(scheme)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK) :: T_con_ID
  external fmm_store_t_contractor

  if (scheme%phase == NEAR_FIELD) then
    T_con_ID = scheme%T_con%NF_id
  else
    T_con_ID = scheme%T_con%FF_id
  end if

  select case(T_con_ID)
    case(T_CONTRACTOR_DIRECT)
      call fmm_store_t_contractor(fmm_T_con_DIRECT)
    case(T_CONTRACTOR_BOUNDARY)
      call fmm_store_t_contractor(fmm_T_con_BOUNDARY)
    case(T_CONTRACTOR_TREE)
      call fmm_store_t_contractor(fmm_T_con_TREE)
    case(T_CONTRACTOR_SCALE_TREE)
      call fmm_store_t_contractor(fmm_T_con_SCALE_TREE)
    case(T_CONTRACTOR_SCALE)
      call fmm_store_t_contractor(fmm_T_con_SCALE)
    case(T_CONTRACTOR_MULTI)
      call fmm_store_t_contractor(fmm_T_con_MULTI)
    case(T_CONTRACTOR_FULL)
      call fmm_store_t_contractor(fmm_T_con_FULL)
    case default
      call fmm_quit('invalid T_contractor requested!')
  end select
  ! initialise diagnostics
  T_con_stat = 'initialised'
  fmm_lock_T_con = .false.

end subroutine fmm_select_T_con

!-------------------------------------------------------------------------------

end module fmm_T_contractors
