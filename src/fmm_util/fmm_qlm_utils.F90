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

module fmm_qlm_utils

use fmm_global_paras, only: INTK, REALK, raw_mm_paras, raw_mm_data, id_list, id_node, One, Two, Half
use fmm_stats, only: stat_pkd_moms_LHS, stat_pkd_moms_rHS, stat_screened_moms_RHS
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_renormalise_qlm, &
          fmm_sort_paras_wrt_centre, &
          fmm_assign_batches, &
          fmm_factor_in_dens, &
          fmm_get_T_sym_qlm, &
          fmm_pack_raw_moments, &
          fmm_pack_raw_parameters

contains

!-------------------------------------------------------------------------------

subroutine fmm_renormalise_qlm(LMAX,qlm)

  implicit none
  integer(INTK), intent(in)  :: LMAX
  real(REALK), intent(inout) :: qlm(:,:)

  integer(INTK) :: i, L, M, p, pp
  real(REALK) :: pref

  ! prefactor to symmetrize T-matrix
  do i=1,size(qlm,2)
    do L=0,LMAX
      pp = L*(L+1)+1
      do M=-L,-1
        pref = -one/(sqrt(two*FACTORIAL(l-m)*FACTORIAL(l+m)))
        p = pp+M
        qlm(p,i) = pref*qlm(p,i)
      end do
      pref = one/FACTORIAL(L)
      p = pp  ! M=0
      qlm(p,i) = pref*qlm(p,i)
      do M=1,L
        pref = ((-1)**m)/sqrt(two*FACTORIAL(l-m)*FACTORIAL(l+m))
        p = pp+M
        qlm(p,i) = pref*qlm(p,i)
      end do
    end do
  end do

contains

real(REALK) function FACTORIAL(n)
  implicit none
  integer(INTK), intent(in) :: n
  integer(INTK) :: i
  FACTORIAL = 1
  do i=n,2,-1
    FACTORIAL = FACTORIAL*i
  end do
end function FACTORIAL

end subroutine fmm_renormalise_qlm

!-------------------------------------------------------------------------------

recursive subroutine fmm_sort_paras_wrt_centre(xyz,paras)

  use fmm_sort_paras, only: fmm_quicksort_wrt_vector

  implicit none
  integer(INTK), intent(in)         :: xyz
  type(raw_mm_paras), intent(inout) :: paras(:)

  integer(INTK) :: i, lo, hi
  real(REALK) :: q1, q2

  if (size(paras) == 1) return

  ! sort only if needed
  q1 = paras(1)%cntr(xyz)
  do i=2,size(paras)
    q2 = paras(i)%cntr(xyz)
    if (q2 < q1) then
      call fmm_quicksort_wrt_vector(paras,xyz)
      exit
    end if
    q1 = q2
  end do

  ! sub-sort next cartesian component
  lo = 1
  do i=2,size(paras)
    q1 = paras(i-1)%cntr(xyz)
    q2 = paras(i)%cntr(xyz)
    if (q2 /= q1) then
      hi = i-1
      if (xyz == 3) then
        return
      else
        call fmm_sort_paras_wrt_centre(xyz+1,paras(lo:hi))
      end if
      lo = i
    end if
  end do

  ! do last batch
  hi = size(paras)
  if (xyz == 3) then
    return
  else
    call fmm_sort_paras_wrt_centre(xyz+1,paras(lo:hi))
  end if

end subroutine fmm_sort_paras_wrt_centre

!-------------------------------------------------------------------------------
! Identify batches of moments with same centre assuming sorted order and
! assign batch index to parameter list

subroutine fmm_assign_batches(paras)

  implicit none
  type(raw_mm_paras), intent(inout) :: paras(:)

  integer(INTK) :: i, batch
  !FIXME
  real(REALK), parameter :: TOLERANCE = 1.0e-20_REALK

  batch = 1
  paras(1)%batch = batch
  do i=2,size(paras)
    if (paras(i)%cntr(3)-paras(i-1)%cntr(3) > TOLERANCE) then
      batch = batch+1
    else if (paras(i)%cntr(2)-paras(i-1)%cntr(2) > TOLERANCE) then
      batch = batch+1
    else if (paras(i)%cntr(1)-paras(i-1)%cntr(1) > TOLERANCE) then
      batch = batch+1
    end if
    paras(i)%batch = batch
  end do

end subroutine fmm_assign_batches

!-------------------------------------------------------------------------------

subroutine fmm_factor_in_dens(dens,qlm)

  implicit none
  real(REALK), intent(in)    :: dens(:)
  real(REALK), intent(inout) :: qlm(:,:)

  integer(INTK) :: i

  do i=1,size(qlm,2)
    qlm(:,i) = qlm(:,i)*dens(i)
  end do

end subroutine fmm_factor_in_dens

!-------------------------------------------------------------------------------
! Prefactorising moments to symmetrize modified T-matrix

subroutine fmm_get_T_sym_qlm(LMAX,qlm_in,qlm_out)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: qlm_in(:,:)
  real(REALK), intent(out)  :: qlm_out(:,:)

  integer(INTK) :: i, L, u, hi, lo
  real(REALK) :: pref

  do i=1,size(qlm_in,2)
    do L=0,LMAX
      u = L*(L+1)+1     ! m=0
      hi = u+L
      lo = u-L
      pref = two*((-1)**L)
      qlm_out(lo:hi,i) = pref*qlm_in(lo:hi,i)
      qlm_out(u,i) = half*pref*qlm_in(u,i)
    end do
  end do

end subroutine fmm_get_T_sym_qlm

!-------------------------------------------------------------------------------
! Get number of unique batches;
! Also check raw_paras are sorted by batch;
! Do this by checking the batch ID is always increasing

subroutine get_nbatch(paras,nbatch)

  implicit none
  type(raw_mm_paras), intent(in) :: paras(:)
  integer(INTK), intent(out)     :: nbatch
  integer(INTK) :: i, ndim

  ndim = size(paras)
  nbatch = 1
  do i=2,ndim
    if (paras(i)%batch < paras(i-1)%batch) then
      call fmm_quit('batches of packed moments not sorted!')
    end if
    ! some batch numbers are skipped so need to take care of this
    if (paras(i)%batch /= paras(i-1)%batch) then
      nbatch = nbatch+1
    end if
  end do

end subroutine get_nbatch

!-------------------------------------------------------------------------------

subroutine get_pkd_data(add_dens,raw_mms,pkd_paras,pkd_qlm)

  implicit none
  logical, intent(in)             :: add_dens
  type(raw_mm_data), intent(in)   :: raw_mms
  type(raw_mm_paras), intent(out) :: pkd_paras(:)
  real(REALK), intent(out)        :: pkd_qlm(:,:)
  integer(INTK) :: i, j, k, last_batch

  j = 0
  last_batch = -1
  do i=1,size(raw_mms%paras)
    ! Use mapping of raw parameters to raw moments
    k = raw_mms%paras(i)%id
    if (raw_mms%paras(i)%batch == last_batch) then
      ! Element same batch as previous
      if (add_dens) then
        pkd_qlm(:,j) = pkd_qlm(:,j)+raw_mms%qlm(:,k)*raw_mms%dens(k)
      else
        pkd_qlm(:,j) = pkd_qlm(:,j)+raw_mms%qlm(:,k)
      end if
    else
      ! Element in new batch
      j = j+1
      pkd_paras(j) = raw_mms%paras(i)
      if (add_dens) then
        pkd_qlm(:,j) = raw_mms%qlm(:,k)*raw_mms%dens(k)
      else
        pkd_qlm(:,j) = raw_mms%qlm(:,k)
      end if
    end if
    last_batch = raw_mms%paras(i)%batch
  end do

end subroutine get_pkd_data

!-------------------------------------------------------------------------------

subroutine get_screened_nmom(screen_thr,qlm,skip,nskip)

  implicit none
  real(REALK), intent(in)    :: screen_thr
  real(REALK), intent(in)    :: qlm(:,:)
  logical, intent(out)       :: skip(:)
  integer(INTK), intent(out) :: nskip
  integer(INTK) :: i, j

  nskip = 0
  batches: do i=1,size(qlm,2)
    skip(i) = .true.
    nskip = nskip+1
    lm_loop: do j=1,size(qlm,1)
      if (abs(qlm(j,i)) > screen_thr) then
        skip(i) = .false.
        nskip = nskip-1
        exit lm_loop
      end if
    end do lm_loop
  end do batches

end subroutine get_screened_nmom

!-------------------------------------------------------------------------------
! Here we squeeze all the significant batches of moments to be sequential
! at the top of the array, with all the insignificant moments overwritten;
! we use skip(:) as a logical mask to direct the skipping.

subroutine squeeze_significant_batches(skip,paras,qlm)

  implicit none
  logical, intent(in)               :: skip(:)
  type(raw_mm_paras), intent(inout) :: paras(:)
  real(REALK), intent(inout)        :: qlm(:,:)
  integer(INTK) :: i, j

  if (size(paras) /= size(qlm,2)) call fmm_quit('paras and qlm should be same size!')
  if (size(paras) /= size(skip)) call fmm_quit('paras and skip should be same size!')

  j = 0
  do i=1,size(paras)
    if (skip(i)) cycle
    j = j+1
    paras(j) = paras(i)
    qlm(:,j) = qlm(:,i)
  end do

end subroutine squeeze_significant_batches

!-------------------------------------------------------------------------------
! Routine to pack a set of raw moments by batches,
! where members of a batch share a common centre and extent,
! including density factoring and screening if requested.
! Note that no record is kept of the packing for later "unpacking".

subroutine fmm_pack_raw_moments(mm_data,dens,dens_thr)

  implicit none
  type(raw_mm_data), intent(inout) :: mm_data
  logical, intent(in)              :: dens
  real(REALK), intent(in)          :: dens_thr

  type(raw_mm_paras), allocatable :: pkd_paras(:)
  real(REALK), allocatable :: pkd_qlm(:,:)
  logical, allocatable :: skip(:)
  integer(INTK) :: nbatch, nskip, foo, ndim

  ! Get number of unique batches
  call get_nbatch(mm_data%paras,nbatch)

  ! Now build packed data by summing raw data in the same batch
  allocate(pkd_paras(nbatch))
  foo = size(mm_data%qlm,1)
  allocate(pkd_qlm(foo,nbatch))
  call get_pkd_data(dens,mm_data,pkd_paras,pkd_qlm)
  ndim = nbatch

  if (dens) then
    ! Perform density-based screening;
    ! Note we do not reallocate whole array, but push all the
    ! significant terms to the top, and only this array section
    ! of significant moments is then pointed to.
    allocate(skip(nbatch))
    call get_screened_nmom(dens_thr,pkd_qlm,skip,nskip)
    ndim = nbatch-nskip  ! number of significant moments
    call squeeze_significant_batches(skip,pkd_paras,pkd_qlm)
    deallocate(skip)
  end if

  ! This is just for statistics
  stat_pkd_moms_RHS = nbatch
  stat_screened_moms_RHS = ndim

  ! Reallocate the original moments and copy across the packed ones
  deallocate(mm_data%paras,mm_data%qlm)
  nullify(mm_data%paras,mm_data%qlm)
  allocate(mm_data%paras(ndim),mm_data%qlm(foo,ndim))
  mm_data%paras(:) = pkd_paras(:ndim)
  mm_data%qlm(:,:) = pkd_qlm(:,:ndim)

  deallocate(pkd_paras,pkd_qlm)

end subroutine fmm_pack_raw_moments

!-------------------------------------------------------------------------------
! Routine to drive the packing of a set of raw mm parameters by batches,
! where members of a batch share a common centre and extent,
! including the build of a mapping between the packed paras and the
! original raw paras.  This map can be used for later "unpacking".

subroutine fmm_pack_raw_parameters(mm_data)

  implicit none
  type(raw_mm_data), intent(inout) :: mm_data

  type(raw_mm_paras), allocatable :: pkd_paras(:)
  integer(INTK) :: i, j, nbatch, last_batch

  ! Get number of unique batches;
  call get_nbatch(mm_data%paras,nbatch)
  stat_pkd_moms_LHS = nbatch

  ! Initialise packed paras and batch map
  allocate(pkd_paras(nbatch))
  allocate(mm_data%batch_map(nbatch))
  do i=1,nbatch
    mm_data%batch_map(i)%occ = 0
    nullify(mm_data%batch_map(i)%head)
  end do

  ! Now build packed paras by compressing raw paras in same batch
  j = 0
  last_batch = -1
  do i=1,size(mm_data%paras)
    if (mm_data%paras(i)%batch == last_batch) then
      ! Add raw parameter mapping to existing linked-list for this batch
      call add_batch_item(mm_data%batch_map(j),mm_data%paras(i)%id)
    else
      ! Element in new batch
      j = j+1
      pkd_paras(j) = mm_data%paras(i)
      ! Linked-list for this batch is empty, so start one
      mm_data%batch_map(j)%occ = 1
      allocate(mm_data%batch_map(j)%head)
      ! Maintain mapping between parameters and moments
      mm_data%batch_map(j)%head%id = mm_data%paras(i)%id
      nullify(mm_data%batch_map(j)%head%next) ! rest of list is empty
    end if
    last_batch = mm_data%paras(i)%batch
  end do

  ! Reallocate the original parameters and copy across the packed ones
  deallocate(mm_data%paras)
  nullify(mm_data%paras)
  allocate(mm_data%paras(nbatch))
  mm_data%paras(:) = pkd_paras(:)

  deallocate(pkd_paras)

contains

subroutine add_batch_item(batch_list,raw_id)

  implicit none
  type(id_list), intent(inout) :: batch_list
  integer(INTK), intent(in)   :: raw_id
  type(id_node), pointer :: new_node

  batch_list%occ = batch_list%occ+1
  allocate(new_node)
  new_node%id = raw_id
  if (associated(batch_list%head%next)) then
    ! More than one entry in list (including head)
    ! so point new_node to old second entry
    new_node%next => batch_list%head%next
    ! Point head to new_node
    nullify(batch_list%head%next)
    batch_list%head%next => new_node
  else
    ! Only head so far; make new_node our second entry
    batch_list%head%next => new_node
    nullify(new_node%next)   ! end of list
  end if

end subroutine add_batch_item

end subroutine fmm_pack_raw_parameters

!-------------------------------------------------------------------------------

end module fmm_qlm_utils
