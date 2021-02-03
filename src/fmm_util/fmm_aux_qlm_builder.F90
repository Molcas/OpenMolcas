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

module fmm_aux_qlm_builder

use fmm_global_paras, only: INTK, REALK, LUPRI, scheme_paras, raw_mm_data, USE_RAW_QLM, USE_T_SYM_QLM
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_get_aux_qlm

contains

!-------------------------------------------------------------------------------

subroutine fmm_get_aux_qlm(scheme,LHS_mms,RHS_mms)

  use fmm_qlm_utils, only: fmm_renormalise_qlm

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  type(raw_mm_data), intent(inout) :: LHS_mms,RHS_mms

  ! We use the SCALED solid harmonic formulation throughout
  call fmm_renormalise_qlm(scheme%raw_LMAX,LHS_mms%qlm)
  call fmm_renormalise_qlm(scheme%raw_LMAX,RHS_mms%qlm)

  ! Identify batches of unique centre/extent and sort them in
  ! preparation for packing
  call sort_centres_and_assign_batches(scheme,LHS_mms,RHS_mms)

  ! Get parameter mappings and preconditioned W and T moments
  call get_RHS_data(scheme,RHS_mms)
  call get_LHS_data(scheme,LHS_mms)

  ! All relevant data should now be held in %qlm_T and %qlm_W
  deallocate(LHS_mms%qlm,RHS_mms%qlm)
  nullify(LHS_mms%qlm,RHS_mms%qlm)

end subroutine fmm_get_aux_qlm

!-------------------------------------------------------------------------------

subroutine sort_centres_and_assign_batches(scheme,LHS,RHS)

  use fmm_qlm_utils, only: fmm_sort_paras_wrt_centre, fmm_assign_batches

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  type(raw_mm_data), intent(inout) :: LHS,RHS

  if (scheme%pack_LHS) then
    call fmm_sort_paras_wrt_centre(1_INTK,LHS%paras)
    call fmm_assign_batches(LHS%paras)
  end if
  if (scheme%pack_RHS) then
    call fmm_sort_paras_wrt_centre(1_INTK,RHS%paras)
    call fmm_assign_batches(RHS%paras)
  end if

end subroutine sort_centres_and_assign_batches

!-------------------------------------------------------------------------------

subroutine get_LHS_data(scheme,LHS)

  use fmm_qlm_utils, only: fmm_factor_in_dens, fmm_pack_raw_parameters

  implicit none
  type(scheme_paras), intent(in)   :: scheme
  type(raw_mm_data), intent(inout) :: LHS
  integer(INTK) :: i, ndim, qlm_dim, alloc_error

  ! Get LHS parameters from raw data
  if (scheme%pack_LHS) then
    call fmm_pack_raw_parameters(LHS)
  end if

  ! Sync (packed)parameters:(packed)moments mappings
  do i=1,size(LHS%paras)
    LHS%paras(i)%id = i
  end do

  ! LHS preconditioning for T-matrix
  select case(scheme%T_con%LHS_mm_type)
    case(USE_RAW_QLM)
      ! Note that dimensions of paras and qlm,2 may be different
      qlm_dim = size(LHS%qlm,2)
      ndim = size(LHS%qlm,1)
      write(LUPRI,*) 'LHS%qlm_T: Attempting to allocate',max(1,qlm_dim*ndim*8/1000000),'MB of memory...'
      allocate(LHS%qlm_T(ndim,qlm_dim),stat=alloc_error)
      if (alloc_error /= 0) write(LUPRI,*) '... Failed!'
      LHS%qlm_T(:,:) = LHS%qlm(:,:)
    case default
     call fmm_quit('cannot reconcile LHS_mm_type')
  end select

  ! Factorise in density if required
  if (scheme%LHS_dens) then
    call fmm_factor_in_dens(LHS%dens,LHS%qlm_T)
    deallocate(LHS%dens)
    nullify(LHS%dens)
  end if

end subroutine get_LHS_data

!-------------------------------------------------------------------------------

subroutine get_RHS_data(scheme,RHS)

  use fmm_qlm_utils, only: fmm_get_T_sym_qlm, &
                           fmm_pack_raw_moments, &
                           fmm_factor_in_dens

  implicit none
  type(scheme_paras),intent(in)   :: scheme
  type(raw_mm_data),intent(inout) :: RHS
  integer(INTK) :: LMAX, i, qlm_dim, ndim, alloc_error
  logical :: dens
  real(REALK) :: thr

  LMAX = scheme%raw_LMAX

  if (scheme%pack_RHS) then
    ! We now pack all moments belonging to
    ! the same batch (same centre,extent);
    ! Also factor in density and perform
    ! density-based screening if flagged
    dens = scheme%RHS_dens
    thr = scheme%dens_screen_thr
    call fmm_pack_raw_moments(RHS,dens,thr)
  end if

  ndim = (1+LMAX)**2
  qlm_dim = size(RHS%qlm,2)
  write(LUPRI,*) 'RHS%qlm_W: Attempting to allocate',max(1,qlm_dim*ndim*8/1000000),'MB of memory...'
  allocate(RHS%qlm_W(ndim,qlm_dim),staT=alloc_error)
  if (alloc_error /= 0) write(LUPRI,*) '... Failed!'
  RHS%qlm_W(:,:) = RHS%qlm(:,:)

  if (.not.scheme%pack_RHS) then
    if (scheme%RHS_dens) call fmm_factor_in_dens(RHS%dens,RHS%qlm_W)
  end if

  ! Density nolonger required
  if (scheme%RHS_dens) then
    deallocate(RHS%dens)
    nullify(RHS%dens)
  end if

  if (qlm_dim /= size(RHS%paras)) call fmm_quit('error in RHS data')
  ! Resync parameters:moments mappings
  do i=1,qlm_dim
    RHS%paras(i)%id = i
  end do

  ! RHS preconditioning for T-matrix
  select case(scheme%T_con%RHS_mm_type)
    case(USE_RAW_QLM)
      RHS%qlm_T => RHS%qlm_W(:,:)
    case(USE_T_SYM_QLM)
      allocate(RHS%qlm_T(ndim,qlm_dim))
      ! build %qlm_T by rescaling significant %qlm_W
      call fmm_get_T_sym_qlm(LMAX,RHS%qlm_W,RHS%qlm_T)
    case default
      call fmm_quit('cannot reconcile RHS_mm_type')
  end select

end subroutine get_RHS_data

!-------------------------------------------------------------------------------

end module fmm_aux_qlm_builder
