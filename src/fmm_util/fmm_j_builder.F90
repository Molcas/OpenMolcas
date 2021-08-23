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

module fmm_J_builder

use fmm_global_paras, only: INTK, REALK, scheme_paras, raw_mm_data, id_node, ELECTRONIC_ONLY, NUCLEAR_ONLY, ALL_MOMENTS, Half
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_get_E_from_Vff, &
          fmm_get_E_from_pkd_Vff, &
          fmm_get_J_from_Vff, &
          fmm_get_J_from_pkd_Vff

contains

!-------------------------------------------------------------------------------
! Trivial check that the number of moments matches the number of potentials

subroutine fmm_verify_Vff_input(scheme,LHS_mms,Vff,J_or_E)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  type(raw_mm_data), intent(in)  :: LHS_mms
  real(REALK), intent(in)        :: Vff(:,:)
  character, intent(in)          :: J_or_E
  logical :: A, B, C

  A = (size(LHS_mms%paras) /= size(Vff,2))
  if (A) call fmm_quit('incompatible SIZE of Vff and LHS moments!')

  if (J_or_E == 'J') then
    A = (scheme%LHS_mm_range == NUCLEAR_ONLY)
    B = (scheme%RHS_mm_range == NUCLEAR_ONLY)
    C = (scheme%LHS_mm_range == ALL_MOMENTS)
    if ((A .and. B) .or. C) call fmm_quit('mm_ranges invalid')
  end if

end subroutine fmm_verify_Vff_input

!-------------------------------------------------------------------------------

subroutine fmm_get_E_with_text(scheme,energy,text)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  real(REALK), intent(inout)     :: energy
  character(len=*), intent(out)  :: text
  logical :: A, B, C, D

  A = (scheme%LHS_mm_range == ELECTRONIC_ONLY)
  B = (scheme%RHS_mm_range == ELECTRONIC_ONLY)
  C = (scheme%LHS_mm_range == NUCLEAR_ONLY)
  D = (scheme%RHS_mm_range == NUCLEAR_ONLY)

  if (scheme%LHS_mm_range == scheme%RHS_mm_range) then
    energy = half*energy
    text = 'total classical Coulomb energy'
    if (A) text = 'classical Coulomb electronic energy'
    if (C) text = 'classical Coulomb nuclear repulsion'
  else if (A .or. B) then
    if (C .or. D) then
      text = 'classical Coulomb nuclear attraction'
    else
      text = 'e-n + 2*(e-e) energy'
    end if
  else
    ! range is different and neither is ELECTRONIC_ONLY
    text = 'e-n + 2*(n-n) energy'
  end if

end subroutine fmm_get_E_with_text

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

subroutine fmm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,text)

  implicit none
  type(scheme_paras), intent(IN) :: scheme
  type(raw_mm_data), intent(IN)  :: LHS_mms
  real(REALK), intent(IN)        :: Vff(:,:)
  real(REALK), intent(OUT)       :: energy
  character(len=*), intent(OUT)  :: text

  real(REALK) :: g
  integer(INTK) :: u, v, lm_max

  call fmm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

  ! although Vff should be the same size, we test for generality
  if (size(LHS_mms%qlm_T,1) /= size(Vff,1)) call fmm_quit('mm_get_E_from_Vff:2')
  lm_max = min(size(LHS_mms%qlm_T,1),size(Vff,1))
  do u=1,size(LHS_mms%paras)
    v = LHS_mms%paras(u)%id
    g = dot_product(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
    energy = energy+g
  end do

  call fmm_get_E_with_text(scheme,energy,text)

end subroutine fmm_get_E_from_Vff

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
! on the basis that the LHS parameters and Vff may have been packed
! into batches (and must be "expanded").
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

subroutine fmm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,text)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  type(raw_mm_data), intent(in)  :: LHS_mms
  real(REALK), intent(in)        :: Vff(:,:)
  real(REALK), intent(out)       :: energy
  character(len=*), intent(out)  :: text

  real(REALK) :: g
  integer(INTK) :: u, v, w, lm_max
  type(id_node), pointer :: batch_map

  call fmm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

  lm_max = min(size(LHS_mms%qlm_T,1),size(Vff,1))
  packed_loop: do u=1,size(LHS_mms%paras)

    v = LHS_mms%paras(u)%id ! LHS packed (batch) parameters
    batch_map => LHS_mms%batch_map(v)%head

    batch_members2: do ! over batch list until pointer disassociated
      w = batch_map%id ! raw LHS moment ID
      g = dot_product(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
      energy = energy+g
      ! only do next raw item in batch list if it exists
      if (.not. associated(batch_map%next)) exit batch_members2
      batch_map => batch_map%next
    end do batch_members2

  end do packed_loop

  call fmm_get_E_with_text(scheme,energy,text)

end subroutine fmm_get_E_from_pkd_Vff

!-------------------------------------------------------------------------------
! Build J_matrix components from contracion of the LHS moments and potentials
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! This choice of Vff is fine here since the e-e interactions require a
! factor of 2 in the J-matrix, but NOT the e-n interactions.

subroutine fmm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  type(raw_mm_data), intent(in)  :: LHS_mms
  real(REALK), intent(in)        :: Vff(:,:)
  real(REALK), intent(out)       :: J_matrix(:,:)

  real(REALK) :: g
  integer(INTK) :: u, v, i, j, lm_max

  call fmm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

  ! although Vff should be the same size, we test for generality
  if (size(LHS_mms%qlm_T,1) /= size(Vff,1)) call fmm_quit('mm_get_J_from_Vff:2')
  lm_max = min(size(LHS_mms%qlm_T,1),size(Vff,1))
  do u=1,size(LHS_mms%paras)
    v = LHS_mms%paras(u)%id
    g = dot_product(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
    i = LHS_mms%J_indices(v)%i_indx
    j = LHS_mms%J_indices(v)%j_indx
    J_matrix(i,j) = J_matrix(i,j)+g
    if (i /= j) J_matrix(j,i) = J_matrix(j,i)+g
  end do

end subroutine fmm_get_J_from_Vff

!-------------------------------------------------------------------------------
! Build J_matrix components from contracion of the LHS moments and potentials
! This routine recognises that the LHS parameters may have been packed
! and thus expands the potential over LHS batches.
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments

subroutine fmm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  type(raw_mm_data), intent(in)  :: LHS_mms
  real(REALK), intent(in)        :: Vff(:,:)
  real(REALK), intent(out)       :: J_matrix(:,:)

  real(REALK) :: g
  integer(INTK) :: u, v, w, i, j, lm_max
  type(id_node), pointer :: batch_map

  call fmm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

  ! Vff should now be the same size as the raw LHS parameters;
  ! However, the size of LHS %qlm_T will be larger;
  ! We need to "expand" Vff using the batch mapping.

  lm_max = min(size(LHS_mms%qlm_T,1),size(Vff,1))
  packed_loop: do u=1,size(LHS_mms%paras)

    v = LHS_mms%paras(u)%id ! LHS packed (batch) moments
    batch_map => LHS_mms%batch_map(v)%head

    batch_members: do ! over batch list until pointer disassociated
      w = batch_map%id ! raw LHS moment ID
      g = dot_product(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
      i = LHS_mms%J_indices(w)%i_indx
      j = LHS_mms%J_indices(w)%j_indx
      J_matrix(i,j) = J_matrix(i,j)+g
      !if (i /= j) J_matrix(j,i) = J_matrix(j,i)+g
      ! only do next raw item in batch list if it exists
      if (.not. associated(batch_map%next)) exit batch_members
      batch_map => batch_map%next
    end do batch_members
  end do packed_loop

end subroutine fmm_get_J_from_pkd_Vff

!-------------------------------------------------------------------------------

end module fmm_J_builder
