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
! Copyright (C) 2002, Mark A. Watson                                   *
!***********************************************************************
! fmm_driver:
! General module for driving multipole method calculations.

module fmm_driver

use fmm_global_paras, only: INTK, REALK, LUPRI, raw_mm_data, scheme_paras, Zero, ELECTRONIC_ONLY, NUCLEAR_ONLY, ALL_MOMENTS, GFC_FMM
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_build_J_matrix, &
          fmm_get_multipole_potential

! The multipole potential
real(REALK), pointer, save :: Vff(:,:)
! The (modified) multipole expansions
type(raw_mm_data), save :: LHS_mms, RHS_mms

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_driver(scheme,dens)

  use fmm_qlm_builder, only: fmm_get_raw_qlm
  use fmm_aux_qlm_builder, only: fmm_get_aux_qlm

  implicit none
  type(scheme_paras), intent(inout) :: scheme
  real(REALK), intent(in)           :: dens(:,:)

  nullify(Vff)

  ! Get basic multipole data
  call fmm_get_raw_qlm(scheme,dens,LHS_mms,RHS_mms)
  ! Now ensure appropriate prefactors, normalisation and packing
  call fmm_get_aux_qlm(scheme,LHS_mms,RHS_mms)
  ! Allocate the far field potential
  call fmm_allocate_Vff(scheme)

end subroutine fmm_init_driver

!-------------------------------------------------------------------------------

subroutine fmm_free_driver()

  use fmm_qlm_builder, only: fmm_deallocate_qlm

  implicit none
  deallocate(Vff)
  nullify(Vff)
  call fmm_deallocate_qlm(LHS_mms,RHS_mms)

end subroutine fmm_free_driver

!-------------------------------------------------------------------------------

subroutine fmm_allocate_Vff(scheme)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK) :: lm_dim, mms_dim, alloc_error

  if (.not. associated(LHS_mms%paras)) call fmm_quit('mms ptrs not set in fmm_driver!')
  if (associated(Vff)) call fmm_quit('Vff should NOT be allocated already!')

  mms_dim = size(LHS_mms%paras)
  ! Note we shouldn't use this array for translated, BOXED potentials
  lm_dim = (1+scheme%raw_LMAX)**2
  if (scheme%job_type == GFC_FMM) lm_dim = 1

  write(LUPRI,*) 'Vff: Attempting to allocate',max(1,lm_dim*mms_dim*8/1000000),'MB of memory...'
  allocate(Vff(lm_dim,mms_dim),stat=alloc_error)
  if (alloc_error /= 0) write(LUPRI,*) '... Failed!'

  ! Must zero out since Vff is built additively in general
  Vff(:,:) = zero

end subroutine fmm_allocate_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_J_via_raw_potentials(scheme,dens,J_matrix,energy,txt)

  use fmm_Vff_driver, only: fmm_get_Vff
  use fmm_J_builder, only: fmm_get_J_from_Vff, &
                           fmm_get_J_from_pkd_Vff, &
                           fmm_get_E_from_pkd_Vff, &
                           fmm_get_E_from_Vff
  use fmm_qlm_utils, only: fmm_factor_in_dens

  implicit none
  type(scheme_paras), intent(inout) :: scheme
  real(REALK), intent(in)           :: dens(:,:)
  real(REALK), intent(out)          :: J_matrix(:,:)
  real(REALK), intent(out)          :: energy
  character(len=*), intent(out)     :: txt

  ! We only have the density on RHS when getting J-matrix
  ! via far-field potential
  scheme%LHS_dens = .false.
  scheme%RHS_dens = .true.

  ! Prepare moments and allocate potential
  call fmm_init_driver(scheme,dens)
  ! Get potential
  call fmm_get_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

  ! Get J-matrix
  J_matrix = zero
  energy = zero
  if (scheme%pack_LHS) then
    call fmm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)
    ! Get energy after factoring in density to LHS
    call fmm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
    call fmm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,txt)
  else
    call fmm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
    ! Get energy after factoring in density to LHS
    call fmm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
    call fmm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,txt)
  end if

  call fmm_free_driver()

end subroutine fmm_get_J_via_raw_potentials

!-------------------------------------------------------------------------------

subroutine fmm_build_J_matrix(n_el,dens,J_matrix)

   use fmm_stats, only: fmm_print_stats
   use fmm_scheme_builder, only: fmm_get_scheme
   use fmm_utils, only: fmm_second, TIMTXT

   implicit none
   character(len=6), intent(in)  :: n_el
   real(REALK), intent(in)       :: dens(:,:)
   real(REALK), intent(out)      :: J_matrix(:,:)

   type(scheme_paras), pointer :: scheme
   character(len=36) :: E_text
   real(REALK) :: energy, T0, TTOT

   T0 = fmm_second()

   call fmm_get_scheme(scheme)

  select case(n_el)
    case('ONE_EL')
      call fmm_quit('nuclear moments not available!')
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY
    case('TWO_EL')
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
    case('FULL_J')
      call fmm_quit('nuclear moments not available!')
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ALL_MOMENTS
    case default
      call fmm_quit('require 1, 2, or full J_matrix build!')
  end select

  call fmm_get_J_via_raw_potentials(scheme,dens,J_matrix,energy,E_text)
  write(LUPRI,'(1X,A," = ",E20.12)') E_text,energy

  TTOT = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_J_matrix',TTOT,LUPRI)
  call fmm_print_stats()

end subroutine fmm_build_J_matrix

!-------------------------------------------------------------------------------
! Note that the potential returned by this routine can only be
! contracted with appropriately scaled multipole moments

subroutine fmm_get_multipole_potential(mode,dens,potential)

  use fmm_stats, only: fmm_print_stats
  use fmm_scheme_builder, only: fmm_get_scheme
  use fmm_boundary, only: fmm_opt_near_field
  use fmm_Vff_driver, only: fmm_get_Vff
  use fmm_utils, only: fmm_second, TIMTXT

  implicit none
  integer(INTK), intent(in) :: mode
  real(REALK), intent(in)   :: dens(:,:)
  real(REALK), intent(out)  :: potential(:,:)

  type(scheme_paras), pointer :: scheme
  real(REALK) :: T0, TTOT
  integer(INTK) :: lmdim

  T0 = fmm_second()

  call fmm_get_scheme(scheme)

  scheme%LHS_mm_range = NUCLEAR_ONLY
  scheme%RHS_mm_range = ELECTRONIC_ONLY
  scheme%LHS_dens = .false.
  scheme%RHS_dens = .true.
  scheme%pack_LHS = .false.

  ! Prepare moments and allocate potential
  call fmm_init_driver(scheme,dens)

  ! Test if we can skip the near-field interactions
  if (mode == GFC_FMM) then
    call fmm_opt_near_field(scheme,LHS_mms%paras,RHS_mms%paras)
  end if

  ! Get potential
  call fmm_get_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

  ! Note we assume here that the LHS hasn't been rearranged
  lmdim = size(potential,1)
  if (size(potential,2) /= size(Vff,2)) call fmm_quit('bounds: potential')
  potential(:,:) = Vff(1:lmdim,:)

  call fmm_free_driver()

  TTOT = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_multipole_potential',TTOT,LUPRI)
  call fmm_print_stats()

end subroutine fmm_get_multipole_potential

!-------------------------------------------------------------------------------

end module fmm_driver
