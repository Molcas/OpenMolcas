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

module fmm_T_pair_tests

use fmm_global_paras, only: INTK, REALK, scheme_paras, gen_mm_paras, LHS_RHS_type, T_pair_single, box_mm_paras, DO_FQ, DO_BQ, &
                            DO_NlogN, DO_FMM, ZERO_DIST_TOL, NEAR_FIELD
use fmm_box_utils, only: fmm_NF_boxes, fmm_RFF_boxes
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_init_T_pair_tests, &
          fmm_close_T_pair_tests, &
          fmm_test_and_buffer_T_pair

! Flag to test initialisation
character(len=11), save :: init_tests

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_T_pair_tests(scheme)

  implicit none
  type(scheme_paras), intent(in) :: scheme

  if (scheme%phase == NEAR_FIELD) then
    call fmm_store_test(fmm_test_NF_ext)
  else
    select case (scheme%algorithm)
      case (DO_FQ)
        call fmm_store_test(fmm_test_raw_FF)
      case (DO_BQ)
        call fmm_store_test(fmm_test_FF)
      case (DO_NlogN)
        call fmm_store_test(fmm_test_LFF)
      case (DO_FMM)
        call fmm_store_test(fmm_test_LFF)
      case default
        call fmm_quit('unable to initialise T_pair_tests')
    end select
  end if

  init_tests = 'initialised'

end subroutine fmm_init_T_pair_tests

!-------------------------------------------------------------------------------

subroutine fmm_close_T_pair_tests()

  implicit none
  if (init_tests /= 'initialised') call fmm_quit('must initialise pair_tests!')
  init_tests = ' '

end subroutine fmm_close_T_pair_tests

!-------------------------------------------------------------------------------
! routine to perform T-pair test on two raw or boxed moments, and if they
! pass the pair is sent for evaluation.

subroutine fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)

  use fmm_T_buffer, only: fmm_add_to_T_buffer

  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  integer(INTK), intent(in)      :: weight
  type(T_pair_single) :: T_pair
  logical, external :: fmm_included_pair
  external :: fmm_stored_t_pair_mould

  if (fmm_included_pair(LHS,RHS,id)) then
    ! pour all relevant info into the T-pair mould
    call fmm_stored_t_pair_mould(LHS,RHS,id,weight,T_pair)
    ! pass single T_pair entity to contractors via buffer
    call fmm_add_to_T_buffer(T_pair)
  end if

end subroutine fmm_test_and_buffer_T_pair

!-------------------------------------------------------------------------------
! fmm_test_pass is always TRUE, bypassing normal classical/nonclassical tests

function fmm_test_pass(LHS,RHS,id)

# include "macros.fh"

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_pass

  unused_var(LHS)
  unused_var(RHS)
  unused_var(id)

  fmm_test_pass = .true.

end function fmm_test_pass

!-------------------------------------------------------------------------------
! Logical function to test if two gaussian overlaps are separated by
! more than the sum of their extents.
! fmm_test_ext is TRUE if sufficiently well-separated.

function fmm_test_ext(LHS,RHS,id)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_ext
  real(REALK) :: r_ij(3), ext_ij, r_mod, r_zero

  ext_ij = RHS%raw_paras(id%RHS)%ext+LHS%raw_paras(id%LHS)%ext
  r_ij = RHS%raw_paras(id%RHS)%cntr-LHS%raw_paras(id%LHS)%cntr
  r_mod = r_ij(1)*r_ij(1)+r_ij(2)*r_ij(2)+r_ij(3)*r_ij(3)
  ext_ij = ext_ij*ext_ij
  r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL
  fmm_test_ext = ((r_mod-ext_ij) > r_zero)

end function fmm_test_ext

!-------------------------------------------------------------------------------
! Interaction pair test routine for near-field phase.
! An interaction should be treated by multipoles if in NF boxes and
! sufficiently well-separated.

function fmm_test_NF_ext(LHS,RHS,id)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_NF_ext
  type(box_mm_paras) :: LHStmp, RHStmp

  LHStmp%box = LHS%raw_paras(id%LHS)%box
  RHStmp%box = RHS%raw_paras(id%RHS)%box
  LHStmp%bra = LHS%raw_paras(id%LHS)%bra
  RHStmp%bra = RHS%raw_paras(id%RHS)%bra
  LHStmp%level = 0
  RHStmp%level = 0

  fmm_test_NF_ext = fmm_NF_boxes(LHStmp,RHStmp) .and. fmm_test_ext(LHS,RHS,id)

end function fmm_test_NF_ext

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Far Field"
! based on separation.  Used in simple full quadratic (FQ) algorithm
! when we wish to compare unboxed moments (c.f. fmm_test_FF function)

function fmm_test_raw_FF(LHS,RHS,id)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_raw_FF
  type(box_mm_paras) :: LHStmp, RHStmp

  LHStmp%box = LHS%raw_paras(id%LHS)%box
  RHStmp%box = RHS%raw_paras(id%RHS)%box
  LHStmp%bra = LHS%raw_paras(id%LHS)%bra
  RHStmp%bra = RHS%raw_paras(id%RHS)%bra
  LHStmp%level = 0
  RHStmp%level = 0

  fmm_test_raw_FF = .not.(fmm_NF_boxes(LHStmp,RHStmp))

end function fmm_test_raw_FF

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Far Field"
! based on separation.  Used by Boxed Quadratic (BQ) algorithm.

function fmm_test_FF(LHS,RHS,id)

  implicit none
  type(gen_mm_paras), intent(in) :: LHS, RHS
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_FF

  fmm_test_FF = .not.(fmm_NF_boxes(LHS%box_paras(id%LHS),RHS%box_paras(id%RHS)))

end function fmm_test_FF

!-------------------------------------------------------------------------------
! Logical function to test if two boxes are in each others "Local Far Field"
! based on separation;  used by FMM and NlogN algorithms;
! note that in the NlogN scheme boxes are at different levels, but the
! definition of LFF is such that we can then translate the box at the
! deeper level to the higher level and apply the same criteria as for
! two boxes at the deeper level (assuming branches join at higher levels)
! FIXME: how about branch scheme where branches don't join at higher levels??

function fmm_test_LFF(LHS_in,RHS_in,id)

  use fmm_box_utils, only: fmm_translate_to_common_grid

  implicit none
  type(gen_mm_paras), intent(in) :: LHS_in, RHS_in
  type(LHS_RHS_type), intent(in) :: id
  logical :: fmm_test_LFF

  ! introduce these so we can translate to common grid if needed
  type(box_mm_paras) :: LHS, RHS

  LHS = LHS_in%box_paras(id%LHS)
  RHS = RHS_in%box_paras(id%RHS)
  if (LHS%level /= RHS%level) then
    ! All NF, LFF and RFF definitions are based on boxes in ONE grid.
    ! We must therefore translate the boxes into a common grid first.
    ! Of course, this can be avoided by passing the transformed paras.
    call fmm_translate_to_common_grid(LHS,RHS)
  end if

  fmm_test_LFF = .false.
  if (fmm_NF_boxes(LHS,RHS)) return
  if (fmm_RFF_boxes(LHS,RHS)) return
  fmm_test_LFF = .true.

end function fmm_test_LFF

!-------------------------------------------------------------------------------

end module fmm_T_pair_tests
