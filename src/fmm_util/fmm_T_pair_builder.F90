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

module fmm_T_pair_builder

use fmm_global_paras, only: INTK, scheme_paras, gen_mm_paras, LHS_RHS_type, LHS_raw_RHS_raw, LHS_box_RHS_box
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_T_pair_builder, &
          fmm_close_T_pair_builder, &
          fmm_gen_local_T_pairs, &
          fmm_gen_nonlocal_T_pairs

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_T_pair_builder(scheme,pair_type)

  use fmm_T_pair_mould, only: fmm_init_T_pair_mould
  use fmm_T_pair_tests, only: fmm_init_T_pair_tests
  use fmm_T_buffer, only: fmm_open_T_buffer

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK), intent(in)      :: pair_type

  call fmm_init_T_pair_mould(scheme,pair_type)
  call fmm_init_T_pair_tests(scheme)
  call fmm_open_T_buffer(scheme)

end subroutine fmm_init_T_pair_builder

!-------------------------------------------------------------------------------

subroutine fmm_close_T_pair_builder()

  use fmm_T_pair_mould, only: fmm_close_T_pair_mould
  use fmm_T_pair_tests, only: fmm_close_T_pair_tests
  use fmm_T_buffer, only: fmm_close_T_buffer

  implicit none

  call fmm_close_T_pair_mould()
  call fmm_close_T_pair_tests()
  call fmm_close_T_buffer()

end subroutine fmm_close_T_pair_builder

!-------------------------------------------------------------------------------

subroutine fmm_gen_local_T_pairs(LHS,RHS,pair_type)

  use fmm_T_pair_tests, only: fmm_test_and_buffer_T_pair
  use fmm_local_search, only: fmm_get_local_paras

  implicit none
  type(gen_mm_paras), intent(inout) :: LHS, RHS
  integer(INTK), intent(in)         :: pair_type
  type(gen_mm_paras) :: RHS_local
  type(LHS_RHS_type) :: id
  integer(INTK) :: weight, i, j, ndim

  weight = 1
  nullify(RHS_local%box_paras,RHS_local%raw_paras)

  select case (pair_type)
    case (LHS_raw_RHS_raw)
      do j=1,size(LHS%raw_paras)
        id%LHS = j
        do i=1,size(RHS%raw_paras)
          id%RHS = i
          call fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
        end do
      end do
    case (LHS_box_RHS_box)
      do j=1,size(LHS%box_paras)
        call fmm_get_local_paras(j,RHS,pair_type,RHS_local,ndim)
        if (ndim == 0) cycle
        id%LHS = j
        do i=1,size(RHS_local%box_paras)
          id%RHS = i
          call fmm_test_and_buffer_T_pair(LHS,RHS_local,id,weight)
        end do
        if (associated(RHS_local%box_paras)) then
          deallocate(RHS_local%box_paras)
          nullify(RHS_local%box_paras)
        end if
      end do
    case default
      call fmm_quit('cannot reconcile requested T_pair type!')
  end select

end subroutine fmm_gen_local_T_pairs

!-------------------------------------------------------------------------------

subroutine fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)

  use fmm_T_pair_tests, only: fmm_test_and_buffer_T_pair

  implicit none
  type(gen_mm_paras), intent(inout) :: LHS, RHS
  integer(INTK), intent(in)         :: pair_type
  type(LHS_RHS_type) :: id
  integer(INTK) :: weight, i, j

  weight = 1

  select case (pair_type)
    case (LHS_raw_RHS_raw)
      do j=1,size(LHS%raw_paras)
        id%LHS = j
        do i=1,size(RHS%raw_paras)
          id%RHS = i
          call fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
        end do
      end do
    case (LHS_box_RHS_box)
      do j=1,size(LHS%box_paras)
        id%LHS = j
        do i=1,size(RHS%box_paras)
          id%RHS = i
          call fmm_test_and_buffer_T_pair(LHS,RHS,id,weight)
        end do
      end do
    case default
      call fmm_quit('cannot reconcile requested T_pair type!')
  end select

end subroutine fmm_gen_nonlocal_T_pairs

!-------------------------------------------------------------------------------

end module fmm_T_pair_builder
