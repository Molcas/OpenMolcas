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

module fmm_T_pair_mould

use fmm_global_paras, only: INTK, gen_mm_paras, LHS_RHS_type, T_pair_single, scheme_paras, LHS_raw_RHS_raw, LHS_box_RHS_box, One
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_init_T_pair_mould, &
          fmm_close_T_pair_mould

integer(INTK), save :: LHS_LMAX, RHS_LMAX
! flag to test initialisation
character(len=11), save :: fmm_init_mould

contains

!-------------------------------------------------------------------------------

subroutine fmm_close_T_pair_mould()

  implicit none
  if (fmm_init_mould /= 'initialised') call fmm_quit('mm_T_pair_mould init')
  fmm_init_mould = ' '
  LHS_LMAX = 0
  RHS_LMAX = 0

end subroutine fmm_close_T_pair_mould

!-------------------------------------------------------------------------------

subroutine fmm_set_T_pair_basics(LHS,RHS,id,weight,T_pair)

# include "macros.fh"

  implicit none
  type(gen_mm_paras), intent(in)   :: LHS, RHS
  type(LHS_RHS_type), intent(in)   :: id
  integer(INTK), intent(in)        :: weight
  type(T_pair_single), intent(out) :: T_pair

  unused_var(LHS)
  unused_var(RHS)
  unused_var(id)

  T_pair%N_or_T = 'N' ! not used for T contractions, so null
  ! init. normalisation scalar for other options
  T_pair%paras%ratio = one
  ! here, T_pair%r_ab is the unnormalised vector

  ! weighting to account for double-counting in T-pair generation
  T_pair%paras%weight = weight

  T_pair%LMAX = max(T_pair%paras%LHS_LMAX,T_pair%paras%RHS_LMAX)
  T_pair%lm_max = (1+T_pair%LMAX)**2

end subroutine fmm_set_T_pair_basics

!-------------------------------------------------------------------------------

subroutine fmm_set_RR_paras(LHS,RHS,id,T_pair)

  implicit none
  type(gen_mm_paras), intent(in)   :: LHS, RHS
  type(LHS_RHS_type), intent(in)   :: id
  type(T_pair_single), intent(out) :: T_pair

  ! interaction vector for building T-matrix
  T_pair%r_ab = RHS%raw_paras(id%RHS)%cntr-LHS%raw_paras(id%LHS)%cntr
  ! indices to map back to actual moments in separate array
  T_pair%paras%LHS_id = LHS%raw_paras(id%LHS)%id
  T_pair%paras%RHS_id = RHS%raw_paras(id%RHS)%id
  ! check that paras:moments mapping was built
  if (T_pair%paras%LHS_id == 0) call fmm_quit('LHS paras:moments mapping')
  if (T_pair%paras%RHS_id == 0) call fmm_quit('RHS paras:moments mapping')

end subroutine fmm_set_RR_paras

!-------------------------------------------------------------------------------

subroutine fmm_set_BB_paras(LHS,RHS,id,T_pair)

  implicit none
  type(gen_mm_paras), intent(in)   :: LHS, RHS
  type(LHS_RHS_type), intent(in)   :: id
  type(T_pair_single), intent(out) :: T_pair

  ! interaction vector for building T-matrix
  T_pair%r_ab = RHS%box_paras(id%RHS)%cntr-LHS%box_paras(id%LHS)%cntr
  ! indices to map back to actual moments in separate array
  T_pair%paras%LHS_id = LHS%box_paras(id%LHS)%id
  T_pair%paras%RHS_id = RHS%box_paras(id%RHS)%id
  ! check that paras:moments mapping was built
  if (T_pair%paras%LHS_id == 0) call fmm_quit('LHS paras:moments mapping')
  if (T_pair%paras%RHS_id == 0) call fmm_quit('RHS paras:moments mapping')

end subroutine fmm_set_BB_paras

!-------------------------------------------------------------------------------

subroutine fmm_set_LHS_LMAX(X,Y,T_pair)

# include "macros.fh"

  implicit none
  type(gen_mm_paras), intent(in)   :: X        ! dummy variables
  type(LHS_RHS_type), intent(in)   :: Y        ! dummy variables
  type(T_pair_single), intent(out) :: T_pair

  unused_var(X)
  unused_var(Y)

  T_pair%paras%LHS_LMAX = LHS_LMAX

end subroutine fmm_set_LHS_LMAX

!-------------------------------------------------------------------------------

subroutine fmm_set_RHS_LMAX(X,Y,T_pair)

# include "macros.fh"

  implicit none
  type(gen_mm_paras), intent(in)   :: X        ! dummy variables
  type(LHS_RHS_type), intent(in)   :: Y        ! dummy variables
  type(T_pair_single), intent(out) :: T_pair

  unused_var(X)
  unused_var(Y)

  T_pair%paras%RHS_LMAX = RHS_LMAX

end subroutine fmm_set_RHS_LMAX

!-------------------------------------------------------------------------------
! This routine directs the saving of functions in fmm_proc_selector.c
! relevant to the formation of a T-pair entity.
! There are 4 parts to the making of a T-pair once raw data has been
! suppplied from the classical interaction search algorithm.
! These are all called consecutively from the C-code via
! fmm_stored_t_pair_mould.
! The 4 functions saved are selected here and depend on whether
! pure boxed moments, or unboxed moments, are being contracted.

subroutine fmm_init_T_pair_mould(scheme,pair_type)

  implicit none
  type(scheme_paras), intent(in) :: scheme
  integer(INTK), intent(in)      :: pair_type
  external fmm_store_t_pair_mould1   ! raw/box dependent variables
  external fmm_store_t_pair_mould2   ! set LMAX (LHS)
  external fmm_store_t_pair_mould3   ! set LMAX (RHS)
  external fmm_store_t_pair_mould4   ! common to all T-pair builds

  call fmm_store_t_pair_mould2(fmm_set_LHS_LMAX)
  call fmm_store_t_pair_mould3(fmm_set_RHS_LMAX)
  call fmm_store_t_pair_mould4(fmm_set_T_pair_basics)
  select case(pair_type)
    case(LHS_raw_RHS_raw)
      LHS_LMAX = scheme%raw_LMAX
      RHS_LMAX = scheme%raw_LMAX
      call fmm_store_t_pair_mould1(fmm_set_RR_paras)
    case(LHS_box_RHS_box)
      LHS_LMAX = scheme%trans_LMAX
      RHS_LMAX = scheme%trans_LMAX
      call fmm_store_t_pair_mould1(fmm_set_BB_paras)
    case default
      call fmm_quit('cannot recognise T_pair type!')
  end select

  fmm_init_mould = 'initialised'

end subroutine fmm_init_T_pair_mould

!-------------------------------------------------------------------------------

end module fmm_T_pair_mould
