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

module fmm_Vff_driver

use fmm_global_paras, only: INTK, REALK, LUPRI, raw_mm_data, raw_mm_paras, scheme_paras, gen_mm_paras, box_mm_paras, &
                            LHS_raw_RHS_raw, LHS_box_RHS_box, DO_NULL, DO_FQ, DO_BQ, DO_NlogN, DO_FMM, TOP_LEVEL, FAR_FIELD, &
                            NEAR_FIELD, Zero
use fmm_stats, only: stat_NF_not_FF
use fmm_T_contractors, only: fmm_set_T_con_ptrs, &
                             fmm_select_T_con, &
                             fmm_init_T_contractors, &
                             fmm_free_T_contractors
use fmm_W_contractors, only: fmm_init_W_contractors, &
                             fmm_free_W_contractors
use fmm_T_pair_builder, only: fmm_init_T_pair_builder, &
                              fmm_close_T_pair_builder, &
                              fmm_gen_local_T_pairs, &
                              fmm_gen_nonlocal_T_pairs

use fmm_box_utils, only: fmm_deepest_level
use fmm_box_builder, only: fmm_get_box_qlm_at_level, &
                           fmm_get_box_paras_at_level
use fmm_W_pair_builder, only: fmm_get_raw_Vff_from_boxed_Vff
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_get_Vff

contains

!-------------------------------------------------------------------------------

subroutine fmm_verify_data_in(RHS_in)

  implicit none
  type(raw_mm_data), intent(in) :: RHS_in
  logical :: fail

  fail = .false.
  if (.not. associated(RHS_in%paras)) fail = .true.
  if (.not. associated(RHS_in%qlm_T)) fail = .true.
  if (.not. associated(RHS_in%qlm_W)) fail = .true.
  if (size(RHS_in%paras) /= size(RHS_in%qlm_T,2)) fail = .true.
  if (size(RHS_in%paras) /= size(RHS_in%qlm_W,2)) fail = .true.
  if (fail) call fmm_quit('mms pointers sent incorrectly to fmm_Vff_driver')

end subroutine fmm_verify_data_in

!-------------------------------------------------------------------------------

subroutine fmm_get_Vff(scheme,LHS_paras,RHS_mms,Vff)

  use fmm_box_builder, only: fmm_init_box_builder, fmm_free_box_builder
  use fmm_local_search, only: fmm_init_local_search
  use fmm_local_search, only: fmm_free_local_search

  implicit none
  type(raw_mm_paras), intent(inout) :: LHS_paras(:)
  type(raw_mm_data), intent(inout)  :: RHS_mms
  type(scheme_paras), intent(inout) :: scheme
  real(REALK), target, intent(out)  :: Vff(:,:)

  write(LUPRI,*) 'Computing multipole potential...'
  call fmm_verify_data_in(RHS_mms)

  call fmm_init_box_builder(LHS_paras,RHS_mms,scheme)
  if (scheme%algorithm == DO_FMM) call fmm_init_local_search(scheme)

  scheme%phase = FAR_FIELD
  stat_NF_not_FF = .false.
  select case(scheme%algorithm)
    case(DO_NULL)
      continue
    case(DO_FQ)
      call fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)
    case(DO_BQ)
      call fmm_get_BQ_Vff(scheme,LHS_paras,Vff)
    case(DO_NlogN)
      call fmm_get_NlogN_Vff(scheme,LHS_paras,Vff)
    case(DO_FMM)
      call fmm_get_FMM_Vff(scheme,LHS_paras,Vff)
    case default
      call fmm_quit('invalid algorithm requested!')
  end select

  if (scheme%include_near_field) then
    ! Also compute multipole interactions within NF boxes
    scheme%phase = NEAR_FIELD
    stat_NF_not_FF = .true.
    call fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)
  end if

  call fmm_free_box_builder()
  call fmm_free_local_search()

end subroutine fmm_get_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)

  use fmm_utils, only: fmm_second, TIMTXT

  implicit none
  type(raw_mm_paras), target, intent(in) :: LHS_paras(:)
  type(raw_mm_data), intent(in)          :: RHS_mms
  type(scheme_paras), intent(in)         :: scheme
  real(REALK), target, intent(inout)     :: Vff(:,:)

  integer(INTK), parameter :: pair_type = LHS_raw_RHS_raw
  type(gen_mm_paras) :: LHS, RHS
  real(REALK) :: T0, TTOTFQ

  T0 = fmm_second()
  nullify(LHS%raw_paras,LHS%box_paras)
  nullify(RHS%raw_paras,RHS%box_paras)

  LHS%raw_paras => LHS_paras(:)
  RHS%raw_paras => RHS_mms%paras(:)

  ! select the T-contractor to be stored/called via C code;
  call fmm_select_T_con(scheme)
  ! set T_contractor pointers
  call fmm_set_T_con_ptrs(Vff,RHS_mms%qlm_T)
  call fmm_init_T_contractors(scheme)

  ! generate full multipole potential
  call fmm_init_T_pair_builder(scheme,pair_type)
  call fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
  call fmm_close_T_pair_builder()
  call fmm_free_T_contractors()

  TTOTFQ = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_FQ_Vff',TTOTFQ,LUPRI)

end subroutine fmm_get_FQ_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_BQ_Vff(scheme,LHS_raw_paras,Vff)

  use fmm_utils, only: fmm_second, TIMTXT

  implicit none
  type(raw_mm_paras), target, intent(in) :: LHS_raw_paras(:)
  type(scheme_paras), intent(in)         :: scheme
  real(REALK), target, intent(inout)     :: Vff(:,:)

  type(gen_mm_paras) :: LHS, RHS
  real(REALK), pointer :: qlm_T(:,:)
  real(REALK), pointer :: boxed_Vff(:,:)
  integer(INTK) :: level, LMAX, lm_dim, foo
  integer(INTK), parameter :: pair_type = LHS_box_RHS_box
  real(REALK) :: T0, TTOT

  T0 = fmm_second()
  nullify(qlm_T,boxed_Vff)
  nullify(LHS%raw_paras,LHS%box_paras)
  nullify(RHS%raw_paras,RHS%box_paras)

  call fmm_init_W_contractors(scheme%trans_LMAX)

  level = fmm_deepest_level(scheme)
  if (level < TOP_LEVEL) return
  LMAX = scheme%trans_LMAX
  lm_dim = (1+LMAX)**2

  ! set up LHS and RHS boxed parameters
  call fmm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
  call fmm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
  call fmm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','free')

  ! allocate temporary boxed potential
  foo = size(LHS%box_paras)
  allocate(boxed_Vff(lm_dim,foo))
  boxed_Vff = zero

  ! select the T-contractor to be stored/called via C code
  call fmm_select_T_con(scheme)
  ! set T_contractor pointers
  call fmm_set_T_con_ptrs(boxed_Vff,qlm_T)
  call fmm_init_T_contractors(scheme)

  ! generate FF Vff
  call fmm_init_T_pair_builder(scheme,pair_type)
  call fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
  call fmm_close_T_pair_builder()
  call fmm_free_T_contractors()

  ! translate boxed potential to raw LHS centres
  call fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
  deallocate(boxed_Vff)

  TTOT = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_BQ_Vff',TTOT,LUPRI)

  call fmm_free_W_contractors()

end subroutine fmm_get_BQ_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_NlogN_Vff(scheme,LHS_raw_paras,Vff)

  use fmm_utils, only: fmm_second, TIMTXT

  implicit none
  type(raw_mm_paras), target, intent(in) :: LHS_raw_paras(:)
  type(scheme_paras), intent(in)         :: scheme
  real(REALK), target, intent(inout)     :: Vff(:,:)

  type(gen_mm_paras) :: LHS, RHS
  real(REALK), pointer :: qlm_T(:,:)
  real(REALK), pointer :: boxed_Vff(:,:)
  integer(INTK) :: LHS_level, RHS_level, deepest_level
  integer(INTK) :: LMAX, lm_dim, foo
  integer(INTK), parameter :: pair_type = LHS_box_RHS_box
  real(REALK) :: T0, TTOT

  T0 = fmm_second()
  nullify(qlm_T,boxed_Vff)
  nullify(LHS%raw_paras,LHS%box_paras)
  nullify(RHS%raw_paras,RHS%box_paras)

  call fmm_init_W_contractors(scheme%trans_LMAX)

  deepest_level = fmm_deepest_level(scheme)
  if (deepest_level < TOP_LEVEL) return
  LMAX = scheme%trans_LMAX
  lm_dim = (1+LMAX)**2

  ! select the T-contractor to be stored/called via C code
  call fmm_select_T_con(scheme)
  call fmm_init_T_contractors(scheme)
  ! get LHS boxed MM parameters
  LHS_level = deepest_level
  call fmm_get_box_paras_at_level(LHS_level,scheme,LHS%box_paras,'LHS')
  ! allocate boxed potential based on number of LHS boxes
  foo = size(LHS%box_paras)
  allocate(boxed_Vff(lm_dim,foo))
  boxed_Vff = zero

  ! generate far-field potential at deepest box centres
  do RHS_level=deepest_level,TOP_LEVEL,-1
    ! get packed RHS MM paras and translated moments
    call fmm_get_box_paras_at_level(RHS_level,scheme,RHS%box_paras,'RHS')
    call fmm_get_box_qlm_at_level(RHS_level,scheme,qlm_T,'RHS','free')
    ! set T_contractor pointers
    call fmm_set_T_con_ptrs(boxed_Vff,qlm_T)
    ! generate LFF potential at prescribed centres
    call fmm_init_T_pair_builder(scheme,pair_type)
    !FIXME: would be nice to add NlogN to local search algorithm
    call fmm_gen_nonlocal_T_pairs(LHS,RHS,pair_type)
    call fmm_close_T_pair_builder()
  end do
  call fmm_free_T_contractors()

  ! translate boxed potential to raw LHS centres
  call fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
  deallocate(boxed_Vff)

  TTOT = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_NlogN_Vff',TTOT,LUPRI)

  call fmm_free_W_contractors()

end subroutine fmm_get_NlogN_Vff

!-------------------------------------------------------------------------------

subroutine fmm_get_FMM_Vff(scheme,LHS_raw_paras,Vff)

  use fmm_W_pair_builder, only: fmm_translate_parents_Vff
  use fmm_utils, only: fmm_second, TIMTXT

  implicit none
  type(raw_mm_paras), target, intent(in) :: LHS_raw_paras(:)
  type(scheme_paras), intent(in)         :: scheme
  real(REALK), intent(inout)             :: Vff(:,:)

  type(gen_mm_paras) :: LHS, RHS
  type(box_mm_paras), pointer :: p_box_paras(:)   ! parents data
  real(REALK), pointer :: qlm_T(:,:)
  real(REALK), pointer :: Vff_tmp1(:,:)
  real(REALK), pointer :: Vff_tmp2(:,:)
  real(REALK), pointer :: Vff_ptr(:,:)
  integer(INTK) :: level, l_up, deepest_level
  integer(INTK) :: LMAX, lm_dim, foo
  integer(INTK) :: use_Vff_tmp
  integer(INTK), parameter :: pair_type = LHS_box_RHS_box
  real(REALK) :: T0, TTOT

  T0 = fmm_second()
  nullify(qlm_T,Vff_tmp1,Vff_tmp2,Vff_ptr)
  nullify(LHS%raw_paras,LHS%box_paras)
  nullify(RHS%raw_paras,RHS%box_paras)

  call fmm_init_W_contractors(scheme%trans_LMAX)

  use_Vff_tmp = 0
  deepest_level = fmm_deepest_level(scheme)
  if (deepest_level < TOP_LEVEL) return
  LMAX = scheme%trans_LMAX
  lm_dim = (LMAX+1)**2

  ! select the T-contractor to be stored/called via C code
  call fmm_select_T_con(scheme)
  call fmm_init_T_contractors(scheme)
  ! generate far-field potential using whole box hierarchy
  do level=TOP_LEVEL,deepest_level
    l_up = level-1
    ! initialise T-contractor and allocate T-matrix
    ! get LHS boxed MM parameters
    call fmm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
    ! get packed RHS MM paras and translated moments
    call fmm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
    call fmm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','keep')
    ! allocate temporary boxed potentials
    if (use_Vff_tmp /= 2) then
      foo = size(LHS%box_paras)
      allocate(Vff_tmp1(lm_dim,foo))
      Vff_tmp1 = zero
      Vff_ptr => Vff_tmp1(:,:)
    else
      foo = size(LHS%box_paras)
      allocate(Vff_tmp2(lm_dim,foo))
      Vff_tmp2 = zero
      Vff_ptr => Vff_tmp2(:,:)
    end if
    ! set T_contractor pointers
    call fmm_set_T_con_ptrs(Vff_ptr,qlm_T)
    ! generate LFF contribution to Vff at this level
    call fmm_init_T_pair_builder(scheme,pair_type)
    call fmm_gen_local_T_pairs(LHS,RHS,pair_type)
    call fmm_close_T_pair_builder()
    ! get RFF contribution from parent level
    select case(use_Vff_tmp)
      case(0)
        use_Vff_tmp = 2
      case(1)
        ! RFF contribution: parents' translated FF ( parent | child )
        call fmm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')
        call fmm_translate_parents_Vff(level,scheme,Vff_tmp2,Vff_tmp1,LHS%box_paras)
        ! deallocate parent's space
        deallocate(Vff_tmp2)
        nullify(Vff_tmp2)
        use_Vff_tmp = 2
      case(2)
        ! RFF contribution: parents' translated FF ( parent | child )
        call fmm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')
        call fmm_translate_parents_Vff(level,scheme,Vff_tmp1,Vff_tmp2,LHS%box_paras)
        ! deallocate parent's space
        deallocate(Vff_tmp1)
        nullify(Vff_tmp1)
        use_Vff_tmp = 1
    end select
  end do
  call fmm_free_T_contractors()

  ! translate boxed potential to raw LHS centres for final contraction
  if (use_Vff_tmp == 1) then  ! Vff_tmp2 holds final boxed potential
    call fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp2,Vff)
  else  ! Vff_tmp1 holds final boxed potential
    call fmm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp1,Vff)
  end if
  if (associated(Vff_tmp1)) deallocate(Vff_tmp1)
  if (associated(Vff_tmp2)) deallocate(Vff_tmp2)

  TTOT = fmm_second()-T0
  call TIMTXT('>>> TIME USED in fmm_get_FMM_Vff',TTOT,LUPRI)

  call fmm_free_W_contractors()

end subroutine fmm_get_FMM_Vff

!-------------------------------------------------------------------------------

end module fmm_Vff_driver
