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

module fmm_global_paras

use Constants, only: Zero, One, Two, Half, Pi
use Definitions, only: wp, iwp, u5, u6

implicit none
private
public :: INTK, REALK, LUPRI, LURD, LUINTM, fmm_stats_printed, Zero, One, Two, Half, Pi
public :: BRFREE_DF, EXTENT_MIN_DF, PACK_RHS_DF, PACK_LHS_DF, WS_MIN, START_LEN, TREE_LENGTH, MAX_AVG_PER_NODE, ZERO_DIST_TOL, &
          TMATM_DF, DISTINCT_T_TOL, ZERO_VECT_TOL, TOP_LEVEL, MAX_LEVEL
public :: GFC_FMM, MD4_FMM, FE_FMM, NEAR_FIELD, FAR_FIELD, DO_NULL, DO_FQ, DO_BQ, DO_NlogN, DO_FMM, LHS_raw_RHS_raw, &
          LHS_box_RHS_box, T_CONTRACTOR_DIRECT, T_CONTRACTOR_BOUNDARY, T_CONTRACTOR_TREE, T_CONTRACTOR_SCALE, &
          T_CONTRACTOR_SCALE_TREE, T_CONTRACTOR_MULTI, T_CONTRACTOR_FULL, W_CONTRACTOR_DIRECT, W_CONTRACTOR_X, W_CONTRACTOR_FAST, &
          W_CONTRACTOR_BOUNDARY, ELECTRONIC_ONLY, NUCLEAR_ONLY, ALL_MOMENTS, USE_RAW_QLM, USE_T_SYM_QLM, NULL_T_BUFFER, &
          NULL_W_BUFFER, TREE_T_BUFFER, TREE_W_BUFFER, SKIP_T_BUFFER, SKIP_W_BUFFER, MULTI_T_BUFFER, SCALE_T_BUFFER, &
          SORT_BY_SCALE, SORT_BY_RHS_MMS, NSPACE
public :: fmm_counters, fmm_planes, LHS_RHS_type, scheme_paras, raw_mm_paras, raw_mm_data, box_mm_paras, box_mm_data, id_node, &
          id_list, gen_mm_paras, T_paras, T_pair_single, T_pair_list, T_pair_batch, old_new, fmm_sh_pairs, fmm_basis, fmm_prim_batch

!------------------------------------------------------------------------------

! We select here the REAL and INTEGER precisions
!integer(kind=iwp), parameter :: INTK = selected_int_kind(9), &
!                                REALK = selected_real_kind(XXX)
integer(kind=iwp), parameter :: INTK = iwp, &
                                REALK = wp

!==============================================================================
! Program-wide global variables |
!===============================!

! Unit numbers for IO
integer(INTK), save :: LUPRI = u6
integer(INTK), save :: LURD = u5
integer(INTK), save :: LUINTM = 77

! Flag for verbose printing of run-time statistics
logical, save :: fmm_stats_printed = .true.

!==============================================================================
! Global parameters |
!===================!

! Branch-free algorithm
logical, parameter :: BRFREE_DF = .true.
! minimum charge separation for qTq evaluation
real(REALK), parameter :: EXTENT_MIN_DF = 0.001_REALK

! Same centre/extent packing flags
logical, parameter :: PACK_RHS_DF = .true.
logical, parameter :: PACK_LHS_DF = .false.

! Parameters for branch and box structure
!integer(INTK), parameter :: WS_MIN = 2
integer(INTK), parameter :: WS_MIN = 1

! Parameters for T-list tree sorting
integer(INTK), parameter :: START_LEN = 8
integer(INTK), parameter :: TREE_LENGTH = 50000
! MAX_AVG_PER_NODE limits the maximal packing ratio.
! there will be never more interactions in the tree than
! tree size times MAX_AVG_PER_NODE.
integer(INTK), parameter :: MAX_AVG_PER_NODE = 15

! Noise buffer when testing for classical interactions
!FIXME: have this passed from clsfmm.h ... used for nothing else!!
real(REALK), parameter :: ZERO_DIST_TOL = 1.0e-14_REALK

! Parameters for contractions and translations
! Number of simultaneous T-matrices when using T_CONTRACTOR_MULTI
integer(INTK), parameter :: TMATM_DF = 25

!FIXME: clarify the use of these tolerance parameters!!
real(REALK), parameter :: DISTINCT_T_TOL = 1.0e-15_REALK
real(REALK), parameter :: ZERO_VECT_TOL = 1.0e-10_REALK

!------------------------------------------------------------------------------

! Parameters to define level depth used in box hierarchy
integer(INTK), parameter :: TOP_LEVEL = 2
integer(INTK), parameter :: MAX_LEVEL = 15

!------------------------------------------------------------------------------
!
! Named parameters for coding use (avoidance of "magic" numbers etc.)
!

! Named parameters for potential or J-matrix request
integer(INTK), parameter :: GFC_FMM = 1, &
                            MD4_FMM = 2, &
                            FE_FMM = 3

! Named parameters for (classical) near-field or far-field phase
integer(INTK), parameter :: NEAR_FIELD = 1, &
                            FAR_FIELD = 2

! Named parameters for different MM methods implemented
integer(INTK), parameter :: DO_NULL = 0, &
                            DO_FQ = 1, &
                            DO_BQ = 3, &
                            DO_NlogN = 4, &
                            DO_FMM = 5

! Named parameters for T-pair interaction types
integer(INTK), parameter :: LHS_raw_RHS_raw = 1, &
                            !LHS_raw_RHS_box = 2, &
                            !LHS_box_RHS_raw = 3, &
                            LHS_box_RHS_box = 4

! Named parameters for different contractors implemented
integer(INTK), parameter :: T_CONTRACTOR_DIRECT = 100, &
                            T_CONTRACTOR_BOUNDARY = 101, &
                            T_CONTRACTOR_TREE = 102, &
                            T_CONTRACTOR_SCALE = 104, &
                            T_CONTRACTOR_SCALE_TREE = 1041, &
                            T_CONTRACTOR_MULTI = 105, &
                            T_CONTRACTOR_FULL = 106, &
                            W_CONTRACTOR_DIRECT = 206, &
                            W_CONTRACTOR_X = 207, &
                            W_CONTRACTOR_FAST = 208, &
                            W_CONTRACTOR_BOUNDARY = 209

! Named parameters for moments used in contraction
integer(INTK), parameter :: ELECTRONIC_ONLY = 1, &
                            NUCLEAR_ONLY = 2, &
                            ALL_MOMENTS = 3

! Named parameters for auxiliary moments
integer(INTK), parameter :: USE_RAW_QLM = 0, &
                            USE_T_SYM_QLM = 1

! Named parameters for lists available
integer(INTK), parameter :: NULL_T_BUFFER = 1, &
                            NULL_W_BUFFER = 2, &
                            TREE_T_BUFFER = 3, &
                            TREE_W_BUFFER = 4, &
                            SKIP_T_BUFFER = 5, &
                            SKIP_W_BUFFER = 6, &
                            MULTI_T_BUFFER = 7, &
                            SCALE_T_BUFFER = 8

! Sort orders expected by interaction evaluators.
integer(INTK), parameter :: SORT_BY_SCALE = 1, &
                            SORT_BY_RHS_MMS = 2

! Named parameters for memory management
character(len=7), parameter :: NSPACE(5) = ['raw_qlm', &
                                            'raw_Vff', &
                                            'box_Vff', &
                                            'Vff_tmp', &
                                            'box_qlm']

!===============================================================================
! Derived types |
!================

! Structure for counting raw multipole moments by type
type fmm_counters
  integer(INTK) :: tot, nuc, elec
end type fmm_counters

!------------------------------------------------------------------------------

! Structure for defining boundary planes (assuming // to axes)
type fmm_planes
  real(REALK) :: xmin, xmax, ymin, ymax, zmin, zmax
end type fmm_planes

!------------------------------------------------------------------------------

type LHS_RHS_type
  integer(INTK) :: LHS, RHS
end type LHS_RHS_type

!------------------------------------------------------------------------------

! Complete prescription of MM execution held in "scheme"
type T_contract_schm
  integer(INTK) :: NF_ID, FF_ID
  integer(INTK) :: NF_T_buffer, FF_T_buffer
  integer(INTK) :: NF_sort_para, FF_sort_para
  integer(INTK) :: LHS_mm_type, RHS_mm_type
end type T_contract_schm

type W_contract_schm
  integer(INTK) :: ID         ! general translator
  integer(INTK) :: BR_ID      ! box-to-raw translator
  integer(INTK) :: W_buffer
  integer(INTK) :: sort_para
end type W_contract_schm

type scheme_paras
  integer(INTK)         :: job_type
  logical               :: include_near_field
  integer(INTK)         :: algorithm, phase
  integer(INTK)         :: NF_T_searcher, FF_T_searcher
  type(T_contract_schm) :: T_con
  type(W_contract_schm) :: W_con
  logical               :: branch_free
  integer(INTK)         :: raw_lmax, trans_lmax
  integer(INTK)         :: LHS_mm_range, RHS_mm_range
  logical               :: pack_LHS, pack_RHS
  logical               :: LHS_dens, RHS_dens
  real(REALK)           :: extent_min
  real(REALK)           :: grain
  real(REALK)           :: dens_screen_thr
  integer(INTK)         :: FEdim
  ! Order of Lagrange interpolating polynomials for finite elements
  integer(INTK)         :: lipn

end type scheme_paras

!------------------------------------------------------------------------------

! Raw (not boxed) multipole moment data structure
type raw_mm_paras
  real(REALK)   :: cntr(3), ext
  integer(INTK) :: batch            ! batch id of unique centres/extents
  integer(INTK) :: id               ! map to raw moments array
  integer(INTK) :: map_up           ! map to boxed multipole moment
  integer(INTK) :: box(3), bra
  real(REALK)   :: box_cntr(3)
end type raw_mm_paras

type J_index_type
  integer(INTK) :: i_indx, j_indx   ! map to J_matrix elements
end type J_index_type

type raw_mm_data
  type(raw_mm_paras), pointer :: paras(:)
  real(REALK), pointer        :: dens(:)
  real(REALK), pointer        :: qlm(:,:)
  real(REALK), pointer        :: qlm_T(:,:)
  real(REALK), pointer        :: qlm_W(:,:)
  type(J_index_type), pointer :: J_indices(:)
  type(id_list), pointer      :: batch_map(:) ! maps raw paras in one batch
end type raw_mm_data

!------------------------------------------------------------------------------

! Boxed multipole moment data structure
type box_mm_paras
  integer(INTK) :: box(3)
  real(REALK)   :: cntr(3)
  integer(INTK) :: bra, level
  integer(INTK) :: map_up            ! map to moment at next level
  real(REALK)   :: cntr_up(3)        ! centre of parent box
  integer(INTK) :: id                ! map to moment at this level
end type box_mm_paras
type box_mm_data
  type(box_mm_paras), pointer :: LHS_paras(:)
  type(box_mm_paras), pointer :: RHS_paras(:)
  ! Just RHS moments (only RHS is built up in box hierarchy)
  real(REALK), pointer        :: qlm_T(:,:)
  real(REALK), pointer        :: qlm_W(:,:)
end type box_mm_data

!------------------------------------------------------------------------------

! Map structure to link packed parameters with all raw moment members
type id_node
  integer(INTK)          :: id
  integer(INTK)          :: i_indx, j_indx
  type(id_node), pointer :: next
end type id_node

type id_list
  integer(INTK)          :: occ
  type(id_node), pointer :: head
end type id_list

!------------------------------------------------------------------------------

! Generalised multipole moment parameter structure
type gen_mm_paras
  type(raw_mm_paras), pointer :: raw_paras(:)
  type(box_mm_paras), pointer :: box_paras(:)
end type gen_mm_paras

!------------------------------------------------------------------------------

type T_paras
  integer(INTK) :: LHS_lmax, LHS_id
  integer(INTK) :: RHS_lmax, RHS_id
  integer(INTK) :: weight
  real(REALK)   :: ratio
end type T_paras

type T_pair_single
  type(T_paras) :: paras
  real(REALK)   :: r_ab(3)
  integer(INTK) :: lmax, lm_max
  ! Used only for W_pairs when translating to distinguish qlm and Vff modes
  character     :: N_or_T
end type T_pair_single

type T_pair_list
  type(T_paras), pointer :: paras(:)
  real(REALK)            :: r_ab(3)
  integer(INTK)          :: lmax, lm_max
  integer(INTK)          :: LHS_lmax, RHS_lmax
  ! Used only for W_pairs when translating to distinguish qlm and Vff modes
  character              :: N_or_T
end type T_pair_list

type T_pair_batch
  type(T_pair_single), pointer :: items(:)
  integer(INTK)                :: ndim
end type T_pair_batch

! Index type for OLD and NEW W_translations
type old_new
  integer(INTK) :: old, new
end type old_new

!------------------------------------------------------------------------------

! Linked-list shell-pair data structure
!type fmm_shell_pair_node
!  integer(INTK)                      :: I, J
!  integer(INTK)                      :: box(3)
!  real(REALK)                        :: extent
!  real(REALK)                        :: centre(3)
!  type(fmm_shell_pair_node), pointer :: next
!end type fmm_shell_pair_node

!-----------------------------------------------------------------------------

! O(N) array of non-vanishing shell paris
type fmm_sh_pairs
  integer(INTK) :: I, J
  real(REALK)   :: extent
  real(REALK)   :: centre(3)
end type fmm_sh_pairs

!------------------------------------------------------------------------------

! Basis set information from host program
type fmm_basis
  integer(INTK)          :: nbas
  integer(INTK)          :: nshells
  integer(INTK)          :: maxsgm2
  integer(INTK)          :: maxangl
  integer(INTK), pointer :: KAtom(:)
  integer(INTK), pointer :: KType(:)
  integer(INTK), pointer :: KStart(:)
  integer(INTK), pointer :: KontG(:)
  integer(INTK), pointer :: KLoc_Car(:)
  integer(INTK), pointer :: LtuvMin_Car(:)
  integer(INTK), pointer :: LtuvMax_Car(:)
  integer(INTK), pointer :: Lt(:)
  integer(INTK), pointer :: Lu(:)
  integer(INTK), pointer :: Lv(:)
  real(REALK), pointer   :: Centr(:,:)
  real(REALK), pointer   :: Expnt(:)
  real(REALK), pointer   :: CCoef(:)
end type fmm_basis

!------------------------------------------------------------------------------

! Batch of primitive data for one shell-pair
type fmm_prim_batch
  real(REALK) :: P(3), PA(3), PB(3), PC(3)
  real(REALK) :: lo(3), hi(3)
  real(REALK) :: ExpntP
  real(REALK) :: ExpPHalf
  real(REALK) :: PreFactAB
  real(REALK) :: CCoefAB
end type fmm_prim_batch

!------------------------------------------------------------------------------

end module fmm_global_paras
