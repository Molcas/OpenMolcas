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
MODULE fmm_global_paras

   IMPLICIT NONE
   PUBLIC

!------------------------------------------------------------------------------

   ! We select here the REAL and INTEGER precisions
   !   INTEGER, PARAMETER :: INTK  = SELECTED_INT_KIND(9),           &
   !                         REALK = SELECTED_REAL_KIND(XXX)
   ! Use standard integers and REAL*8 double precision floats
   INTEGER, PARAMETER :: INTK  = KIND(1),                         &
                         REALK = KIND(1D0)

!==============================================================================
! Program-wide global variables |
!===============================!

    ! Unit numbers for IO
   INTEGER(INTK), SAVE :: LUPRI = 6
   INTEGER(INTK), SAVE :: LUINTM = 77
    ! Interface file name (contains moments and co-ordinate data)
   CHARACTER(LEN=7),  SAVE :: INPUT_FILE = 'MM_DATA'
   CHARACTER(LEN=14), SAVE :: INPUT_FILE_HEADER = 'MM_DATA_HEADER'

    ! Flag for verbose printing of run-time statistics
   LOGICAL, SAVE :: fmm_stats_printed = .TRUE.

!==============================================================================
! Global parameters |
!===================!

   REAL(REALK), PARAMETER :: zero   = 0.0_REALK,                        &
                             one    = 1.0_REALK,                        &
                             two    = 2.0_REALK,                        &
                             half   = 0.5_REALK,                        &
                             PI     = 3.14159265358979323846264_REALK,  &
                             ROOTPI = 1.77245385090552_REALK


!------------------------------------------------------------------------------

   ! Branch-free algorithm
   LOGICAL,     PARAMETER :: BRFREE_DF = .TRUE.
   ! minimum charge separation for qTq evaluation
   REAL(REALK), PARAMETER :: EXTENT_MIN_DF = 0.001_REALK

   ! Same centre/extent packing flags
   LOGICAL, PARAMETER :: PACK_RHS_DF = .TRUE.
   LOGICAL, PARAMETER :: PACK_LHS_DF = .FALSE.

    ! Parameters for branch and box structure
!   INTEGER(INTK), PARAMETER :: WS_MIN = 2
   INTEGER(INTK), PARAMETER :: WS_MIN = 1
   LOGICAL,       PARAMETER :: JOIN_BRANCHES = .TRUE.

    ! Parameters for T-list tree sorting
   INTEGER(INTK), PARAMETER :: START_LEN   = 8
   INTEGER(INTK), PARAMETER :: TREE_LENGTH = 50000
   ! MAX_AVG_PER_NODE limits the maximal packing ratio.
   ! there will be never more interactions in the tree than
   ! tree size times MAX_AVG_PER_NODE.
   INTEGER(INTK), PARAMETER :: MAX_AVG_PER_NODE = 15

    ! Noise buffer when testing for classical interactions
 !FIXME: have this passed from clsfmm.h ... used for nothing else!!
   REAL(REALK), PARAMETER :: ZERO_DIST_TOL = 1e-14_REALK

    ! Parameters for contractions and translations
   ! Number of simultaneous T-matrices when using T_CONTRACTOR_MULTI
   INTEGER(INTK), PARAMETER :: TMATM_DF = 25

 !FIXME: clarify the use of these tolerance parameters!!
   REAL(REALK), PARAMETER :: DISTINCT_T_TOL  = 1e-15_REALK
   REAL(REALK), PARAMETER :: EXTENT_TEST_TOL = 1e-15_REALK
   REAL(REALK), PARAMETER :: ZERO_VECT_TOL   = 1e-10_REALK

!------------------------------------------------------------------------------

    ! Parameters to define level depth used in box hierarchy
   INTEGER(INTK), PARAMETER :: TOP_LEVEL = 2
   INTEGER(INTK), PARAMETER :: MAX_LEVEL = 15

!------------------------------------------------------------------------------
!
! Named parameters for coding use (avoidance of "magic" numbers etc.)
!
    ! Named parameters for potential or J-matrix request
   INTEGER(INTK), PARAMETER :: GFC_FMM = 1,  &
                               MD4_FMM = 2,  &
                               FE_FMM = 3

    ! Named parameters for (classical) near-field or far-field phase
   INTEGER(INTK), PARAMETER :: NEAR_FIELD = 1 ,    &
                               FAR_FIELD  = 2

    ! Named parameters for different MM methods implemented
   INTEGER(INTK), PARAMETER :: DO_NULL  = 0 ,     &
                               DO_FQ    = 1 ,     &
                               DO_BQ    = 3 ,     &
                               DO_NlogN = 4 ,     &
                               DO_FMM   = 5

    ! Named parameters for T-pair interaction types
   INTEGER(INTK), PARAMETER :: LHS_raw_RHS_raw = 1,     &
                               LHS_raw_RHS_box = 2,     &
                               LHS_box_RHS_raw = 3,     &
                               LHS_box_RHS_box = 4

    ! Named parameters for T-pair interaction search algorithms
   INTEGER(INTK), PARAMETER :: FQ_LOOPS    = 1,     &
                               GRID_SEARCH = 2,     &
                               TREE_SEARCH = 3

    ! Named parameters for different contractors implemented
   INTEGER(INTK), PARAMETER :: T_CONTRACTOR_DIRECT = 100,        &
                               T_CONTRACTOR_BOUNDARY = 101,      &
                               T_CONTRACTOR_TREE = 102,          &
                               T_CONTRACTOR_SCALE = 104,         &
                               T_CONTRACTOR_SCALE_TREE = 1041,   &
                               T_CONTRACTOR_MULTI = 105,         &
                               T_CONTRACTOR_FULL = 106,          &
                               W_CONTRACTOR_DIRECT = 206,        &
                               W_CONTRACTOR_X = 207,             &
                               W_CONTRACTOR_FAST = 208,          &
                               W_CONTRACTOR_BOUNDARY = 209

    ! Named parameters for different J-builders (contractors) implemented
   INTEGER(INTK), PARAMETER :: J_BUILDER_DIRECT = 1,        &
                               J_BUILDER_TREE = 2

    ! Named parameters for moments used in contraction
   INTEGER(INTK), PARAMETER :: ELECTRONIC_ONLY = 1,        &
                               NUCLEAR_ONLY    = 2,        &
                               ALL_MOMENTS     = 3

    ! Named parameters for T-pair interaction tests
   INTEGER(INTK), PARAMETER :: TEST_EXT     = 1,        &
                               TEST_R_IJ    = 2,        &
                               TEST_NN_EXT  = 3,        &
                               TEST_NN_R_IJ = 4,        &
                               TEST_FF      = 5,        &
                               TEST_LFF     = 6

    ! Named parameters for auxiliary moments
   INTEGER(INTK), PARAMETER :: USE_RAW_QLM   = 0 ,        &
                               USE_T_SYM_QLM = 1

    ! Named parameters for lists available
   INTEGER(INTK), PARAMETER :: NULL_T_BUFFER  = 1 ,        &
                               NULL_W_BUFFER  = 2 ,        &
                               TREE_T_BUFFER  = 3 ,        &
                               TREE_W_BUFFER  = 4 ,        &
                               SKIP_T_BUFFER  = 5 ,        &
                               SKIP_W_BUFFER  = 6 ,        &
                               MULTI_T_BUFFER = 7 ,        &
                               SCALE_T_BUFFER = 8

    ! Sort orders expected by interaction evaluators.
   INTEGER(INTK), PARAMETER :: NO_SORT         = 0 ,       &
                               SORT_BY_SCALE   = 1,        &
                               SORT_BY_RHS_MMS = 2

    ! Named parameters for memory management
   INTEGER(INTK), PARAMETER :: NMEMDIVS = 5
   CHARACTER(LEN=7),  PARAMETER :: NSPACE(NMEMDIVS) = (/ 'raw_qlm',   &
                                                     'raw_Vff',   &
                                                     'box_Vff',   &
                                                     'Vff_tmp',   &
                                                     'box_qlm' /)
   INTEGER(INTK), PARAMETER :: MEM_RAW_QLM = 1,        &
                               MEM_RAW_VFF = 2,        &
                               MEM_BOX_VFF = 3,        &
                               MEM_VFF_TMP = 4,        &
                               MEM_BOX_QLM = 5

!===============================================================================
! Derived types |
!================

    ! Structure for counting raw multipole moments by type
   TYPE fmm_counters
      INTEGER(INTK) :: tot, nuc, elec
   END TYPE fmm_counters

!------------------------------------------------------------------------------

    ! Structure for defining boundary planes (assuming // to axes)
   TYPE fmm_planes
      REAL(REALK) :: xmin,xmax, ymin,ymax, zmin,zmax
   END TYPE fmm_planes

!------------------------------------------------------------------------------

   TYPE LHS_RHS_type
      INTEGER(INTK) :: LHS, RHS
   END TYPE LHS_RHS_type

!------------------------------------------------------------------------------

   TYPE fmm_range
      INTEGER(INTK) :: hi, lo, tot
   END TYPE fmm_range

!------------------------------------------------------------------------------

    ! Complete prescription of MM execution held in "scheme"
   TYPE T_contract_schm
      INTEGER(INTK) :: NF_ID, FF_ID
      INTEGER(INTK) :: NF_T_buffer, FF_T_buffer
      INTEGER(INTK) :: NF_sort_para, FF_sort_para
      INTEGER(INTK) :: LHS_mm_type, RHS_mm_type
   END TYPE T_contract_schm

   TYPE W_contract_schm
      INTEGER(INTK) :: ID         ! general translator
      INTEGER(INTK) :: BR_ID      ! box-to-raw translator
      INTEGER(INTK) :: W_buffer
      INTEGER(INTK) :: sort_para
   END TYPE W_contract_schm

   TYPE scheme_paras
      INTEGER(INTK)         :: job_type
      LOGICAL               :: include_near_field
      INTEGER(INTK)         :: algorithm, phase
      INTEGER               :: NF_T_searcher, FF_T_searcher
      TYPE(T_contract_schm) :: T_con
      TYPE(W_contract_schm) :: W_con
      LOGICAL               :: branch_free
      INTEGER(INTK)         :: raw_lmax, trans_lmax
      INTEGER(INTK)         :: LHS_mm_range, RHS_mm_range
      LOGICAL               :: pack_LHS, pack_RHS
      LOGICAL               :: LHS_dens, RHS_dens
      REAL(REALK)           :: extent_min
      REAL(REALK)           :: grain
      REAL(REALK)           :: dens_screen_thr
      INTEGER(INTK)         :: FEdim
      ! Order of Lagrange interpolating polynomials for finite elements
      INTEGER(INTK)         :: lipn

   END TYPE scheme_paras

!------------------------------------------------------------------------------

    ! Raw (not boxed) multipole moment data structure
   TYPE raw_mm_paras
      REAL(REALK)   :: cntr(3), ext
      INTEGER(INTK) :: batch            ! batch id of unique centres/extents
      INTEGER(INTK) :: id               ! map to raw moments array
      INTEGER(INTK) :: map_up           ! map to boxed multipole moment
      INTEGER(INTK) :: box(3), bra
      REAL(REALK)   :: box_cntr(3)
   END TYPE raw_mm_paras

   TYPE J_index_type
      INTEGER(INTK) :: i_indx, j_indx   ! map to J_matrix elements
   END TYPE J_index_type

   TYPE raw_mm_data
      TYPE(raw_mm_paras), POINTER :: paras(:)
      REAL(REALK),        POINTER :: dens(:)
      REAL(REALK),        POINTER :: qlm(:,:)
      REAL(REALK),        POINTER :: qlm_T(:,:)
      REAL(REALK),        POINTER :: qlm_W(:,:)
      TYPE(J_index_type), POINTER :: J_indices(:)
      TYPE(id_list),      POINTER :: batch_map(:) ! maps raw paras in one batch
   END TYPE raw_mm_data

!------------------------------------------------------------------------------

    ! Boxed multipole moment data structure
   TYPE box_mm_paras
      INTEGER(INTK) :: box(3)
      REAL(REALK)   :: cntr(3)
      INTEGER(INTK) :: bra, level
      INTEGER(INTK) :: map_up            ! map to moment at next level
      REAL(REALK)   :: cntr_up(3)        ! centre of parent box
      INTEGER(INTK) :: id                ! map to moment at this level
   END TYPE box_mm_paras
   TYPE box_mm_data
      TYPE(box_mm_paras), POINTER :: LHS_paras(:)
      TYPE(box_mm_paras), POINTER :: RHS_paras(:)
      ! Just RHS moments (only RHS is built up in box hierarchy)
      REAL(REALK),        POINTER :: qlm_T(:,:)
      REAL(REALK),        POINTER :: qlm_W(:,:)
   END TYPE box_mm_data

!------------------------------------------------------------------------------

    ! Map structure to link packed parameters with all raw moment members
   TYPE id_node
      INTEGER(INTK) :: id
      INTEGER(INTK) :: i_indx, j_indx
      TYPE(id_node), POINTER :: next
   END TYPE id_node

   TYPE id_list
      INTEGER(INTK) :: occ
      TYPE(id_node), POINTER :: head
   END TYPE id_list

!------------------------------------------------------------------------------

    ! Generalised multipole moment parameter structure
   TYPE gen_mm_paras
      TYPE(raw_mm_paras), POINTER :: raw_paras(:)
      TYPE(box_mm_paras), POINTER :: box_paras(:)
   END TYPE gen_mm_paras

!------------------------------------------------------------------------------

   TYPE T_paras
      INTEGER(INTK) :: LHS_lmax, LHS_id
      INTEGER(INTK) :: RHS_lmax, RHS_id
      INTEGER(INTK) :: weight
      REAL(REALK)   :: ratio
   END TYPE T_paras

   TYPE T_pair_single
      TYPE(T_paras) :: paras
      REAL(REALK)   :: r_ab(3)
      INTEGER(INTK) :: lmax, lm_max
      ! Used only for W_pairs when translating to distinguish qlm and Vff modes
      CHARACTER(LEN=1)  :: N_or_T
   END TYPE T_pair_single

   TYPE T_pair_list
      TYPE(T_paras), POINTER :: paras(:)
      REAL(REALK)   :: r_ab(3)
      INTEGER(INTK) :: lmax, lm_max
      INTEGER(INTK) :: LHS_lmax, RHS_lmax
      ! Used only for W_pairs when translating to distinguish qlm and Vff modes
      CHARACTER(LEN=1)  :: N_or_T
   END TYPE T_pair_list

   TYPE T_pair_batch
      TYPE(T_pair_single), POINTER :: items(:)
      INTEGER(INTK) :: ndim
   END TYPE T_pair_batch

    ! Index type for OLD and NEW W_translations
   TYPE old_new
      INTEGER(INTK) :: old, new
   END TYPE old_new

!------------------------------------------------------------------------------
!
!    ! Linked-list shell-pair data structure
!   TYPE fmm_shell_pair_node
!      INTEGER(INTK) :: I,J
!      INTEGER(INTK) :: box(3)
!      REAL(REALK)   :: extent
!      REAL(REALK)   :: centre(3)
!      TYPE(fmm_shell_pair_node), POINTER :: next
!   END TYPE fmm_shell_pair_node
!
!!-----------------------------------------------------------------------------
!
!    ! Shell-pairs are collected into boxes
!   TYPE fmm_boxed_sh_pairs
!      TYPE(box_mm_paras) :: paras
!      TYPE(fmm_shell_pair_node), POINTER :: sh_pairs(:)
!   END TYPE fmm_boxed_sh_pairs
!
!------------------------------------------------------------------------------

    ! O(N) array of non-vanishing shell paris
   TYPE fmm_sh_pairs
      INTEGER(INTK) :: I, J
      REAL(REALK)   :: extent
      REAL(REALK)   :: centre(3)
   END TYPE fmm_sh_pairs

!------------------------------------------------------------------------------

    ! Basis set information from host program
   TYPE fmm_basis
      INTEGER(INTK) :: nbas
      INTEGER(INTK) :: nshells
      INTEGER(INTK) :: maxsgm2
      INTEGER(INTK) :: maxangl
      INTEGER(INTK), POINTER :: KAtom(:)
      INTEGER(INTK), POINTER :: KType(:)
      INTEGER(INTK), POINTER :: KStart(:)
      INTEGER(INTK), POINTER :: KontG(:)
      INTEGER(INTK), POINTER :: KLoc_Car(:)
      INTEGER(INTK), POINTER :: LtuvMin_Car(:)
      INTEGER(INTK), POINTER :: LtuvMax_Car(:)
      INTEGER(INTK), POINTER :: Lt(:)
      INTEGER(INTK), POINTER :: Lu(:)
      INTEGER(INTK), POINTER :: Lv(:)
      REAL(REALK), POINTER :: Centr(:,:)
      REAL(REALK), POINTER :: Expnt(:)
      REAL(REALK), POINTER :: CCoef(:)
   END TYPE fmm_basis

!------------------------------------------------------------------------------

    ! Batch of primitive data for one shell-pair
   TYPE fmm_prim_batch
      REAL(REALK) :: P(3), PA(3), PB(3), PC(3)
      REAL(REALK) :: lo(3), hi(3)
      REAL(REALK) :: ExpntP
      REAL(REALK) :: ExpPHalf
      REAL(REALK) :: PreFactAB
      REAL(REALK) :: CCoefAB
   END TYPE fmm_prim_batch

!------------------------------------------------------------------------------

   TYPE fmm_box_id_type
      INTEGER(INTK) :: id
      REAL(REALK)   :: centre(3)
   END TYPE fmm_box_id_type

!------------------------------------------------------------------------------

END MODULE fmm_global_paras


