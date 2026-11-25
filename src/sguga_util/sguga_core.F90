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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
! The core of SGUGA: types, states, shared utilities (coupling-block
! lookups, CSF encoding, level maps, ASTA segment values), SG_Free.
! Default public: the other SGUGA_* modules use its entities; the
! interface of the library is the module SGUGA (sguga.F90).
module SGUGA_Core

use Molcas, only: MxLev, MxSym
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6
use CSF_Graph, only: CSFG_ATab, CSFG_BTab, CSFG_Count, CSFG_Graph, CSFG_LTab, CSFG_Make, CSFG_Mk_IOW, CSFG_MkCOT, CSFG_nPack, &
                     CSFG_Walks

#include "intent.fh"

implicit none

integer(kind=iwp) :: iq, NCP_Max

! Slater Determinant reordering to fit the Lucia machinary
! This might actually get redundant once the Lucia codes have been removed!
type :: SDReoInfo
  integer(kind=iwp), allocatable :: List(:)
end type SDReoInfo

! Configuration-to-Determinant transformation informtion

type CDCInfo
  integer(kind=iwp) :: nDet = 0
  integer(kind=iwp) :: nCSF = 0
  ! Lucia spin-combination dimension
  integer(kind=iwp) :: nDetSC = 0

  integer(kind=iwp) :: nClosed = 0
  integer(kind=iwp) :: nClosedOCP = 0

  integer(kind=iwp), allocatable :: DOC(:,:)
  integer(kind=iwp), allocatable :: SOC(:,:)
  integer(kind=iwp), allocatable :: SOCSym(:)
  integer(kind=iwp), allocatable :: nOCSym(:)
  integer(kind=iwp), allocatable :: OCDOC(:,:)
  integer(kind=iwp), allocatable :: OCSOC(:,:)
  integer(kind=iwp), allocatable :: OCReverse(:,:,:)

  integer(kind=iwp), allocatable :: nOpenOCP(:)
  integer(kind=iwp), allocatable :: nOCP(:)
  integer(kind=iwp), allocatable :: DetSpin(:,:)
  integer(kind=iwp), allocatable :: Twos(:,:)

  ! Full determinant representation
  real(kind=wp), allocatable :: CDC(:,:)

  ! Lucia spin-combination representation
  real(kind=wp), allocatable :: DTOC(:,:)

  ! Pairing information
  integer(kind=iwp), allocatable :: Pair(:)

end type CDCInfo

! Split-Graph descriptor, sizes, addresses...
!
! Ism(:): Symmetry label of each graph level.
! DRT: Final compact distinct-row table with columns CSFG_LTab,CSFG_NTab,CSFG_ATab,CSFG_BTab,CSFG_CTab.
! DRT0: Unrestricted DRT before pruning.
! Down: Final downchar arc table; child vetrex reached by each of the four step types
! Down0: unrestricted Down table
! Up: Reverse connectivity; for a lower vertex and a step typ, gives the matching upper vertex
! Ver: temporary mask and old->new renumbering map used during pruning
! MAW: modified split-graph are weights used in split-graph numbering. A hybrid of the DAW and the RAW tables.
! LTV: level-to-first-vertex pointer table
! DAW: Direct arc weights; cumulative numbering of downward continuation
! RAW: reverse direct arc weights; cumulative numbering of upward continuation
! L2ACT: map from graph level to active orbital index
! LEVEL: external ordering of graphs levels
! nRas: number of orbitals in each symmetry and RAS/GAS partionining
! nRasEl: Electron limit for each RAS/GAS partioning
! The graph (CSFG_Graph of CSF_Graph) and the parts of SGUGA (CSF/SD codec)
type, extends(CSFG_Graph) :: SGStruct
  type(CDCInfo), allocatable :: CDCInfo(:,:)
  integer(kind=iwp), allocatable :: OC(:)
  type(SDReoInfo), allocatable :: SDREO(:)
end type SGStruct

! the DRT column indices CSFG_LTab, CSFG_NTab, CSFG_ATab, CSFG_BTab, CSFG_CTab come from CSF_Graph

! CI Structures, addresses,..
!
! nIpWlk: number of integer per walk to store the packed stepvector of a half walk
! NOW(:,:,:): number of upper (1) and lower (2) half-walks by symmetry and midvertex
! IOW(:,:,:): offsets of the NOW blocks in packed walk storage
! NCSF(:): Total number of CSFs per total symmetry
! NOCSF(:,:,:): number of CSFs for a given upper symmetry, midvertex, and total symmetry.
! IOCSF(:,:,:): Offsets corresponding to NOCSF blocks
! ICASE(:): packed step vectors for all upper and lower half-walks
! IVR(:,:): partner index one the same level for delta(b)=-1/+1 segments tops.
! ISGM(:,:): segment connectivity: lower-left destination vertex for each valid segment, zero otherwise.
! VSGM(:,): numerical value of each valid segment-
! The walks (CSFG_Walks of CSF_Graph) and the segment tables of SGUGA
type, extends(CSFG_Walks) :: CIStruct
  integer(kind=iwp), allocatable :: IVR(:,:), ISGM(:,:)
  real(kind=wp), allocatable :: VSGM(:,:)
end type CIStruct

! Excitation operators, coupling coefficients,...
!
! NOCP(:,:,:): number of compressed coupling coefficients in each operator/symmetry/midvertex block
! IOCP(:,:,:): offset of each NOCP block in ICoup
! ICoup(:,:): Compressed coupling tuples: left half-walk label, righ half-walk label. value-pool index
! MVL(:,:) : Left partner midvertex for delta(b)=-1/+1
! MVR(:,:) : Right partner midvertex for delta(b)=-1/+1
! USGN(:,:): map from upper reverese-direct-arc-weight sums to upper walk numbers
! LSGN(:,:): map from lower direct-arc-weight sums to lower walk numbers
! VTAB(:): pool of distinct numerical coupling values.
! IOrd(k,:): index of the coupling tuples by walk (Mk_Coup_Index): within each
!            NOCP block (offset IOCP), IOrd(k,IOCP+i), i=1..NOCP, is the local
!            index of the tuple with the i-th smallest ICoup(k,:) (k=1,2).
! SGTMP(:): intermediates of SG_Epq_Psi for the operators crossing the
!           midlevel; with Reuse_SGTMP they are kept for the next call with
!           the same upper-half orbital (Save_U), the same symmetry of the
!           lower-half orbital (Save_LSym) and the same CI vector (Save_CISym,
!           Save_nCSF, and the tag Save_Tag the caller gave the vector, or for
!           untagged calls (Save_Tag = 0) the fingerprints Save_ID); Save_Valid
!           tells whether SGTMP holds the intermediates of that key,
!           Save_Owner which routine wrote them (1: SG_Epq_Psi,
!           2: SG_EpqEqp_Diag; different layouts).
type EXStruct
  integer(kind=iwp) :: MxEO, nICoup
  integer(kind=iwp), allocatable :: NOCP(:,:,:), IOCP(:,:,:), ICoup(:,:), MVL(:,:), MVR(:,:), USGN(:,:), LSGN(:,:), IOrd(:,:)
  real(kind=wp), allocatable :: VTab(:), SGTMP(:)
  logical(kind=iwp) :: Reuse_SGTMP
  logical(kind=iwp) :: Save_Valid = .false., Save_Swap = .false.
  integer(kind=iwp) :: Save_Owner = 0, Save_U = 0, Save_LSym = 0, Save_CISym = 0, Save_nCSF = 0, Save_Tag = 0
  ! the symmetry for which USGN and LSGN (MkSgNum) were built (0: none)
  integer(kind=iwp) :: Num_Sym = 0
  real(kind=wp) :: Save_ID(2) = 0.0_wp
  integer(kind=iwp), allocatable :: I1list(:), I2list(:)
  real(kind=wp), allocatable :: XList(:)
end type EXStruct

! Configurations (orbital occupations) of the CSFs of one symmetry, built
! from the graph alone (SG_Mk_Cnf_List), in a canonical order: increasing
! packed occupation Key (2 bits per level, (bit_size-2)/2 levels per word,
! the first level in the lowest bits of word 1, compared from the last word
! down).
!   nCnf      : number of configurations
!   iOff(k)   : the CSFs of configuration k are iCSF(iOff(k)+1:iOff(k+1)),
!               split-graph CSF indices in increasing order
!   Key(:,k)  : packed occupation of configuration k (SG_Cnf_Occ unpacks it)
type SGCnfList
  integer(kind=iwp) :: LSM = 0, nCnf = 0, nWord = 0, nLev = 0
  integer(kind=iwp), allocatable :: iOff(:), iCSF(:), Key(:,:)
end type SGCnfList

type TRStruct
  ! Metadata
  integer(kind=iwp) :: nClass = 0   ! number of topological classes
  integer(kind=iwp) :: nOpenBand = 0   ! number of open-family bands
  integer(kind=iwp) :: nTrans = 0   ! total number of transitions

  ! Bucket structure:
  !   NTR(IVLT,ICLASS) = number of transitions from source vertex IVLT
  !                      given input class ICLASS
  !   ITR0(IVLT,ICLASS) = offset into flat transition arrays
  integer(kind=iwp), allocatable :: NTR(:,:)   ! (nVert,0:nClass-1)
  integer(kind=iwp), allocatable :: ITR0(:,:)  ! (nVert,0:nClass-1)

  ! Flat transition arrays, indexed by ITR = 1..nTrans
  integer(kind=iwp), allocatable :: ISGT(:)    ! segment type
  integer(kind=iwp), allocatable :: IVLT(:)    ! source left-upper vertex
  integer(kind=iwp), allocatable :: IVLB(:)    ! destination left-lower vertex

  integer(kind=iwp), allocatable :: ITOP(:)    ! required input class
  integer(kind=iwp), allocatable :: IBOT(:)    ! resulting output class

  integer(kind=iwp), allocatable :: ICL(:)     ! left step code
  integer(kind=iwp), allocatable :: ICR(:)     ! right step code

  ! Right-partner handling:
  ! IPRT = compact partner-slot index, 0 if no distinct right partner needed
  integer(kind=iwp), allocatable :: IPRT(:)

  ! Symmetry handling:
  ! ISYM can store the local symmetry factor label (often 1 or SGS%ISm(LEV))
  integer(kind=iwp), allocatable :: ISYM(:)

  ! Open-band index for coupling-family numbering:
  ! 0 if not an open-family transition
  integer(kind=iwp), allocatable :: IOBAND(:)

  ! Optional flags for transition type
  integer(kind=iwp), allocatable :: IFLAG(:)

  ! Optional precomputed MAW increments
  integer(kind=iwp), allocatable :: MAWL(:)
  integer(kind=iwp), allocatable :: MAWR(:)

  ! Segment value
  real(kind=wp), allocatable :: VSEG(:)
end type TRStruct

! This lists nSegTot different types of segments, i=1,...,nSegTot
!  For the operator E_ij
!  1- 4: segments of the head walk from the loop head to the graph head
!  5- 8: head segments
!  9-13: intermediate segments for the case delta(b)=-1
! 14-18: intermediate segments for the case delta(b)=+1
! 19-22: tails segments
! 23-26: segments of the tail walk from the loop tail to the graph tail
! For the operator E_ii
! 27-30: upper walk diagonal operator
! 31-34: lower walk diagonal operator

! Segment values according to ASTA.

! Vector descriptions:
! IC1(i) and IC2(i): each segment, i, is described by the pair of step vector (IC1(i),IC2(I)), where IC1(i) is the step vector of
! the bra CSF and iC2(i) is the step vector of the ket CSF.
!  ISVC(i): the index ISVC(i), tells which formula to use to compute the segment value of the associated segment
! ITVPT(i): denotes the class of the top two vertices, the upper boundary state.
! IBVPT(i): denotes the class of the bottom two vertices, the lower boundary state.
!           Classes: for the top vertices
!           0: a top walk segment, or a head segment
!           1: an intermediate segment or a tail segment with delta(b)=-1
!           2: an intermediate segment or a tail segment with delta(b)=+1
!           3: a tail walk segment
!           Classes: for the bottom vertices
!           0: a top walk segment
!           1: an intermediate segment or a head segment with delta(b)=-1
!           2: an intermediate segment or a head segment with delta(b)=+1
!           3: a tail segment or a tail walk segment
!
!           When segments are matched together there tail class, or an upper segments, must match the head class of the
!           lower segment. Matching upper and lower boundaries must be in the same state.
integer(kind=iwp), parameter :: nSegBase = 26, nSegWeight = 8, nSegTot = nSegBase+nSegWeight, nSeg = nSegBase

! Offsets for diagonal segments  (27:30) and (31:34)
integer(kind=iwp), parameter :: iSegWUpBeg = nSegBase+1
integer(kind=iwp), parameter :: iSegWUpEnd = nSegBase+4
integer(kind=iwp), parameter :: iSegWLoBeg = nSegBase+5
integer(kind=iwp), parameter :: iSegWLoEnd = nSegBase+8

! The segments of the triplet operator T_pq = a+_pa a_qa - a+_pb a_qb (ASTA,
! paper II), appended to the catalog (35:62). Three further classes of the
! vertex pairs above the loop head, where the two walks of T_pq may differ:
!           4: db = 0 (b /= 0), 5: db = -2, 6: db = +2 (db = b(right) -
!           b(left), the right/ket vertex of a class-5/6 pair is IVR(:,3/4))
! 35-38: head segments, db = 0 at the top (0u, 0d, u2, d2)
! 39-42: head segments, db = -2 (u2, 0d) and db = +2 (d2, 0u)
! 43-58: top segments, db = 0 (00, uu, dd, 22, ud, du), db = -2 (00, uu, dd,
!        22, ud), db = +2 (00, uu, dd, 22, du)
! 59-62: weight segments of T_pp (uu, dd; ud with db = -2, du with db = +2)
! Their values are those of ASTA_Head_T, ASTA_Top_T and ASTA_Weight_T (MKSEG;
! ISVC unused). The singlet code uses only the segments 1:nSeg.
integer(kind=iwp), parameter :: nSegTrip = 28, nSegCat = nSegTot+nSegTrip
integer(kind=iwp), parameter :: iSegTHeadBeg = nSegTot+1, iSegTHeadEnd = nSegTot+8, iSegTTopBeg = nSegTot+9, &
                                iSegTTopEnd = nSegTot+24, iSegTWBeg = nSegTot+25, iSegTWEnd = nSegCat

integer(kind=iwp), parameter :: &
  ITVPT(nSegCat) = &
  [0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,1,1,2,2,3,3,3,3,0,0,0,0,3,3,3,3, &
   4,4,4,4,5,5,6,6,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,4,4,5,6], &
  IBVPT(nSegCat) = &
  [0,0,0,0,1,1,2,2,1,1,2,1,1,2,2,1,2,2,3,3,3,3,3,3,3,3,0,0,0,0,3,3,3,3, &
   1,2,2,1,1,1,2,2,4,4,4,4,6,5,5,5,5,5,4,6,6,6,6,4,3,3,3,3], &
  IC1(nSegCat) = &
  [0,1,2,3,0,2,0,1,0,1,1,2,3,0,1,2,2,3,1,3,2,3,0,1,2,3,0,1,2,3,0,1,2,3, &
   0,0,1,2,1,0,2,0,0,1,2,3,1,2,0,1,2,3,1,0,1,2,3,2,1,2,1,2], &
  IC2(nSegCat) = &
  [0,1,2,3,1,3,2,3,0,1,2,2,3,0,1,1,2,3,0,2,0,1,0,1,2,3,0,1,2,3,0,1,2,3, &
   1,2,3,3,3,2,3,1,0,1,2,3,2,1,0,1,2,3,2,0,1,2,3,1,1,2,2,1], &
  ISVC(nSegCat) = &
  [1,1,1,1,1,7,8,4,1,2,9,10,2,1,2,11,12,2,1,5,6,3,1,1,1,1,0,1,1,13,0,1,1,13, &
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

integer(kind=iwp), parameter :: TR_WALK = 1
integer(kind=iwp), parameter :: TR_OPEN = 2
integer(kind=iwp), parameter :: TR_MID = 4
integer(kind=iwp), parameter :: TR_CLOSE = 8
integer(kind=iwp), parameter :: TR_TAIL = 16
integer(kind=iwp), parameter :: TR_WEIGHT = 32

!#define _TEST_CSFDET_

! CSFG_nPack (cases per integer of the packed walks) comes from CSF_Graph

type(SGStruct), target :: SGS(4)
type(CIStruct) :: CIS(4)
type(EXStruct) :: EXS(4)
logical(kind=iwp) :: State_is_used(4) = .false.

! Sanity checks of the SGUGA tables and coupling coefficients (CI
! development); set by the CI interface before the initialization.
logical(kind=iwp) :: SG_Checks = .false.
! Tags of CI vectors (SG_New_Tag): a caller that applies several E_pq to the
! same vector passes the same tag to SG_Epq_Psi, which may then reuse its
! intermediates; positive, a new one for every vector.
! SG_No_Reuse: the vector is used only once (no reuse, no fingerprints).
integer(kind=iwp) :: SG_Tag_Last = 0
integer(kind=iwp), parameter :: SG_No_Reuse = -1

! State and CSF symmetry of the single-CSF matrix elements (MatEl_Epq,
! MatEl_EpqErs, MatEl_H_Direct and their helpers), set by
! SG_Test_EpqErs_Ref and by SG_MatEl_Init.
integer(kind=iwp) :: ME_iState = 1, ME_LSM = 1
! Decode_ME: the last two decoded CSFs (symmetry, CSF, MV, ISYUP, IUP, IDWN)
! of CIS(ME_iState); reset by Reset_Decode_ME when ME_iState is (re)set.
integer(kind=iwp) :: ME_Dec(6,2) = 0, ME_DecNext = 1
! ME_Blk(:,k,ISYC): the nonempty (ISYUP,MV) CSF blocks of CIS(ME_iState)
! for CSF symmetry ISYC, (offset, ISYUP, MV), sorted by offset (ME_nBlk of
! them), for a binary search in Decode_ME.
integer(kind=iwp), allocatable :: ME_Blk(:,:,:), ME_nBlk(:)
! MatEl_H_Direct: occupation-only parts kept for the CSF pairs of one
! configuration pair, the diagonal part ME_HD (occupation ME_HD_Occ) and the
! single-excitation coefficient ME_HW (E_ab, ket occupation ME_HW_Occ);
! invalidated (Reset_MatEl_H) whenever the integrals may change.
logical(kind=iwp) :: ME_HD_Valid = .false., ME_HW_Valid = .false.
integer(kind=iwp) :: ME_HW_ab(2) = 0
integer(kind=iwp), allocatable :: ME_HD_Occ(:), ME_HW_Occ(:)
real(kind=wp) :: ME_HD = Zero, ME_HW = Zero
! Coup_Value, Path_Value: blocks up to this size are scanned, not searched
integer(kind=iwp), parameter :: nCoup_Lin = 16
! MatEl_EpqErs_OO: the partial sums over the upper and the lower half of the
! intermediate walks (OO_Half). They are coupling coefficients only (no
! integrals) and recur for many CSF pairs, so they are kept in a
! direct-mapped cache, valid as long as the graph of ME_iState: the key
! (bra walk, ket walk, packed levels/orientations/bands/symmetry/half/state,
! bra midvertex) in OO_Key, the sum in OO_Val; cleared by Reset_Decode_ME,
! freed by SG_Free. Empty slot: OO_Key(4,:) = 0. Only sums with a coupling
! block of more than OO_MinN transitions are kept (see OO_Half).
integer(kind=iwp), parameter :: OO_MinN = 4, nOO_Cache = 8192, OO_MaxLev = 127, OO_MaxState = 511
integer(kind=iwp), allocatable :: OO_Key(:,:)
real(kind=wp), allocatable :: OO_Val(:)

! The split-graph coupling tables of the triplet operator (TRP_Setup): the
! tuples (key, left/right half walk, their symmetries, right midvertex,
! value), the half-walk lookup and the work variables of TRP_Gen.
integer(kind=iwp) :: TRP_iState = 0, TRP_n = 0, TRP_Cap = 0, TRP_Half = 0, TRP_MVL0 = 0, TRP_Start = 0, TRP_HW_n = 0
integer(kind=iwp), allocatable :: TRP_Key(:), TRP_IL(:), TRP_IR(:), TRP_SYL(:), TRP_SYR(:), TRP_MVR(:), TRP_SL(:), TRP_SR(:), &
                                  TRP_HW_K(:,:), TRP_HW_SY(:), TRP_HW_ID(:)
real(kind=wp), allocatable :: TRP_V(:)

! ASTA (SG_Xpq_Psi): step vectors of all CSFs of state ASTA_iState
! (symmetry blocks at ASTA_Off), and per symmetry the packed step vectors
! sorted (ASTA_Key) with their CSF index (ASTA_Idx); see ASTA_Setup.
integer(kind=iwp) :: ASTA_iState = 0, ASTA_Off(8) = 0
integer(kind=iwp), allocatable :: ASTA_Steps(:,:), ASTA_Key(:), ASTA_Idx(:)

! A sigma routine H|Psi> of a split graph (module SGUGA_Sigma: sg_h_psi,
! sg_h_psi_ras), passed to SG_Test_EpqErs_Ref as the reference of its
! <I|H|J> check (SGUGA_Sigma uses this module, so SGUGA cannot use it).
abstract interface
  subroutine SG_HPsi_Proc(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)
    import :: SGStruct, CIStruct, EXStruct, iwp, wp
    type(SGStruct), intent(in) :: SGS
    type(CIStruct), intent(in) :: CIS
    type(EXStruct), intent(inout) :: EXS
    integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
    real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
    real(kind=wp), intent(out) :: Sigma(nCSFs)
  end subroutine SG_HPsi_Proc
end interface

contains

subroutine SG_Init_Light(iState,nSym,nActEl,iSpin, &
                         nRas,nRasEl,nRsPrt, &
                         xLevel,xL2Act,xNLEV,xNSM,Do_MkSGUGA,xnRasMx)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xNLEV, xNSM(MxLev), xnRasMx(nRsPrt)
  logical(kind=iwp), optional, intent(in) :: Do_MkSGUGA
  integer(kind=iwp) :: iSym

  if (State_is_used(iState)) then
    write(u6,*) 'SG_Init_Light: State is already initiated'
    write(u6,*) 'iState: ',iState
    call abend()
  end if

  SGS(istate)%IFRAS = 0
  if (nRsPrt == 3) SGS(istate)%IFRAS = 1
  do iSym=1,nSym
    if (sum(nRas(iSym,1:nRsPrt)) /= 0) SGS(istate)%IFRAS = SGS(istate)%IFRAS+1
  end do
  SGS(istate)%nRasEL(1:nRsPrt) = nRasEl(1:nRsPrt)
  SGS(istate)%nRas(1:MxSym,1:nRsPrt) = nRas(1:MxSym,1:nRsPrt)
  SGS(istate)%nRsPrt = nRsPrt
  ! GAS: upper bounds on the electrons in the partitions 1..k
  SGS(istate)%nRasMx(:) = -1
  if (present(xnRasMx)) SGS(istate)%nRasMx(1:nRsPrt) = xnRasMx(1:nRsPrt)

  ! Make sure that we start from a clean slate.
  call SG_Free(istate)
  State_is_used(iState) = .true.

  if ((nSym < 1) .or. (nSym > 8)) then
    write(u6,*) ' SG_Init_Light: illegal nSym value:',nSym
    call Abend()
  end if
  if (iSpin < 1) then
    write(u6,*) ' SG_Init_Light: illegal iSpin value:',iSpin
    call Abend()
  end if
  if (nActEl < 0) then
    write(u6,*) ' SG_Init_Light: illegal nActEl value:',nActEl
    call Abend()
  end if

  SGS(istate)%nSym = nSym
  SGS(istate)%iSpin = iSpin
  SGS(istate)%nActEl = nActEl

  if (present(xLevel)) SGS(istate)%Level(:) = xLevel(:)
  if (present(xL2Act)) SGS(istate)%L2Act(:) = xL2Act(:)
! Initiate if not already set externally.
  if (SGS(istate)%LEVEL(1) == 0) SGS(istate)%LEVEL(1:SGS(istate)%nLev) = [(iq,iq=1,SGS(istate)%nLev)]
  if (SGS(istate)%L2Act(1) == 0) SGS(istate)%L2Act(1:SGS(istate)%nLev) = [(iq,iq=1,SGS(istate)%nLev)]

! CREATE THE SYMMETRY INDEX VECTOR

  SGS(istate)%NLEV = xnLEV
! Allocate Level to Symmetry table ISm:
  call mma_allocate(SGS(istate)%ISM,SGS(istate)%nLev,Label='SGS%ISM')
  SGS(istate)%ISM(1:SGS(istate)%nLev) = xNSM(1:SGS(istate)%nLev)

  if (present(Do_MkSGUGA)) then
    if (Do_MkSGUGA) then
      call CSFG_Make(SGS(istate),CIS(istate))

    else
      SGS(istate)%iSpin = 0
      SGS(istate)%nActEl = 0

      ! INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
      call mma_allocate(CIS(istate)%NCSF,SGS(istate)%nSym,Label='CIS%NCSF')
      CIS(istate)%NCSF(:) = 0
    end if
  else
    call CSFG_Make(SGS(istate),CIS(istate))
  end if

  if (NActEl == 0) CIS(istate)%NCSF(:) = 1

end subroutine SG_Init_Light

subroutine SG_Free(iState)
! PURPOSE: FREE THE SGUGA TABLES

  integer(kind=iwp), intent(in) :: iState

  integer(kind=iwp) :: iOpen, mult, iRefSm

  if (.not. State_is_used(iState)) return
  if (ASTA_iState == iState) call ASTA_Free()
  if (TRP_iState == iState) call TRP_Free()

  associate(SGS => SGS(iState),CIS => CIS(iState),EXS => EXS(iState))

    call mma_deallocate(SGS%ISM,safe='*')
    call mma_deallocate(SGS%DRT0,safe='*')
    call mma_deallocate(SGS%DOWN0,safe='*')
    call mma_deallocate(SGS%DRT,safe='*')
    call mma_deallocate(SGS%DOWN,safe='*')
    call mma_deallocate(SGS%UP,safe='*')
    call mma_deallocate(SGS%MAW,safe='*')
    call mma_deallocate(SGS%LTV,safe='*')
    call mma_deallocate(SGS%DAW,safe='*')
    call mma_deallocate(SGS%RAW,safe='*')
    call mma_deallocate(SGS%Ver,safe='*')

    call mma_deallocate(SGS%OC,safe='*')
    if (allocated(SGS%CDCInfo)) then

      do iOpen=lbound(SGS%CDCInfo,1),ubound(SGS%CDCInfo,1)

        do mult=lbound(SGS%CDCInfo,2),ubound(SGS%CDCInfo,2)

          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOpenOCP,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOCP,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%DOC,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%SOC,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%OCReverse,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%SOCSym,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOCSym,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%OCDOC,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%OCSOC,safe='*')

          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%detspin,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%twos,safe='*')

          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%cdc,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%dtoc,safe='*')
          call mma_deallocate(SGS%CDCInfo(iOpen,mult)%pair,safe='*')

        end do

      end do

      deallocate(SGS%CDCInfo)

    end if

    if (allocated(SGS%SDREO)) then
      do iRefSm=1,size(SGS%SDREO)
        call mma_deallocate(SGS%SDREO(iRefSm)%List,Safe='*')
      end do
      deallocate(SGS%SDREO)
    end if

    nullify(SGS%DRTP,SGS%DOWNP)

    call mma_deallocate(CIS%NOW,safe='*')
    call mma_deallocate(CIS%IOW,safe='*')
    call mma_deallocate(CIS%NCSF,safe='*')
    call mma_deallocate(CIS%NOCSF,safe='*')
    call mma_deallocate(CIS%IOCSF,safe='*')
    call mma_deallocate(ME_Blk,safe='*')
    call mma_deallocate(ME_nBlk,safe='*')
    call mma_deallocate(ME_HD_Occ,safe='*')
    call mma_deallocate(ME_HW_Occ,safe='*')
    call mma_deallocate(OO_Key,safe='*')
    call mma_deallocate(OO_Val,safe='*')
    call Reset_MatEl_H()
    call mma_deallocate(CIS%ICase,safe='*')
    call mma_deallocate(CIS%VSGM,safe='*')
    call mma_deallocate(CIS%IVR,safe='*')
    call mma_deallocate(CIS%ISGM,safe='*')

    call mma_deallocate(EXS%NOCP,safe='*')
    call mma_deallocate(EXS%IOCP,safe='*')
    call mma_deallocate(EXS%ICoup,safe='*')
    call mma_deallocate(EXS%IOrd,safe='*')
    call mma_deallocate(EXS%VTab,safe='*')
    call mma_deallocate(EXS%SGTMP,safe='*')
    EXS%Save_Valid = .false.
    call mma_deallocate(EXS%MVL,safe='*')
    call mma_deallocate(EXS%MVR,safe='*')
    call mma_deallocate(EXS%USGN,safe='*')
    call mma_deallocate(EXS%LSGN,safe='*')
    EXS%Num_Sym = 0
    call mma_deallocate(EXS%I1list,safe='*')
    call mma_deallocate(EXS%I2list,safe='*')
    call mma_deallocate(EXS%Xlist,safe='*')

  end associate

  State_is_used(iState) = .false.

end subroutine SG_Free

!> A new tag for a CI vector (see SG_Epq_Psi, CI_Tag).
function SG_New_Tag() result(iTag)
  integer(kind=iwp) :: iTag
  SG_Tag_Last = SG_Tag_Last+1
  if (SG_Tag_Last >= huge(SG_Tag_Last)) SG_Tag_Last = 1
  iTag = SG_Tag_Last
end function SG_New_Tag

!> The binomial coefficient (n over k), the number of k-subsets of n
!> elements (as binom_coef of Faroald's second_quantization; SGUGA keeps its
!> own copy).
function binomial(k,n) result(b)
  integer(kind=iwp), intent(in) :: k, n
  integer(kind=iwp) :: b
  integer(kind=iwp) :: i
  ! exactly as binom_coef: 0 for k > n, 1 for k <= 0
  b = 0
  if (k > n) return
  b = 1
  if (k <= 0) return
  do i=1,min(k,n-k)
    b = (b*(n-min(k,n-k)+i))/i
  end do
end function binomial

function SG_NUM(SGS,EXS,IWALK)
! PURPOSE: FOR ANY GIVEN WALK (STEP VECTOR) COMPUTE THE
!          LEXICAL NUMBER IN THE SPLIT GUGA REPRESENTATION

  use Definitions, only: iwp

  implicit none
  integer(kind=iwp) :: SG_NUM
  type(SGStruct), intent(in) :: SGS
  type(EXStruct), intent(in) :: EXS
  integer(kind=iwp), intent(in) :: IWALK(SGS%nLev)
  integer(kind=iwp) :: IC, ICASE, ICONF, IDAWSUM, IRAWSUM, IUW, LEV, LV, MIDV

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

  MIDV = 1
  do LEV=SGS%nLev,(SGS%MidLev+1),-1
    ICASE = IWALK(LEV)
    MIDV = SGS%DOWN(MIDV,ICASE)
  end do
  MIDV = MIDV-SGS%MVSta+1

! FIND REVERSE ARC WEIGHT FOR THE UPPER WALK

  IRAWSUM = 1
  LV = 1
  do LEV=SGS%nLev,(SGS%MidLev+1),-1
    IC = IWALK(LEV)
    LV = SGS%DOWN(LV,IC)
    IRAWSUM = IRAWSUM+SGS%RAW(LV,IC)
  end do
  IUW = EXS%USGN(IRAWSUM,MIDV)

! FIND DIRECT ARC WEIGHT FOR THE LOWER WALK

  IDAWSUM = 1
  LV = SGS%nVert
  do LEV=1,SGS%MidLev
    IC = IWALK(LEV)
    LV = SGS%UP(LV,IC)
    IDAWSUM = IDAWSUM+SGS%DAW(LV,IC)
  end do
  ICONF = EXS%LSGN(IDAWSUM,MIDV)

! COMPUTE LEXICAL ORDERING NUMBER

  SG_NUM = IUW+ICONF

end function SG_NUM

function SG_PHASE(SGS,IWALK)
! PURPOSE: THE SYMMETRIC GROUP APPROACH AND THE UNITARY GROUP
!          APPROACH DIFFER IN THE PHASE CONVENTION. FIND THE
!          PHASE FACTOR RELATING THE CSFS IN EITHER BASIS.

  use Definitions, only: iwp

  implicit none
  integer(kind=iwp) :: SG_PHASE
  type(SGStruct), intent(in) :: SGS
  integer(kind=iwp), intent(in) :: IWALK(SGS%NLEV)
  integer(kind=iwp) :: ICASE, ISGN, IVERT, LEV

! FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY

  SG_PHASE = 1
  IVERT = SGS%NVERT
  do LEV=1,SGS%NLEV
    ICASE = IWALK(LEV)
    IVERT = SGS%UP(IVERT,ICASE)
    select case (iCase)
      case (2,3)
        ISGN = (-1)**SGS%DRT(IVERT,4)
      case default
        ISGN = 1
    end select
    SG_PHASE = SG_PHASE*ISGN
  end do

end function SG_PHASE

!> Decode_CSF for CIS(ME_iState), remembering the last two CSFs (least
!> recently used slot replaced): the matrix-element routines decode the same
!> bra and ket for every operator string of a matrix element. A new CSF is
!> located by a binary search over the sorted block offsets ME_Blk.
subroutine Decode_ME(ISYCI,iCSF,MV,ISYUP,IUP,IDWN)

  use Definitions, only: iwp

  integer(kind=iwp), intent(in) :: ISYCI, iCSF
  integer(kind=iwp), intent(out) :: MV, ISYUP, IUP, IDWN
  integer(kind=iwp) :: k, iLo, iHi, iMid, ILOC, NUP

  do k=1,2
    if ((ME_Dec(1,k) == ISYCI) .and. (ME_Dec(2,k) == iCSF)) then
      MV = ME_Dec(3,k)
      ISYUP = ME_Dec(4,k)
      IUP = ME_Dec(5,k)
      IDWN = ME_Dec(6,k)
      ME_DecNext = 3-k
      return
    end if
  end do
  if (.not. allocated(ME_Blk)) call Reset_Decode_ME()
! last block with offset < iCSF
  iLo = 1
  iHi = ME_nBlk(ISYCI)
  if (iHi < 1) call Abend()
  do while (iLo < iHi)
    iMid = (iLo+iHi+1)/2
    if (ME_Blk(1,iMid,ISYCI) < iCSF) then
      iLo = iMid
    else
      iHi = iMid-1
    end if
  end do
  ISYUP = ME_Blk(2,iLo,ISYCI)
  MV = ME_Blk(3,iLo,ISYCI)
  ILOC = iCSF-ME_Blk(1,iLo,ISYCI)
  if ((ILOC < 1) .or. (ILOC > CIS(ME_iState)%NOCSF(ISYUP,MV,ISYCI))) call Abend()
  NUP = CIS(ME_iState)%NOW(1,ISYUP,MV)
  IDWN = (ILOC-1)/NUP+1
  IUP = ILOC-(IDWN-1)*NUP
  k = ME_DecNext
  ME_Dec(:,k) = [ISYCI,iCSF,MV,ISYUP,IUP,IDWN]
  ME_DecNext = 3-k

end subroutine Decode_ME

!> Forget the decoded CSFs and rebuild ME_Blk for CIS(ME_iState).
subroutine Reset_Decode_ME()

  integer(kind=iwp) :: nSy, nMV, ISYC, ISYU, MV, n, k, Tmp(3)

  ME_Dec(:,:) = 0
  ME_DecNext = 1
  call Reset_MatEl_H()
  if (allocated(OO_Key)) OO_Key(4,:) = 0
  call mma_deallocate(ME_Blk,safe='*')
  call mma_deallocate(ME_nBlk,safe='*')
  if (.not. allocated(CIS(ME_iState)%NOCSF)) return
  nSy = size(CIS(ME_iState)%NOCSF,1)
  nMV = CIS(ME_iState)%nMidV
  call mma_allocate(ME_Blk,3,max(1,nSy*nMV),size(CIS(ME_iState)%NOCSF,3),Label='ME_Blk')
  call mma_allocate(ME_nBlk,size(CIS(ME_iState)%NOCSF,3),Label='ME_nBlk')
  do ISYC=1,size(CIS(ME_iState)%NOCSF,3)
    n = 0
    do MV=1,nMV
      do ISYU=1,nSy
        if (CIS(ME_iState)%NOCSF(ISYU,MV,ISYC) <= 0) cycle
        n = n+1
        ME_Blk(:,n,ISYC) = [CIS(ME_iState)%IOCSF(ISYU,MV,ISYC),ISYU,MV]
        ! insertion sort by offset
        k = n
        do while (k > 1)
          if (ME_Blk(1,k-1,ISYC) <= ME_Blk(1,k,ISYC)) exit
          Tmp(:) = ME_Blk(:,k,ISYC)
          ME_Blk(:,k,ISYC) = ME_Blk(:,k-1,ISYC)
          ME_Blk(:,k-1,ISYC) = Tmp(:)
          k = k-1
        end do
      end do
    end do
    ME_nBlk(ISYC) = n
  end do

end subroutine Reset_Decode_ME

!> Invalidate the occupation-only parts kept by MatEl_H_Direct.
subroutine Reset_MatEl_H()
  ME_HD_Valid = .false.
  ME_HW_Valid = .false.
end subroutine Reset_MatEl_H

subroutine Keep_Occ(Kept,Occ,nL)
  integer(kind=iwp), allocatable, intent(inout) :: Kept(:)
  integer(kind=iwp), intent(in) :: Occ(:), nL
  if (allocated(Kept)) then
    if (size(Kept) < nL) call mma_deallocate(Kept)
  end if
  if (.not. allocated(Kept)) call mma_allocate(Kept,nL,Label='ME_H_Occ')
  Kept(1:nL) = Occ(1:nL)
end subroutine Keep_Occ

subroutine Decode_CSF( &
  CIS,ISYCI, &
  iCSF, &
  MV,ISYUP,IUP,IDWN)

  use Definitions, only: iwp, u6

  type(CIStruct), intent(in) :: CIS

  integer(kind=iwp), intent(in) :: ISYCI
  integer(kind=iwp), intent(in) :: iCSF

  integer(kind=iwp), intent(out) :: MV
  integer(kind=iwp), intent(out) :: ISYUP
  integer(kind=iwp), intent(out) :: IUP
  integer(kind=iwp), intent(out) :: IDWN

  integer(kind=iwp) :: IOFF
  integer(kind=iwp) :: NCI
  integer(kind=iwp) :: NUP
  integer(kind=iwp) :: ILOC

  do MV=1,CIS%nMidV

    do ISYUP=1,size(CIS%NOCSF,1)

      IOFF = CIS%IOCSF(ISYUP,MV,ISYCI)
      NCI = CIS%NOCSF(ISYUP,MV,ISYCI)

      ! CSF iCSF lies in this (MV,ISYUP) block
      if ((NCI > 0) .and. (iCSF > IOFF) .and. (iCSF <= IOFF+NCI)) then

        ILOC = iCSF-IOFF

        NUP = CIS%NOW(1,ISYUP,MV)

        IDWN = (ILOC-1)/NUP+1
        IUP = ILOC-(IDWN-1)*NUP

        return

      end if

    end do

  end do

  write(u6,*) 'Decode_CSF failed'
  write(u6,*) 'iCSF=',iCSF
  call Abend()

end subroutine Decode_CSF

subroutine MK_FLAT_Epq(SGS)
  integer(kind=iwp) p, q, iSym
  type(SGStruct), intent(inout) :: SGS

! The purpose of this routine is to remove all the logic from sg_Epq_Psi.
! We note that there are two modes in which the routine is used:
! 1) it generates all entities, and
! 2) it reuse some of the sigma vectors.
! This results in two different sequences of calls, due to different logic.
! We have to set up flat tables corresponding to both of these cases.
! ToDo: Develop this further.
  return
  do iSym=1,SGS%nSym
    do p=1,SGS%nLev
      do q=1,SGS%nLev

      end do
    end do
  end do

end subroutine MK_FLAT_Epq

! PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
!          AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
!          RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
!          VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
!          CORRESPONDING STEP VECTOR.
subroutine STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)

  use Definitions, only: iwp

  implicit none

  integer(iwp), intent(in) :: NCL, NOP, NORB
  integer(iwp), intent(in) :: ICL(NCL), IOP(NOP)
  integer(iwp), intent(in) :: ISPIN(NOP)
  integer(iwp), intent(out) :: IWALK(NORB)

  integer(iwp) :: iClosed, iOpen
!
!  Start from an empty occupation pattern
!
  IWALK(:) = 0
!
!  Mark doubly occupied orbitals
!
  do iClosed=1,NCL
    IWALK(ICL(iClosed)) = 3
  end do
!
!  Mark singly occupied orbitals
!
  do iOpen=1,NOP
    if (ISPIN(iOpen) == 1) then
      IWALK(IOP(iOpen)) = 1
    else
      IWALK(IOP(iOpen)) = 2
    end if
  end do

end subroutine STEPVEC

! iTri of Index_functions, here in the module so that the compiler can
! inline it in the integral lookups of MatEl_H_Direct.
pure function Tri_SG(i,j) result(ij)

  integer(kind=iwp), intent(in) :: i, j
  integer(kind=iwp) :: ij

  ij = merge(j*(j-1)/2+i,i*(i-1)/2+j,j > i)

end function Tri_SG

!
! Orbital occupations of CSF I of symmetry ISY, unpacked from its step
! vector: the upper and lower walks are stored in CIS%ICase, 2 bits per
! level (step d = 0,1,2,3 gives occupation 0,1,1,2), CSFG_nPack levels per
! integer. Same layout as read in the diagonal case of MatEl_Epq.
!
subroutine CSF_Occupation(ISY,I,OccVec)
  integer(kind=iwp), intent(in) :: ISY, I
  integer(kind=iwp), intent(out) :: OccVec(:)
  integer(kind=iwp) :: MV, ISYU, ISYD, IUP, IDN, iLev, IOW0, IPSHFT, IC, ICS
  call Decode_ME(ISY,I,MV,ISYU,IUP,IDN)
  ISYD = Mul(ISY,ISYU)
  do iLev=1,SGS(ME_iState)%nLev
    if (iLev > SGS(ME_iState)%MidLev) then
      IOW0 = CIS(ME_iState)%IOW(1,ISYU,MV)
      IPSHFT = 2*(iLev-1-SGS(ME_iState)%MidLev)
      IC = CIS(ME_iState)%ICase(1+IOW0-CIS(ME_iState)%nIpWlk+IPSHFT/(2*CSFG_nPack)+IUP*CIS(ME_iState)%nIpWlk)
    else
      IOW0 = CIS(ME_iState)%IOW(2,ISYD,MV)
      IPSHFT = 2*(iLev-1)
      IC = CIS(ME_iState)%ICase(1+IOW0-CIS(ME_iState)%nIpWlk+IPSHFT/(2*CSFG_nPack)+IDN*CIS(ME_iState)%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*CSFG_nPack)
    ICS = mod(IC/2**IPSHFT,4)
    OccVec(iLev) = (1+ICS)/2
  end do
end subroutine CSF_Occupation

!***********************************************************************
! Matrix elements of the singlet E_pq and of the triplet operator
! T_pq = a+_pa a_qa - a+_pb a_qb (sqrt(2) T^{1,0}_pq, M = 0) between
! two CSFs of the same graph, as products of ASTA segment values (the
! algebraic spin-tensor approach: Fdez. Galvan, Rooein, Lindh, papers I
! and II). Steps 0,1,2,3 = 0,u,d,2, level 1 at the graph tail; b = 2S at
! the top of a level. As in the papers, U is the bra and T the ket of a
! raising operator (p < q: tail at p, head at q); every segment value
! refers to the b of the top vertex of the ket walk, and db = bT - bU is
! taken at the top of the segment. A lowering operator is the transpose,
! <I|X_pq|J> = <J|X_qp|I>. Above the loop head the walks of E_pq coincide;
! those of T_pq may differ by db = 0 (b /= 0), -2 or +2 (top segments).
!***********************************************************************

!> Step vector of CSF I of symmetry ISY (0,1,2,3 = 0,u,d,2 per level,
!> level 1 at the graph tail), from the packed upper and lower walks.
subroutine SG_CSF_StepVec(SGS,CIS,ISY,I,Step)
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: ISY, I
  integer(kind=iwp), intent(out) :: Step(SGS%nLev)
  integer(kind=iwp) :: MV, ISYU, ISYD, IUP, IDN, iLev, IOW0, IPSHFT, IC
  call Decode_CSF(CIS,ISY,I,MV,ISYU,IUP,IDN)
  ISYD = Mul(ISY,ISYU)
  do iLev=1,SGS%nLev
    if (iLev > SGS%MidLev) then
      IOW0 = CIS%IOW(1,ISYU,MV)
      IPSHFT = 2*(iLev-1-SGS%MidLev)
      IC = CIS%ICase(1+IOW0-CIS%nIpWlk+IPSHFT/(2*CSFG_nPack)+IUP*CIS%nIpWlk)
    else
      IOW0 = CIS%IOW(2,ISYD,MV)
      IPSHFT = 2*(iLev-1)
      IC = CIS%ICase(1+IOW0-CIS%nIpWlk+IPSHFT/(2*CSFG_nPack)+IDN*CIS%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*CSFG_nPack)
    Step(iLev) = mod(IC/2**IPSHFT,4)
  end do
end subroutine SG_CSF_StepVec

!> Map(I) = the index in state iB of the CSF I (symmetry ISY) of state iA,
!> identified by its step vector, or 0 if state iB does not hold it. The
!> two states must have the same levels (e.g. a GAS space and its relaxed
!> space, which contains it).
subroutine SG_Map_CSFs(iA,iB,ISY,nA,Map)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iA, iB, ISY, nA
  integer(kind=iwp), intent(out) :: Map(nA)
  integer(kind=iwp), parameter :: nPerWord = (bit_size(0_iwp)-2)/2
  integer(kind=iwp) :: nL, nW, nB, I, lo, hi, md, c
  integer(kind=iwp), allocatable :: Step(:), KB(:,:), KA(:), Ord(:)

  nL = SGS(iA)%nLev
  if (SGS(iB)%nLev /= nL) then
    write(u6,*) 'SG_Map_CSFs: the states have different levels',nL,SGS(iB)%nLev
    call Abend()
  end if
  if (.not. allocated(CIS(iA)%ICase)) call CSFG_MkCOT(SGS(iA),CIS(iA))
  if (.not. allocated(CIS(iB)%ICase)) call CSFG_MkCOT(SGS(iB),CIS(iB))
  nW = max(1,(nL+nPerWord-1)/nPerWord)
  nB = CIS(iB)%NCSF(ISY)
  call mma_allocate(Step,max(1,nL),Label='Map_Step')
  call mma_allocate(KB,nW,max(1,nB),Label='Map_KB')
  call mma_allocate(KA,nW,Label='Map_KA')
  call mma_allocate(Ord,max(1,nB),Label='Map_Ord')
  do I=1,nB
    call SG_CSF_StepVec(SGS(iB),CIS(iB),ISY,I,Step)
    call Pack_Step(Step,KB(:,I))
  end do
  call Sort_Keys(KB,nW,nB,Ord)
  do I=1,nA
    call SG_CSF_StepVec(SGS(iA),CIS(iA),ISY,I,Step)
    call Pack_Step(Step,KA)
    Map(I) = 0
    lo = 1
    hi = nB
    do while (lo <= hi)
      md = (lo+hi)/2
      c = Cmp_Key(KA,KB(:,Ord(md)))
      if (c == 0) then
        Map(I) = Ord(md)
        exit
      else if (c < 0) then
        hi = md-1
      else
        lo = md+1
      end if
    end do
  end do
  call mma_deallocate(Ord)
  call mma_deallocate(KA)
  call mma_deallocate(KB)
  call mma_deallocate(Step)

contains

  subroutine Pack_Step(St,K)
    integer(kind=iwp), intent(in) :: St(:)
    integer(kind=iwp), intent(out) :: K(:)
    integer(kind=iwp) :: iLev, iW
    K(:) = 0
    do iLev=1,nL
      iW = (iLev-1)/nPerWord+1
      K(iW) = K(iW)+St(iLev)*4_iwp**(iLev-1-(iW-1)*nPerWord)
    end do
  end subroutine Pack_Step

  function Cmp_Key(X,Y) result(c)
    integer(kind=iwp), intent(in) :: X(:), Y(:)
    integer(kind=iwp) :: c, w
    c = 0
    do w=size(X),1,-1
      if (X(w) /= Y(w)) then
        c = merge(-1,1,X(w) < Y(w))
        return
      end if
    end do
  end function Cmp_Key

end subroutine SG_Map_CSFs

!> b = 2S at the top of every level of a walk (B(0) = 0, the graph tail).
pure subroutine Walk_B(nLev,Step,B)
  integer(kind=iwp), intent(in) :: nLev, Step(nLev)
  integer(kind=iwp), intent(out) :: B(0:nLev)
  integer(kind=iwp) :: k
  B(0) = 0
  do k=1,nLev
    select case (Step(k))
      case (1)
        B(k) = B(k-1)+1
      case (2)
        B(k) = B(k-1)-1
      case default
        B(k) = B(k-1)
    end select
  end do
end subroutine Walk_B

pure function Sgn(n) result(s)
  integer(kind=iwp), intent(in) :: n
  real(kind=wp) :: s
  s = One
  if (mod(abs(n),2) == 1) s = -One
end function Sgn

pure function Rt(n,m) result(r)
  ! sqrt(n/m)
  integer(kind=iwp), intent(in) :: n, m
  real(kind=wp) :: r
  r = sqrt(real(n,kind=wp)/real(m,kind=wp))
end function Rt

!> Tail segment of a raising loop (paper I, Table 3, R_), bra u, ket t.
pure function ASTA_Tail(u,t,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, b
  real(kind=wp) :: V
  V = Zero
  select case (4*u+t)
    case (4)   ! u0
      V = One
    case (8)   ! d0
      V = Sgn(b-1)*Rt(b,b+1)
    case (13)  ! 2u
      V = Sgn(b-1)
    case (14)  ! 2d
      V = Rt(b+1,b+2)
  end select
end function ASTA_Tail

!> Intermediate segment of a raising loop (paper I, Table 4, R), db = +-1.
pure function ASTA_Mid(u,t,db,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, db, b
  real(kind=wp) :: V
  V = Zero
  if (abs(db) /= 1) return
  select case (4*u+t)
    case (0)   ! 00
      V = One
    case (5)   ! uu
      V = -One
    case (6)   ! ud
      if (db == -1) V = Sgn(b+1)/sqrt(real((b+2)*(b+1),kind=wp))
    case (9)   ! du
      if (db == 1) V = Sgn(b+1)/sqrt(real((b+1)*b,kind=wp))
    case (10)  ! dd
      if (db == -1) then
        V = sqrt(real((b+3)*(b+1),kind=wp))/real(b+2,kind=wp)
      else
        V = sqrt(real((b+2)*b,kind=wp))/real(b+1,kind=wp)
      end if
    case (15)  ! 22
      V = -One
  end select
end function ASTA_Mid

!> Head segment of a raising loop of E_pq (paper I, Table 2, R^-), db = 0.
pure function ASTA_Head(u,t,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, b
  real(kind=wp) :: V
  V = Zero
  select case (4*u+t)
    case (1)   ! 0u
      V = One
    case (2)   ! 0d
      V = Sgn(b)*Rt(b+2,b+1)
    case (7)   ! u2
      V = Sgn(b+1)
    case (11)  ! d2
      V = Rt(b+2,b+1)
  end select
end function ASTA_Head

!> Head segment of a raising loop of T_pq (paper II, Table 4), db = 0, +-2.
pure function ASTA_Head_T(u,t,db,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, db, b
  real(kind=wp) :: V
  V = Zero
  if ((db == 0) .and. (b == 0)) return
  select case (4*u+t)
    case (7)   ! u2
      if (db == -2) V = Sgn(b)*Two/sqrt(real(b+2,kind=wp))
      if (db == 0) V = Sgn(b)
    case (2)   ! 0d
      if (db == -2) V = Sgn(b+1)*Two*sqrt(real(b+1,kind=wp))/real(b+2,kind=wp)
      if (db == 0) V = Sgn(b+1)*real(b,kind=wp)/sqrt(real((b+1)*(b+2),kind=wp))
    case (11)  ! d2
      if (db == 0) V = real(b,kind=wp)/sqrt(real((b+1)*(b+2),kind=wp))
      if (db == 2) V = -Two*sqrt(real(b-1,kind=wp))/real(b,kind=wp)
    case (1)   ! 0u
      if (db == 0) V = One
      if (db == 2) V = -Two/sqrt(real(b,kind=wp))
  end select
end function ASTA_Head_T

!> Top segment of T_pq above the loop head (paper II, Table 5), db = 0
!> (b /= 0), -2, +2; the same for raising and lowering operators.
pure function ASTA_Top_T(u,t,db,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, db, b
  real(kind=wp) :: V
  V = Zero
  if ((db /= 0) .and. (abs(db) /= 2)) return
  if ((db == 0) .and. (b == 0)) return
  select case (4*u+t)
    case (0,15)  ! 00, 22
      V = One
    case (5)     ! uu
      if (db == -2) V = Rt(b+1,b+2)
      if (db == 0) V = One
      if (db == 2) V = Rt(b-1,b)
    case (10)    ! dd
      if (db == -2) V = sqrt(real((b+1)*(b+4),kind=wp))/real(b+2,kind=wp)
      if (db == 0) V = real(b*(b+3),kind=wp)/real((b+2)*(b+1),kind=wp)
      if (db == 2) V = sqrt(real((b-1)*(b+2),kind=wp))/real(b,kind=wp)
    case (6)     ! ud
      if (db == -2) V = Sgn(b)*Two/(real(b+2,kind=wp)*sqrt(real(b+1,kind=wp)))
      if (db == 0) V = Sgn(b+1)/sqrt(real(b+2,kind=wp))
    case (9)     ! du
      if (db == 0) V = Sgn(b+1)/sqrt(real(b+2,kind=wp))
      if (db == 2) V = Sgn(b)*Two/(real(b,kind=wp)*sqrt(real(b-1,kind=wp)))
  end select
end function ASTA_Top_T

!> Weight segment of T_pp (paper II, Table 2): the diagonal steps and the
!> off-diagonal ud (db = -2) and du (db = +2).
pure function ASTA_Weight_T(u,t,db,b) result(V)
  integer(kind=iwp), intent(in) :: u, t, db, b
  real(kind=wp) :: V
  V = Zero
  select case (4*u+t)
    case (5)   ! uu
      V = One
    case (10)  ! dd
      V = -real(b,kind=wp)/real(b+2,kind=wp)
    case (6)   ! ud
      if (db == -2) V = Sgn(b+1)*Two*sqrt(real(b+1,kind=wp))/real(b+2,kind=wp)
    case (9)   ! du
      if (db == 2) V = Sgn(b-1)*Two*sqrt(real(b-1,kind=wp))/real(b,kind=wp)
  end select
end function ASTA_Weight_T

subroutine ASTA_Free()
  call mma_deallocate(ASTA_Steps,safe='*')
  call mma_deallocate(ASTA_Key,safe='*')
  call mma_deallocate(ASTA_Idx,safe='*')
  ASTA_iState = 0
end subroutine ASTA_Free

subroutine TRP_Free()
  call mma_deallocate(TRP_Key,safe='*')
  call mma_deallocate(TRP_IL,safe='*')
  call mma_deallocate(TRP_IR,safe='*')
  call mma_deallocate(TRP_SYL,safe='*')
  call mma_deallocate(TRP_SYR,safe='*')
  call mma_deallocate(TRP_MVR,safe='*')
  call mma_deallocate(TRP_V,safe='*')
  call mma_deallocate(TRP_HW_K,safe='*')
  call mma_deallocate(TRP_HW_SY,safe='*')
  call mma_deallocate(TRP_HW_ID,safe='*')
  TRP_n = 0
  TRP_Cap = 0
  TRP_iState = 0
end subroutine TRP_Free

!> Index of the coupling tuples by walk, EXS%IOrd (see EXStruct): each NOCP
!> block sorted by ICoup(1,:) and by ICoup(2,:) (heapsort, O(N log N)).
subroutine Mk_Coup_Index(EXS)
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp) :: INDEO, iSym, MV, N

  call mma_allocate(EXS%IOrd,2,max(1,size(EXS%ICoup,2)),Label='EXS%IOrd')
  do MV=1,size(EXS%NOCP,3)
    do iSym=1,size(EXS%NOCP,2)
      do INDEO=1,size(EXS%NOCP,1)
        N = EXS%NOCP(INDEO,iSym,MV)
        if (N > 0) call Index_Coup_Block(EXS,EXS%IOCP(INDEO,iSym,MV),N)
      end do
    end do
  end do

end subroutine Mk_Coup_Index

!> Walk index EXS%IOrd(:,Off+1:Off+N) of one coupling block. Also called by
!> SG_Epq_Psi, which reorders the tuples of a block in place.
subroutine Index_Coup_Block(EXS,Off,N)
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: Off, N
  integer(kind=iwp) :: k, i
  integer(kind=iwp), allocatable :: Key(:), Idx(:)

  call mma_allocate(Key,N,Label='Key')
  call mma_allocate(Idx,N,Label='Idx')
  do k=1,2
    do i=1,N
      Idx(i) = i
      Key(i) = EXS%ICoup(k,Off+i)
    end do
    call Heap_Sort_Index(N,Key,Idx)
    EXS%IOrd(k,Off+1:Off+N) = Idx(:)
  end do
  call mma_deallocate(Key)
  call mma_deallocate(Idx)

end subroutine Index_Coup_Block

!> Sort Key(1:N) in increasing order (heapsort), permuting Idx(1:N) along.
subroutine Heap_Sort_Index(N,Key,Idx)
  integer(kind=iwp), intent(in) :: N
  integer(kind=iwp), intent(inout) :: Key(N), Idx(N)
  integer(kind=iwp) :: i, iEnd

  do i=N/2,1,-1
    call Sift_Down(i,N)
  end do
  do iEnd=N,2,-1
    call Swap_Entries(1,iEnd)
    call Sift_Down(1,iEnd-1)
  end do

contains

  subroutine Sift_Down(iStart,iLast)
    integer(kind=iwp), intent(in) :: iStart, iLast
    integer(kind=iwp) :: iRoot, iChild
    iRoot = iStart
    do
      iChild = 2*iRoot
      if (iChild > iLast) exit
      if (iChild < iLast) then
        if (Key(iChild+1) > Key(iChild)) iChild = iChild+1
      end if
      if (Key(iRoot) >= Key(iChild)) exit
      call Swap_Entries(iRoot,iChild)
      iRoot = iChild
    end do
  end subroutine Sift_Down

  subroutine Swap_Entries(i1,i2)
    integer(kind=iwp), intent(in) :: i1, i2
    integer(kind=iwp) :: iTmp
    iTmp = Key(i1)
    Key(i1) = Key(i2)
    Key(i2) = iTmp
    iTmp = Idx(i1)
    Idx(i1) = Idx(i2)
    Idx(i2) = iTmp
  end subroutine Swap_Entries

end subroutine Heap_Sort_Index

!> Positions iFirst..iLast (within the NOCP block of N tuples at offset Off,
!> order X%IOrd(k,:)) of the tuples with ICoup(k,:) = iWalk (binary search);
!> iFirst > iLast if there are none.
subroutine Coup_Range(X,Off,N,k,iWalk,iFirst,iLast)
  type(EXStruct), intent(in) :: X
  integer(kind=iwp), intent(in) :: Off, N, k, iWalk
  integer(kind=iwp), intent(out) :: iFirst, iLast
  integer(kind=iwp) :: iLo, iHi, iMid

  ! first position with key >= iWalk
  iLo = 1
  iHi = N+1
  do while (iLo < iHi)
    iMid = (iLo+iHi)/2
    if (X%ICoup(k,Off+X%IOrd(k,Off+iMid)) < iWalk) then
      iLo = iMid+1
    else
      iHi = iMid
    end if
  end do
  iFirst = iLo
  ! first position with key > iWalk
  iHi = N+1
  do while (iLo < iHi)
    iMid = (iLo+iHi)/2
    if (X%ICoup(k,Off+X%IOrd(k,Off+iMid)) <= iWalk) then
      iLo = iMid+1
    else
      iHi = iMid
    end if
  end do
  iLast = iLo-1

end subroutine Coup_Range

!> Sum of the coupling coefficients of the tuples (block of N at offset Off)
!> with ICoup(c,:) = iFrom and ICoup(3-c,:) = iTo, i.e. the transition
!> iFrom -> iTo (c = 1 canonical, c = 2 transposed). Uses the walk index
!> X%IOrd when present, else scans the block.
function Coup_Value(X,Off,N,c,iFrom,iTo) result(V)
  type(EXStruct), intent(in) :: X
  integer(kind=iwp), intent(in) :: Off, N, c, iFrom, iTo
  real(kind=wp) :: V
  integer(kind=iwp) :: i, j, ja, jb

  V = Zero
  if (N <= 0) return
  if (allocated(X%IOrd) .and. (N > nCoup_Lin)) then
    call Coup_Range(X,Off,N,c,iFrom,ja,jb)
    do j=ja,jb
      i = X%IOrd(c,Off+j)
      if (X%ICoup(3-c,Off+i) == iTo) V = V+X%VTab(X%ICoup(3,Off+i))
    end do
  else
    do i=1,N
      if ((X%ICoup(c,Off+i) == iFrom) .and. (X%ICoup(3-c,Off+i) == iTo)) V = V+X%VTab(X%ICoup(3,Off+i))
    end do
  end if

end function Coup_Value

!> Sum_K A1(iTo <- K) A2(K <- iFrom) for two blocks (Off2,N2,c2: the first
!> operator, applied to iFrom; Off1,N1,c1: the second, ending at iTo); c as
!> in Coup_Value. Only the transitions starting at iFrom are followed.
function Path_Value(X,Off2,N2,c2,iFrom,Off1,N1,c1,iTo) result(V)
  type(EXStruct), intent(in) :: X
  integer(kind=iwp), intent(in) :: Off2, N2, c2, iFrom, Off1, N1, c1, iTo
  real(kind=wp) :: V
  integer(kind=iwp) :: i2, j2, j2a, j2b

  V = Zero
  if ((N2 <= 0) .or. (N1 <= 0)) return
  if (allocated(X%IOrd) .and. (N2 > nCoup_Lin)) then
    call Coup_Range(X,Off2,N2,c2,iFrom,j2a,j2b)
    do j2=j2a,j2b
      i2 = X%IOrd(c2,Off2+j2)
      V = V+X%VTab(X%ICoup(3,Off2+i2))*Coup_Value(X,Off1,N1,c1,X%ICoup(3-c2,Off2+i2),iTo)
    end do
  else
    do i2=1,N2
      if (X%ICoup(c2,Off2+i2) /= iFrom) cycle
      V = V+X%VTab(X%ICoup(3,Off2+i2))*Coup_Value(X,Off1,N1,c1,X%ICoup(3-c2,Off2+i2),iTo)
    end do
  end if

end function Path_Value

!> Orbital occupations (per level) of the GUGA CSF I (after SG_MatEl_Init).
subroutine SG_CSF_Occ(I,OccVec)
  integer(kind=iwp), intent(in) :: I
  integer(kind=iwp), intent(out) :: OccVec(:)
  call CSF_Occupation(ME_LSM,I,OccVec)
end subroutine SG_CSF_Occ

!> Map graph level -> active orbital in symmetry-blocked order, the order
!> of the RASSCF integrals TU/TUVX. The levels run over the RAS (GAS)
!> partitions and within each over the symmetries; the orbitals of one
!> symmetry keep the order of their levels. For CAS the map is the identity.
subroutine SG_Level_Map(SGS,L2Orb)
  use stdalloc, only: mma_allocate, mma_deallocate
  type(SGStruct), intent(in) :: SGS
  integer(kind=iwp), intent(out) :: L2Orb(SGS%nLev)
  integer(kind=iwp) :: iLev, iSym, iOrb, n
  integer(kind=iwp), allocatable :: nOff(:)
  call mma_allocate(nOff,SGS%nSym,Label='nOff')
  nOff(:) = 0
  do iLev=1,SGS%nLev
    nOff(SGS%ISm(iLev)) = nOff(SGS%ISm(iLev))+1
  end do
  iOrb = 0
  do iSym=1,SGS%nSym
    n = nOff(iSym)
    nOff(iSym) = iOrb
    iOrb = iOrb+n
  end do
  do iLev=1,SGS%nLev
    nOff(SGS%ISm(iLev)) = nOff(SGS%ISm(iLev))+1
    L2Orb(iLev) = nOff(SGS%ISm(iLev))
  end do
  call mma_deallocate(nOff)
end subroutine SG_Level_Map

!> TU/TUVX (triangular, symmetry-blocked active orbitals) reordered to the
!> graph levels, for the codes that index the integrals by level
!> (MatEl_H_Direct, the diagonal codes). Identity for CAS.
subroutine SG_Level_Ints(SGS,TU,nTU,TUVX,nTUVX,TU_L,TUVX_L)
  use stdalloc, only: mma_allocate, mma_deallocate
  type(SGStruct), intent(in) :: SGS
  integer(kind=iwp), intent(in) :: nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: TU_L(nTU), TUVX_L(nTUVX)
  integer(kind=iwp) :: p, q, r, s, pq, rs, PQo, RSo
  integer(kind=iwp), allocatable :: L2Orb(:)
  call mma_allocate(L2Orb,SGS%nLev,Label='L2Orb')
  call SG_Level_Map(SGS,L2Orb)
  TU_L(:) = Zero
  TUVX_L(:) = Zero
  do p=1,SGS%nLev
    do q=1,p
      pq = Tri_SG(p,q)
      PQo = Tri_SG(L2Orb(p),L2Orb(q))
      TU_L(pq) = TU(PQo)
      do r=1,SGS%nLev
        do s=1,r
          rs = Tri_SG(r,s)
          if (rs > pq) cycle
          RSo = Tri_SG(L2Orb(r),L2Orb(s))
          TUVX_L(Tri_SG(pq,rs)) = TUVX(Tri_SG(PQo,RSo))
        end do
      end do
    end do
  end do
  call mma_deallocate(L2Orb)
end subroutine SG_Level_Ints

! Index heapsort of the multi-word keys K(1:nW,1:n), increasing (word nW
! most significant), ties by index.
subroutine Sort_Keys(K,nW,n,Ord)
  integer(kind=iwp), intent(in) :: nW, n, K(nW,*)
  integer(kind=iwp), intent(out) :: Ord(*)
  integer(kind=iwp) :: i, iEnd, t
  do i=1,n
    Ord(i) = i
  end do
  if (n < 2) return
  do i=n/2,1,-1
    call Sift(i,n)
  end do
  do iEnd=n,2,-1
    t = Ord(1)
    Ord(1) = Ord(iEnd)
    Ord(iEnd) = t
    call Sift(1,iEnd-1)
  end do
contains
  subroutine Sift(iStart,iLast)
    integer(kind=iwp), intent(in) :: iStart, iLast
    integer(kind=iwp) :: r, c, s
    r = iStart
    do
      c = 2*r
      if (c > iLast) exit
      if (c < iLast) then
        if (Key_Less(Ord(c),Ord(c+1))) c = c+1
      end if
      if (.not. Key_Less(Ord(r),Ord(c))) exit
      s = Ord(r)
      Ord(r) = Ord(c)
      Ord(c) = s
      r = c
    end do
  end subroutine Sift
  function Key_Less(a,b) result(L)
    integer(kind=iwp), intent(in) :: a, b
    logical(kind=iwp) :: L
    integer(kind=iwp) :: w
    do w=nW,1,-1
      if (K(w,a) /= K(w,b)) then
        L = K(w,a) < K(w,b)
        return
      end if
    end do
    L = a < b
  end function Key_Less
end subroutine Sort_Keys

end module SGUGA_Core
