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
! The core of SGUGA: types, state, graph and coupling tables, E_pq|Psi>.
! Default public: the other SGUGA_* modules use its entities; the
! interface of the library is the module SGUGA (sguga.F90).
module SGUGA_Core

use Molcas, only: MxLev, MxSym
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6
use CSF_Graph, only: SG_Graph, SG_Walks, nPack, LTAB, ATAB, BTAB, MkSGUGA, MkCOT, CSFCOUNT, Mk_IOW

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
   integer(kind=iwp), Allocatable :: OCReverse(:,:,:)

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
! DRT: Final compact distinct-row table with columns LTAB,NTAB,ATAB,BTAB,CTAB.
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
! The graph (SG_Graph of CSF_Graph) and the parts of SGUGA (CSF/SD codec)
type, extends(SG_Graph) :: SGStruct
  type(CDCInfo), allocatable :: CDCInfo(:,:)
  integer(kind=iwp), Allocatable :: OC(:)
  type(SDReoInfo), allocatable :: SDREO(:)
end type SGStruct

! the DRT column indices LTAB, NTAB, ATAB, BTAB, CTAB come from CSF_Graph

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
! The walks (SG_Walks of CSF_Graph) and the segment tables of SGUGA
type, extends(SG_Walks) :: CIStruct
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
  integer(kind=iwp), allocatable:: I1list(:), I2list(:)
  real(kind=wp), allocatable:: XList(:)
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
  integer(kind=iwp) :: nClass    = 0   ! number of topological classes
  integer(kind=iwp) :: nOpenBand = 0   ! number of open-family bands
  integer(kind=iwp) :: nTrans    = 0   ! total number of transitions

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
integer(kind=iwp), parameter :: iSegWUpBeg = nSegBase + 1
integer(kind=iwp), parameter :: iSegWUpEnd = nSegBase + 4
integer(kind=iwp), parameter :: iSegWLoBeg = nSegBase + 5
integer(kind=iwp), parameter :: iSegWLoEnd = nSegBase + 8

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

integer(kind=iwp), parameter ::                                                                                                    &
                                ITVPT(nSegCat) =  &
[ 0, 0, 0, 0,  0, 0, 0, 0,  1, 1, 1, 1, 1,  2, 2, 2, 2, 2,  1, 1, 2, 2,  3, 3, 3, 3,  0, 0, 0, 0,  3, 3, 3, 3, &
  4, 4, 4, 4,  5, 5, 6, 6,  4, 4, 4, 4, 4, 4,  5, 5, 5, 5, 5,  6, 6, 6, 6, 6,  4, 4, 5, 6],&
                                IBVPT(nSegCat) =  &
[ 0, 0, 0, 0,  1, 1, 2, 2,  1, 1, 2, 1, 1,  2, 2, 1, 2, 2,  3, 3, 3, 3,  3, 3, 3, 3,  0, 0, 0, 0,  3, 3, 3, 3, &
  1, 2, 2, 1,  1, 1, 2, 2,  4, 4, 4, 4, 6, 5,  5, 5, 5, 5, 4,  6, 6, 6, 6, 4,  3, 3, 3, 3],&
                                IC1(nSegCat)   =  &
[ 0, 1, 2, 3,  0, 2, 0, 1,  0, 1, 1, 2, 3,  0, 1, 2, 2, 3,  1, 3, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, &
  0, 0, 1, 2,  1, 0, 2, 0,  0, 1, 2, 3, 1, 2,  0, 1, 2, 3, 1,  0, 1, 2, 3, 2,  1, 2, 1, 2],&
                                IC2(nSegCat)   =  &
[ 0, 1, 2, 3,  1, 3, 2, 3,  0, 1, 2, 2, 3,  0, 1, 1, 2, 3,  0, 2, 0, 1,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, &
  1, 2, 3, 3,  3, 2, 3, 1,  0, 1, 2, 3, 2, 1,  0, 1, 2, 3, 2,  0, 1, 2, 3, 1,  1, 2, 2, 1],&
                                ISVC(nSegCat)  =  &
[ 1, 1, 1, 1,  1, 7, 8, 4,  1, 2, 9,10, 2,  1, 2,11,12, 2,  1, 5, 6, 3,  1, 1, 1, 1,  0, 1, 1,13,  0, 1, 1,13, &
  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0]

integer(kind=iwp), parameter :: TR_WALK  = 1
integer(kind=iwp), parameter :: TR_OPEN  = 2
integer(kind=iwp), parameter :: TR_MID   = 4
integer(kind=iwp), parameter :: TR_CLOSE = 8
integer(kind=iwp), parameter :: TR_TAIL  = 16
integer(kind=iwp), parameter :: TR_WEIGHT = 32

!#define _TEST_CSFDET_

! nPack (cases per integer of the packed walks) comes from CSF_Graph


Type (SGStruct), target :: SGS(4)
Type (CIStruct) :: CIS(4)
Type (EXStruct) :: EXS(4)
Logical(kind=iwp) :: State_is_used(4)=.False.


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


subroutine SG_Init(iState,nSym,nActEl,iSpin,                      &
                   nRas,nRasEl,nRsPrt,                            &
                   xLevel,xL2Act,xNLEV,xNSM,xnRasMx)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt),nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xnLev, xNSM(MxLev), xnRasMx(nRsPrt)

  type(TRStruct) :: TRS

  call SG_Init_Light(iState,nSym,nActEl,iSpin,        &
                      nRas,nRasEl,nRsPrt,             &
                      xLevel,xL2Act,xnLev,xNSM,xnRasMx=xnRasMx)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

  call MKMAW(SGS(iState))

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

  call MKSEG(SGS(iState),CIS(iState),EXS(iState))

    !Create the transition infrastructure.
  call MkTrans(SGS(iState),CIS(iState),TRS)

    ! Count coupling coefficients in compressed blocks indexed by excitation/operator type, symmetry and midvertex.
  call MkNRCOUP(SGS(iState),CIS(iState),EXS(iState),TRS)

    ! Explicitly generates the compressed coupling tuples '(left walk, right walk, value index)' and compacts repeated
    ! numerical values into 'VTab'.
  call MKCOUP(SGS(iState),CIS(iState),EXS(iState),TRS)

  call Trans_Free(TRS)

    ! Generate flat tables for sg_Epq_Psi_flat

  call MK_FLAT_Epq(SGS(iState))

    ! Generate the Configuration-to-Determinant Coefficients
    ! Generate Occupation Configurations a la Lucia.
    ! Temporarily just enabled for CAS
    ! (CDCInfo is needed by SG_ReOrd. The determinant reordering SDREO is
    ! only needed by the SD/CSF codec and is built on first use by
    ! SG_Mk_SDREO: building it here made e.g. a CAS(12,12) RASSCF about 75%
    ! slower.)

  if (SGS(iState)%nRsPrt == 1) then
     call Mk_CDC(SGS(iState))
  end if

end subroutine SG_Init

!> Build the determinant reordering (SDREO) for state iState, if not yet
!> done. It is needed only by the SD/CSF codec (CASVB with Faroald, and the
!> MOLCAS_CI_VERIFY=SIGMA checks). Only available for CAS, where SG_Init has
!> already built the CDCInfo tables it relies on.
subroutine SG_Mk_SDREO(iState)
  integer(kind=iwp), intent(in) :: iState
  if (allocated(SGS(iState)%SDREO)) return
  if (SGS(iState)%nRsPrt /= 1) then
     write(u6,*) 'SG_Mk_SDREO: the SD/CSF codec is only available for CASSCF'
     call Abend()
  end if
  if (.not. allocated(SGS(iState)%CDCInfo)) call Mk_CDC(SGS(iState))
  call Mk_SDREO_SGUGA(SGS(iState))
end subroutine SG_Mk_SDREO

subroutine SG_Init_Light(iState,nSym,nActEl,iSpin,              &
                          nRas,nRasEl,nRsPrt,                     &
                          xLevel,xL2Act,xNLEV,xNSM,Do_MkSGUGA,xnRasMx)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xNLEV, xNSM(MxLev), xnRasMx(nRsPrt)
  logical(kind=iwp), optional, intent(in) :: Do_MkSGUGA
  integer(kind=iwp) :: iSym

  If (State_is_used(iState)) Then
     Write (u6,*) 'SG_Init_Light: State is already initiated'
     Write (u6,*) 'iState: ',iState
     Call abend()
  End If

  SGS(istate)%IFRAS=0
  If (nRsPrt==3) SGS(istate)%IFRAS=1
  Do iSym = 1, nSym
     If (Sum(nRas(iSym,1:nRsPrt))/=0) SGS(istate)%IFRAS=SGS(istate)%IFRAS+1
  End Do
  SGS(istate)%nRasEL(1:nRsPrt)=nRasEl(1:nRsPrt)
  SGS(istate)%nRas(1:MxSym,1:nRsPrt)=nRas(1:MxSym,1:nRsPrt)
  SGS(istate)%nRsPrt=nRsPrt
  ! GAS: upper bounds on the electrons in the partitions 1..k
  SGS(istate)%nRasMx(:)=-1
  if (present(xnRasMx)) SGS(istate)%nRasMx(1:nRsPrt)=xnRasMx(1:nRsPrt)

  ! Make sure that we start from a clean slate.
  call SG_Free(istate)
  State_is_used(iState)=.True.

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

  SGS(istate)%nSym=nSym
  SGS(istate)%iSpin=iSpin
  SGS(istate)%nActEl=nActEl

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
    if (Do_MkSGUGA) Then
       call MkSGUGA(SGS(istate),CIS(istate))

  else
       SGS(istate)%iSpin = 0
       SGS(istate)%nActEl = 0

       ! INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
       call mma_allocate(CIS(istate)%NCSF,SGS(istate)%nSym,Label='CIS%NCSF')
       CIS(istate)%NCSF(:) = 0
    endif
  else
    call MkSGUGA(SGS(istate),CIS(istate))
  end if

  if (NActEl == 0) CIS(istate)%NCSF(:) = 1

end subroutine SG_Init_Light

subroutine SG_Free(iState)
! PURPOSE: FREE THE SGUGA TABLES

integer(kind=iwp), intent(in):: iState

integer(kind=iwp) :: iOpen, mult, iRefSm

If (.Not.State_is_used(iState)) Return
if (ASTA_iState == iState) call ASTA_Free()
if (TRP_iState == iState) call TRP_Free()

Associate (SGS=>SGS(iState),CIS=>CIS(iState),EXS=>EXS(iState))

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

Call mma_deallocate(SGS%OC,safe='*')
if (allocated(SGS%CDCInfo)) then


   do iOpen = lbound(SGS%CDCInfo,1), ubound(SGS%CDCInfo,1)

      do mult = lbound(SGS%CDCInfo,2), ubound(SGS%CDCInfo,2)

         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOpenOCP,safe='*')
         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOCP,safe='*')
         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%DOC,safe='*')
         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%SOC,safe='*')
         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%OCReverse,safe='*')
         call mma_deallocate(SGS%CDCInfo(iOpen,mult)%SOCSym,safe='*')
         Call mma_deallocate(SGS%CDCInfo(iOpen,mult)%nOCSym,safe='*')
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

endif

if (allocated(SGS%SDREO)) then
   do iRefSm = 1, size(SGS%SDREO)
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

End Associate

State_is_used(iState)=.False.

end subroutine SG_Free

subroutine Mk_CDC(SGS)


type(SGStruct), intent(inout) :: SGS

integer(kind=iwp) :: iOpen, mult, ms2, nDet, nCSF, MinOpen, MaxOpen, iSym, iDOC
integer(kind=iwp) :: nClosed, nClosedOCP, nSOC, iSOC, nDOC, iOC
integer(kind=iwp) :: NoNew, nConf, nElec, nDetSC
Integer :: i,j
Integer :: iStart
integer(kind=iwp), allocatable :: nSOCSym(:)
integer(kind=iwp), allocatable :: detspin(:,:)
integer(kind=iwp), allocatable :: twos(:,:)

Integer :: DOCOcc(SGS%nLev)
Integer :: SOCOcc(SGS%nLev)
Integer(kind=iwp) :: nDOCList
Integer(kind=iwp) :: nSOCList

allocate(SGS%CDCInfo(0:SGS%nLev,SGS%iSpin:SGS%iSpin))
SGS%CDCInfo(:,:)%nDet = 0


mult = SGS%iSpin

MinOpen = mult - 1
MaxOpen = Min(SGS%nActEl, 2*SGS%nLev-SGS%nActel)
if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen - 1

call mma_allocate( nSOCSym,SGS%nSym,Label='nSOCSym')
do iOpen = MinOpen, MaxOpen, 2

   ms2 = mult-1

   call count_detspin(iOpen,ms2,nDet)
   call count_twos(iOpen,mult,nCSF)

   SGS%CDCInfo(iOpen,mult)%nDet = nDet
   SGS%CDCInfo(iOpen,mult)%nCSF = nCSF

   nClosed = (SGS%nActEl-iOpen)/2
   SGS%CDCInfo(iOpen,Mult)%nClosed    = nClosed

   nClosedOCP = binomial(nClosed,SGS%nLev)
   SGS%CDCInfo(iOpen,Mult)%nClosedOCP = nClosedOCP

   nSOC = binomial(iOpen,SGS%nLev)

   if (nClosed > 0) then
      call mma_allocate( SGS%CDCInfo(iOpen,mult)%DOC, &
                         nClosed, &
                         nClosedOCP, &
                         Label='DOC')
      iDOC = 0
      call Make_DOCs(1,1)
      if (SG_Checks) then  ! sanity check, only when verifying CI developments
      if (iDOC /= SGS%CDCInfo(iOpen,mult)%nClosedOCP) then
         write(u6,*) 'DOC generation failed'
         write(u6,*) 'Generated = ',iDOC
         write(u6,*) 'Expected  = ',SGS%CDCInfo(iOpen,mult)%nClosedOCP
         call Abend()
      endif
      end if
   endif

   nSOCSym(:) = 0

   call mma_allocate( SGS%CDCInfo(iOpen,mult)%nOCP, &
                      SGS%nSym, &
                      Label='nOCP')

   call mma_allocate( SGS%CDCInfo(iOpen,mult)%nOpenOCP, &
                      SGS%nSym, &
                      Label='nOpenOCP')
   SGS%CDCInfo(iOpen,mult)%nOpenOCP(:) = 0

   call Count_OpenShell_OCPs(SGS%ISm,      &
                             SGS%nLev,     &
                             iOpen,        &
                             SGS%nSym,     &
                             SGS%CDCInfo(iOpen,mult)%nOpenOCP )

   SGS%CDCInfo(iOpen,mult)%nOCP(:) = SGS%CDCInfo(iOpen,mult)%nOpenOCP(:) * SGS%CDCInfo(iOpen,mult)%nClosedOCP

   if (iOpen > 0) then
      call mma_allocate( SGS%CDCInfo(iOpen,mult)%SOC, &
                         iOpen, &
                         nSOC, &
                         Label='SOC')

      call mma_allocate( SGS%CDCInfo(iOpen,mult)%SOCSym, &
                         nSOC, &
                         Label='SOCSym')
      iSOC = 0
      nSOCSym(:) = 0
      call Make_SOCs(1,1,1)

      if (SG_Checks) then  ! sanity check, only when verifying CI developments
      if (iSOC /= nSOC) then
         write(u6,*) 'SOC generation failed'
         write(u6,*) 'Generated = ',iSOC
         write(u6,*) 'Expected  = ',nSOC
         call Abend()
      endif
      do iSym=1,SGS%nSym
         if (nSOCSym(iSym) /= SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym)) then
            write(u6,*) 'SOC symmetry count mismatch'
            write(u6,*) 'ISYM     = ',iSym
            write(u6,*) 'GENERATED= ',nSOCSym(iSym)
            write(u6,*) 'EXPECTED = ',SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym)
            call Abend()
         endif
      enddo
      end if

   endif

   ! Allocate and initialize the reverse list: iDOC,iSOC -> iOC

   If (Allocated(SGS%CDCInfo(iOpen,mult)%DOC)) Then
      nDOCList = Size(SGS%CDCInfo(iOpen,mult)%DOC,2)
   Else
      nDOCList = 1
   End If
   If (Allocated(SGS%CDCInfo(iOpen,mult)%SOC)) Then
      nSOCList = Size(SGS%CDCInfo(iOpen,mult)%SOC,2)
   Else
      nSOCList = 1
   End If
   Call mma_allocate(SGS%CDCInfo(iOpen,mult)%OCReverse, &
                     nDOCList,nSOCList,SGS%nSym,         &
                     Label='OCReverse')

   SGS%CDCInfo(iOpen,mult)%OCReverse = 0_iwp

   if (SGS%CDCInfo(iOpen,mult)%nDet > 0 .and. &
       SGS%CDCInfo(iOpen,mult)%nCSF > 0) then

      call generate_detspin(iOpen, mult-1, detspin, nDet)

      call mma_allocate( &
         SGS%CDCInfo(iOpen,mult)%DetSpin, &
         [lbound(detspin,1),ubound(detspin,1)], &
         [lbound(detspin,2),ubound(detspin,2)], &
         Label='DetSpin')

      SGS%CDCInfo(iOpen,mult)%DetSpin(:,:) = detspin(:,:)

      call generate_twos(iOpen, mult, twos, nCSF)
      call mma_allocate( &
         SGS%CDCInfo(iOpen,mult)%Twos, &
         [lbound(twos,1),ubound(twos,1)], &
         [lbound(twos,2),ubound(twos,2)], &
         Label='Twos')

      SGS%CDCInfo(iOpen,mult)%Twos(:,:)    = twos(:,:)

      call mma_allocate( &
         SGS%CDCInfo(iOpen,mult)%CDC, &
         nDet, nCSF, &
         Label='CDC')

      call sg_csfdet( &
         iOpen,  &
         nDet,   &
         detspin,&
         nCSF,   &
         twos,   &
         SGS%CDCInfo(iOpen,mult)%CDC )

! keep this for the time being.
      if (mult == 1) then

         nDetSC = Max(1,nDet/2)
         SGS%CDCInfo(iOpen,mult)%nDetSC = nDetSC

         call mma_allocate(   &
            SGS%CDCInfo(iOpen,mult)%DTOC,   &
            nDetSC,   &
            nCSF,   &
            Label='DTOC')

         call mma_allocate(   &
            SGS%CDCInfo(iOpen,mult)%Pair,   &
            nDet,   &
            Label='Pair')

         if (nDet > 1) then
            call mk_dtoc_singlet( &
                 nDet, &
                 nCSF, &
                 SGS%CDCInfo(iOpen,mult)%DetSpin, &
                 SGS%CDCInfo(iOpen,mult)%CDC, &
                 SGS%CDCInfo(iOpen,mult)%Pair, &
                 SGS%CDCInfo(iOpen,mult)%DTOC )
         else if (nDet == 1) Then
            SGS%CDCInfo(iOpen,mult)%Pair(1) = 1
            SGS%CDCInfo(iOpen,mult)%DTOC(1,1) = SGS%CDCInfo(iOpen,mult)%CDC(1,1)
         else
            Write (u6,*) 'Mk_CDC: illegal nDet value: nDet<0'
            Call Abend()
         end if
      endif

      call mma_deallocate(detspin)
      call mma_deallocate(twos)

   endif


   if (SG_Checks) then  ! sanity check, only when verifying CI developments
   do iSym=1,SGS%nSym
      if (SGS%CDCInfo(iOpen,mult)%nOCP(iSym) /= &
          SGS%CDCInfo(iOpen,mult)%nOpenOCP(iSym) * &
          SGS%CDCInfo(iOpen,mult)%nClosedOCP) then
         call Abend()
      endif
   enddo

   if (sum(SGS%CDCInfo(iOpen,mult)%nOpenOCP) /= binomial(iOpen,SGS%nLev)) then
      call Abend()
   endif
   end if

   ! Compute the number of orbital occupatation configurations per symmetry.
   ! Here the order doesn't matter, we will deal with that in the next loop
   ! where we copy the Lucia order.

   Call mma_allocate(SGS%CDCInfo(iOpen,Mult)%nOCSym,SGS%nSym,Label='nOCSym')

   SGS%CDCInfo(iOpen,mult)%nOCSym(:) = &
      SGS%CDCInfo(iOpen,mult)%nClosedOCP * &
      SGS%CDCInfo(iOpen,mult)%nOpenOCP(:)

   call mma_allocate( SGS%CDCInfo(iOpen,mult)%OCDOC, &
                      maxval(SGS%CDCInfo(iOpen,mult)%nOCSym), &
                      SGS%nSym, &
                      Label='OCDOC' )

   call mma_allocate( SGS%CDCInfo(iOpen,mult)%OCSOC, &
                      maxval(SGS%CDCInfo(iOpen,mult)%nOCSym), &
                      SGS%nSym, &
                      Label='OCSOC' )

   SGS%CDCInfo(iOpen,mult)%OCDOC(:,:) = 0
   SGS%CDCInfo(iOpen,mult)%OCSOC(:,:) = 0
end do

call mma_deallocate(nSOCSym)

Call mma_allocate( SGS%OC, SGS%nActel, Label='OC' )

! Now loop over the Orbital Configurations a la Lucia.
! For each OC identify the symmtry and store away the iDOC and iSOC indicies.

! Reset and use now as a counter
do iOpen = MinOpen, MaxOpen, 2
   SGS%CDCInfo(iOpen,mult)%nOCSym(:) = 0
end do

nElec = SGS%nActEl
nConf=0
NoNew = 0
Call Nxt_OC(SGS%OC,nElec,SGS%nLev,1,NoNew)

Do While (NoNew == 0)
   nConf=nConf+1

   Call Decode_OC(SGS%OC,nElec, &
                  nDOC,DOCOcc,  &
                  nSOC,SOCOcc)

   If (nSOC < MinOpen .or. nSOC > MaxOpen) Then
      Call Nxt_OC(SGS%OC,nElec,SGS%nLev,0,NoNew)
      Cycle
   End If

   iDOC = 1
   Do i=1,nDOC
      If (i == 1) Then
         iStart = 1
      Else
         iStart = DOCOcc(i-1) + 1
      End If
      Do j=iStart,DOCOcc(i)-1
         iDOC = iDOC + Binomial(nDOC-i,SGS%nLev-j)
      End Do
   End Do

   if (SG_Checks) then  ! sanity check, only when verifying CI developments
   If (nDOC > 0) Then
      Do i=1,nDOC
         If (DOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%DOC(i,iDOC)) Then
            Write(u6,*) 'DOC ranking mismatch'
            Write(u6,'(A,20I4)') 'DOC OCC = ',DOCOcc(1:nDOC)
            Write(u6,'(A,I6)')   'iDOC    = ',iDOC
            Write(u6,'(A,20I4)') 'DOC TAB = ', &
                                 SGS%CDCInfo(nSOC,mult)%DOC(1:nDOC,iDOC)
            Call Abend()
         End If
      End Do
   End If
   end if

   iSOC = 1
   Do i=1,nSOC

      If (i == 1) Then
         iStart = 1
      Else
         iStart = SOCOcc(i-1) + 1
      End If

      Do j=iStart,SOCOcc(i)-1
         iSOC = iSOC + Binomial(nSOC-i,SGS%nLev-j)
      End Do

   End Do

   if (SG_Checks) then  ! sanity check, only when verifying CI developments
   If (nSOC > 0) Then
      Do i=1,nSOC
         If (SOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%SOC(i,iSOC)) Then
            Write(u6,*) 'SOC ranking mismatch'
            Write(u6,'(A,20I4)') 'SOC OCC = ',SOCOcc(1:nSOC)
            Write(u6,'(A,I6)')   'iSOC    = ',iSOC
            Write(u6,'(A,20I4)') 'SOC TAB = ', &
                              SGS%CDCInfo(nSOC,mult)%SOC(1:nSOC,iSOC)
            Call Abend()
         End If
      End Do
   End If
!  Write(u6,'(A,4I8)') 'OC:', nConf,nSOC,iDOC,iSOC
   end if

   if (nSOC == 0) then
      iSym = 1
      iSOC = 1
   else
      iSym = SOC_Sym(SGS%ISm,nSOC,SOCOcc)

      if (SG_Checks) then  ! sanity check, only when verifying CI developments
      do i=1,nSOC
         if (SOCOcc(i) /= SGS%CDCInfo(nSOC,mult)%SOC(i,iSOC)) then
            write(u6,*) 'SOC ranking mismatch'
            write(u6,'(A,20I4)') 'SOC OCC = ',SOCOcc(1:nSOC)
            write(u6,'(A,I6)')   'iSOC    = ',iSOC
            write(u6,'(A,20I4)') 'SOC TAB = ', &
                              SGS%CDCInfo(nSOC,mult)%SOC(1:nSOC,iSOC)
            call Abend()
         end if
      end do

      if (iSym /= SGS%CDCInfo(nSOC,mult)%SOCSym(iSOC)) then
         write(u6,*) 'SOC symmetry mismatch'
         call Abend()
      end if
      end if
   end if

   SGS%CDCInfo(nSOC,mult)%nOCSym(iSym) = SGS%CDCInfo(nSOC,mult)%nOCSym(iSym) + 1

   iOC = SGS%CDCInfo(nSOC,mult)%nOCSym(iSym)

   SGS%CDCInfo(nSOC,mult)%OCDOC(iOC,iSym) = iDOC
   SGS%CDCInfo(nSOC,mult)%OCSOC(iOC,iSym) = iSOC

   SGS%CDCInfo(nSOC,mult)%OCReverse(iDOC,iSOC,iSym) = iOC

   Call Nxt_OC(SGS%OC,nElec,SGS%nLev,0,NoNew)

End Do

contains

subroutine mk_dtoc_singlet(nDet,nCSF,DetSpin,CDC,Pair,DTOC)

use Definitions, only: wp, iwp
implicit none

integer(iwp), intent(in)  :: nDet
integer(iwp), intent(in)  :: nCSF

integer(iwp), intent(in)  :: DetSpin(:,:)

real(wp),    intent(in)   :: CDC(nDet,nCSF)

integer(iwp), intent(out) :: Pair(nDet)

real(wp),    intent(out)  :: DTOC(nDet/2,nCSF)

integer(iwp) :: iDet
integer(iwp) :: jDet
integer(iwp) :: iRow

Pair(:) = 0

!
! Determine spin-reversal partner for every determinant
!
do iDet=1,nDet

   do jDet=1,nDet

      if (all(DetSpin(:,iDet) == 1-DetSpin(:,jDet))) then
         Pair(iDet) = jDet
         exit
      end if

   end do

   if (Pair(iDet) == 0) then
      write(6,*) 'MK_DTOC_SINGLET: partner not found'
      write(6,*) 'IDET =',iDet
      call Abend()
   end if

end do

!
! Generate Lucia DTOC matrix.
!
! Keep one determinant from every spin-reversal pair
! and renormalize by sqrt(2).
!
iRow = 0

do iDet=1,nDet

   ! Lucia/SPNCOM ordering:
   ! keep determinants with first spin = alpha

   if (DetSpin(1,iDet) == 1) then

      iRow = iRow + 1

      DTOC(iRow,:) = sqrt(2.0_wp) * CDC(iDet,:)

   end if

end do

if (iRow /= Max(1,nDet/2)) then
   write(6,*) 'MK_DTOC_SINGLET: incorrect number of pairs'
   write(6,*) 'Generated =',iRow
   write(6,*) 'Expected  =',Max(1,nDet/2)
   call Abend()
end if

end subroutine mk_dtoc_singlet

recursive subroutine Make_SOCs(iPos,iStart,iSym)

   integer(iwp), intent(in) :: iPos
   integer(iwp), intent(in) :: iStart
   integer(iwp), intent(in) :: iSym

   integer(iwp) :: i
   integer(iwp), save :: work(MxLev)

   if (iPos > iOpen) then

      iSOC = iSOC + 1

      SGS%CDCInfo(iOpen,mult)%SOC(:,iSOC) = work(1:iOpen)
      SGS%CDCInfo(iOpen,mult)%SOCSym(iSOC) = iSym

      nSOCSym(iSym) = nSOCSym(iSym) + 1

#ifdef _DEBUGPRINT_
      write(u6,*) 'ISm = ',SGS%ISm(1:SGS%nLev)
      write(u6,*) 'NOPEN=',iOpen, &
                  ' SOC=',iSOC, &
                  ' OCC=',work(1:iOpen), &
                  ' SYM=',iSym
#endif
      return
   endif

   do i=iStart,SGS%nLev-(iOpen-iPos)

      work(iPos) = i

      call Make_SOCs( iPos+1, &
                      i+1, &
                      Mul(iSym,SGS%ISm(i)) )

   enddo

end subroutine Make_SOCs
recursive subroutine Make_DOCs(iPos,iStart)

   integer(iwp), intent(in) :: iPos
   integer(iwp), intent(in) :: iStart

   integer(iwp) :: i
   integer(iwp), save :: work(MxLev)
!
!  Completed DOC
!
   if (iPos > nClosed) then
      iDOC = iDOC + 1
      SGS%CDCInfo(iOpen,mult)%DOC(:,iDOC) = work(1:nClosed)
#ifdef _DEBUGPRINT_
      write(u6,'(A,I4,A,100I4)') 'NOPEN=',iOpen, &
                                 ' DOC=',iDOC, &
                                 SGS%CDCInfo(iOpen,mult)%DOC(:,iDOC)
#endif
      return
   endif
!
!  Lexical combinations
!
   do i=iStart,SGS%nLev-(nClosed-iPos)
      work(iPos) = i
      call Make_DOCs(iPos+1,i+1)
   enddo

end subroutine Make_DOCs

subroutine count_detspin(nopen,ms2,ndet)

use Definitions, only: iwp

implicit none

integer(iwp), intent(in)  :: nopen
integer(iwp), intent(in)  :: ms2
integer(iwp), intent(out) :: ndet

integer(iwp) :: ibit
integer(iwp) :: mx

ndet = 0
mx   = 2**nopen

do ibit=0,mx-1
   if (2*popcnt(ibit)-nopen == ms2) ndet = ndet + 1
enddo

end subroutine count_detspin

subroutine count_twos(nopen,mult,ncsf)

use Definitions, only: iwp

implicit none

integer(iwp), intent(in)  :: nopen
integer(iwp), intent(in)  :: mult
integer(iwp), intent(out) :: ncsf

integer(iwp) :: ibit
integer(iwp) :: i
integer(iwp) :: ms2l
integer(iwp) :: mx

logical :: upper

ncsf = 0
mx   = 2**nopen

do ibit=0,mx-1

   ms2l  = 0
   upper = .true.

   do i=1,nopen

      if (ibits(ibit,i-1,1) == 1) then
         ms2l = ms2l + 1
      else
         ms2l = ms2l - 1
      endif

      if (ms2l < 0) then
         upper = .false.
         exit
      endif

   enddo

   if (.not.upper) cycle
   if (ms2l /= mult-1) cycle

   ncsf = ncsf + 1

enddo

end subroutine count_twos

end subroutine Mk_CDC

!***********************************************************************
! Construct the SD reordering information required by the SGUGA-based
! CAS SD-CSF codec.
!
! Initial version:
!   Loop over the alpha and beta string symmetry pairs compatible with
!   the requested total spatial symmetry.
!***********************************************************************

subroutine Mk_SDREO_SGUGA(SGS)

  use Definitions,   only: iwp
  use Symmetry_Info, only: Mul

  implicit none

  type(SGStruct), intent(inout) :: SGS

  integer(kind=iwp) :: iRefSm
  integer(kind=iwp) :: iASm
  integer(kind=iwp) :: iBSm
  integer(kind=iwp) :: nSym
  integer(kind=iwp) :: iDet

  integer(kind=iwp) :: nAlpha
  integer(kind=iwp) :: nBeta
  integer(kind=iwp) :: nAlphaS
  integer(kind=iwp) :: nBetaS
  integer(kind=iwp) :: iAlphaS
  integer(kind=iwp) :: iBetaS

  integer(kind=iwp) :: nAlphaSSym(SGS%nSym)
  integer(kind=iwp) :: nBetaSSym(SGS%nSym)
  integer(kind=iwp), allocatable :: AlphaStrings(:,:)
  integer(kind=iwp), allocatable :: BetaStrings(:,:)
  Integer(kind=iwp) :: iDetOC(SGS%nActEl)
  Integer(kind=iwp) :: iDetSP(SGS%nActEl)
  Integer(kind=iwp) :: iSign

  Integer(kind=iwp) :: nOpen, i, iOpen, mult, iDOC, iSOC
  Integer(kind=iwp) :: nDouble
  Integer(kind=iwp) :: DOCOcc(SGS%nActel)
  Integer(kind=iwp) :: SOCOcc(SGS%nActel)
  Integer(kind=iwp) :: OpenSP(SGS%nActel)

  Integer(kind=iwp) :: iOC
  Integer(kind=iwp) :: nOC
  Integer(kind=iwp) :: iSpinPattern
  Integer(kind=iwp) :: iOpenOffset(0:SGS%nLev)
  Integer(kind=iwp) :: jOpen
  Integer(kind=iwp) :: MinOpen
  Integer(kind=iwp) :: MaxOpen
  Integer(kind=iwp) :: nOCSector
  Integer(kind=iwp) :: nDetSector
  Integer(kind=iwp) :: iAdr, nDetTotal
  Integer(kind=iwp) :: iAdrSDConfOrder

  nSym = SGS%nSym
  mult = SGS%iSpin

  !
  ! There is one SD reorder list for each total symmetry.
  !
  allocate(SGS%SDREO(1:nSym))

  ! Highest-Ms component:
  !
  !   nAlpha + nBeta = nActEl
  !   nAlpha - nBeta = 2*Ms = iSpin - 1
  !
  nAlpha = (SGS%nActEl + SGS%iSpin - 1_iwp)/2_iwp
  nBeta  = (SGS%nActEl - SGS%iSpin + 1_iwp)/2_iwp

  call Count_Strings_By_Symmetry(SGS,nAlpha,nAlphaSSym)
  call Count_Strings_By_Symmetry(SGS,nBeta,nBetaSSym)

  ! Outer loop over possible total spatial symmetries.
  do iRefSm = 1, nSym

     MinOpen = mult - 1_iwp

     MaxOpen = Min(SGS%nActEl,2_iwp*SGS%nLev-SGS%nActEl)
     If (Mod(MaxOpen-MinOpen,2_iwp) /= 0_iwp) MaxOpen = MaxOpen - 1_iwp

     iOpenOffset(:) = 0_iwp
     iAdr = 0_iwp
     Do jOpen=MinOpen,MaxOpen,2_iwp
        iOpenOffset(jOpen) = iAdr
        nOCSector = SGS%CDCInfo(jOpen,mult)%nOCSym(iRefSm)
        nDetSector = SGS%CDCInfo(jOpen,mult)%nDet

        iAdr = iAdr + nOCSector*nDetSector
     End Do
     nDetTotal = iAdr

     If (Allocated(SGS%SDREO(iRefSm)%List)) Deallocate(SGS%SDREO(iRefSm)%List)

     Call mma_allocate(SGS%SDREO(iRefSm)%List,nDetTotal,Label='List')

     SGS%SDREO(iRefSm)%List(:) = 0_iwp

     ! Loop over symmetry blocks in the same order as LUCIA.
     iDet = 0_iwp

     do iASm = 1, nSym

        iBSm = Mul(iASm,iRefSm)

        nAlphaS = nAlphaSSym(iASm)
        nBetaS  = nBetaSSym(iBSm)

        If (nAlphaS*nBetaS == 0) Cycle

        call mma_allocate(AlphaStrings,Max(1_iwp,nAlpha),nAlphaS,Label='AStrings')
        call mma_allocate(BetaStrings,Max(1_iwp,nBeta),nBetaS,Label='BStrings')

        AlphaStrings(:,:)=0_iwp
        BetaStrings(:,:)=0_iwp

        call GenerateStringSet(nAlpha,iASm,nAlphaS,AlphaStrings)
        call GenerateStringSet(nBeta, iBSm,nBetaS, BetaStrings)

        do iBetaS = 1, nBetaS
           do iAlphaS = 1, nAlphaS


              Call BuildOrderedDeterminant(AlphaStrings(:,iAlphaS),BetaStrings(:,iBetaS), &
                                           nAlpha,nBeta,iDetOC,iDetSP,iSign)

              Call ExtractOccupationPattern(iDetOC,iDetSP,SGS%nActEl, &
                                            DOCOcc,SOCOcc,OpenSP,nDouble,nOpen)

              iOpen = nOpen

if (SG_Checks) then  ! sanity check, only when verifying CI developments
If (2_iwp*nDouble+iOpen /= SGS%nActEl) Then
   Write (u6,*) 'Mk_SDREO_SGUGA_new: inconsistent occupation counts'
   Call Abend()
End If
If (Mod(iOpen,2_iwp) /= Mod(SGS%nActEl,2_iwp)) Then
   Write (u6,*) 'Mk_SDREO_SGUGA_new: invalid number of open orbitals'
   Call Abend()
End If

If (iOpen > 0_iwp) Then
   If (Count(OpenSP(1:iOpen) == 1_iwp) /= nAlpha-nDouble) Then
      Write (u6,*) 'Mk_SDREO_SGUGA_new: incorrect open alpha count'
      Call Abend()
   End If
   If (Count(OpenSP(1:iOpen) == -1_iwp) /= nBeta-nDouble) Then
      Write (u6,*) 'Mk_SDREO_SGUGA_new: incorrect open beta count'
      Call Abend()
   End If
Else
   If (nAlpha /= nDouble .or. nBeta /= nDouble) Then
      Write (u6,*) 'Mk_SDREO_SGUGA_new: invalid closed-shell determinant'
      Call Abend()
   End If
End If
If (nDouble /= 0_iwp .and. .NOT.Allocated(SGS%CDCInfo(nOpen,mult)%DOC) ) Then
   Write (u6,*) 'Mk_SDREO_SGUGA_new: ...%DOC not allocated while nDouble /= 0!'
   Call Abend()
End If
end if
              If (nDouble == 0_iwp) Then
                 iDOC = 1_iwp
              Else
                 Call FindDOCIndex(iDetOC,nAlpha+nBeta, &
                                   SGS%CDCInfo(nOpen,mult)%DOC,iDOC)
              End If

              If (nOpen == 0_iwp) Then
                 iSOC = 1_iwp
              Else
                 Call FindSOCIndex(iDetOC,nAlpha+nBeta, &
                                   SGS%CDCInfo(nOpen,mult)%SOC,iSOC)
              End If

              nOC = SGS%CDCInfo(nOpen,mult)%nOCSym(iRefSm)

              Call FindOCIndex(iDOC,                                               &
                               iSOC,                                               &
                               SGS%CDCInfo(nOpen,mult)%OCDOC(:,iRefSm),            &
                               SGS%CDCInfo(nOpen,mult)%OCSOC(:,iRefSm),            &
                               nOC,                                                &
                               iOC)


              if (SG_Checks) then  ! sanity check, only when verifying CI developments
              If (SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC) Then
                 Write (u6,'(A)') 'SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC'
                 Call Abend()
              End If
              end if

              Call FindSpinIndex(OpenSP,                                               &
                                 iOpen,                                                &
                                 SGS%CDCInfo(iOpen,mult)%DetSpin,                      &
                                 iSpinPattern)

              if (SG_Checks) then  ! sanity check, only when verifying CI developments
              Do i=1,iOpen
                 If (OpenSP(i) == 1_iwp) Then
                    If (SGS%CDCInfo(iOpen,mult)%DetSpin(i,iSpinPattern) /= 1_iwp) Then
                       Write (u6,'(A)') 'Mk_SDREO_SGUGA_New: alpha-spin mismatch'
                       Call Abend()
                    End If
                 Else
                    If (SGS%CDCInfo(iOpen,mult)%DetSpin(i,iSpinPattern) /= 0_iwp) Then
                       Write (u6,'(A)') 'Mk_SDREO_SGUGA_New: beta-spin mismatch'
                       Call Abend()
                    End If
                 End If
              End Do
              end if
              iAdrSDConfOrder = iOpenOffset(iOpen)                  &
                              + (iOC-1_iwp)*SGS%CDCInfo(iOpen,mult)%nDet + iSpinPattern
              iDet = iDet + 1_iwp
              if (SG_Checks) then  ! sanity check, only when verifying CI developments
              If (iAdrSDConfOrder < 1_iwp .or. &
                  iAdrSDConfOrder > nDetTotal) Then
                 Write(u6,'(A,2(1X,I0))') 'ADDRESS RANGE', iAdrSDConfOrder, nDetTotal
                 Write(u6,'(A)') 'Mk_SDREO_SGUGA_New: invalid configuration address'
                 Call Abend()
              End If
              If (SGS%SDREO(iRefSm)%List(iAdrSDConfOrder) /= 0_iwp) Then
                 Write(u6,'(A,3(1X,I0))') 'DUPLICATE SDREO', iRefSm, iAdrSDConfOrder, iDet
                 Write(u6,'(A)') 'Mk_SDREO_SGUGA: duplicate SDREO address'
                 Call Abend()
              End If
              end if

              SGS%SDREO(iRefSm)%List(iAdrSDConfOrder) = iSign*iDet

           end do
        end do

        call mma_deallocate(AlphaStrings)
        call mma_deallocate(BetaStrings)

     end do
     if (SG_Checks) then  ! sanity check, only when verifying CI developments
     If (iDet /= nDetTotal) Then
        Write(u6,'(A,3(1X,I0))') 'SDREO SIZE MISMATCH', iRefSm, iDet, nDetTotal
        Write(u6,'(A)') 'Mk_SDREO_SGUGA: incorrect determinant count'
        Call Abend()
     End If

     If (Any(SGS%SDREO(iRefSm)%List == 0_iwp)) Then
        Write(u6,'(A,1X,I0)') 'INCOMPLETE SDREO FOR SYMMETRY', iRefSm
        Write(u6,'(A)') 'Mk_SDREO_SGUGA: incomplete SDREO vector'
        Call Abend()
     End If
     end if

  end do

contains

Subroutine FindSpinIndex(OpenSP,iOpen,DetSpin,iSpinPattern)

Implicit None

Integer(kind=iwp), Intent(In)  :: iOpen
Integer(kind=iwp), Intent(In)  :: OpenSP(Max(1_iwp,iOpen))
Integer(kind=iwp), Intent(In)  :: DetSpin(:,:)
Integer(kind=iwp), Intent(Out) :: iSpinPattern

Integer(kind=iwp) :: i
Integer(kind=iwp) :: j
Integer(kind=iwp) :: nSpin
Integer(kind=iwp) :: nMatch
Integer(kind=iwp) :: SpinValue
Logical           :: Match

nSpin  = Size(DetSpin,2)
iSpinPattern  = 0_iwp
nMatch = 0_iwp

If (iOpen == 0_iwp) Then

   If (nSpin /= 1_iwp) Then
      Write (u6,'(A)') 'FindSpinIndex: invalid closed-shell spin table'
      Call Abend()
   End If

   iSpinPattern = 1_iwp
   Return

End If

If (Size(DetSpin,1) < iOpen) Then
   Write (u6,'(A)') 'FindSpinIndex: spin table has too few rows'
   Call Abend()
End If

SpinValue = -1_iwp
Do j=1,nSpin

   Match = .True.

   Do i=1,iOpen

      If (OpenSP(i) == 1_iwp) Then
         SpinValue = 1_iwp
      Else If (OpenSP(i) == -1_iwp) Then
         SpinValue = 0_iwp
      Else
         Write (u6,'(A)')  'FindSpinIndex: invalid open-shell spin'
         Call Abend()
      End If

      If (DetSpin(i,j) /= SpinValue) Then
         Match = .False.
         Exit
      End If

   End Do

   If (Match) Then
      iSpinPattern  = j
      nMatch = nMatch + 1_iwp
   End If

End Do

If (nMatch == 0_iwp) Then
   Write(u6,'(A)') 'FindSpinIndex: open-spin pattern not found'
   Write(u6,'(A,2(1X,I0))') 'IOPEN NSPIN',iOpen,nSpin
   Write (u6,'(A)') 'FindSpinIndex: spin pattern not found'
   Call Abend()
End If

If (nMatch > 1_iwp) Then
   Write(u6,'(A,1X,I0)') 'FindSpinIndex: duplicate count',nMatch
   Write (u6,'(A)') 'FindSpinIndex: duplicate spin pattern'
   Call Abend()
End If

End Subroutine FindSpinIndex

Subroutine FindOCIndex(iDOC,iSOC,OCDOC,OCSOC,nOC,iOC)

Implicit None

Integer(kind=iwp), Intent(In)  :: iDOC
Integer(kind=iwp), Intent(In)  :: iSOC
Integer(kind=iwp), Intent(In)  :: nOC
Integer(kind=iwp), Intent(In)  :: OCDOC(:)
Integer(kind=iwp), Intent(In)  :: OCSOC(:)
Integer(kind=iwp), Intent(Out) :: iOC

Integer(kind=iwp) :: i
Integer(kind=iwp) :: nMatch

If (nOC < 1_iwp) Then
   Write(u6,'(A)') 'FindOCIndex: no configurations for symmetry'
   Call Abend()
End If

If (nOC > Size(OCDOC) .or. nOC > Size(OCSOC)) Then
   Write(u6,'(A)') 'FindOCIndex: invalid configuration count'
   Call Abend()
End If

iOC = 0_iwp
nMatch = 0_iwp

Do i=1,nOC

   If (OCDOC(i) /= iDOC) Cycle
   If (OCSOC(i) /= iSOC) Cycle

   iOC = i
   nMatch = nMatch + 1_iwp

End Do

If (nMatch == 0_iwp) Then
   Write(u6,'(A)') 'FindOCIndex: occupation configuration not found'
   Write(u6,'(A,2(1X,I0))') 'IDOC ISOC',iDOC,iSOC
   Write(u6,'(A,1X,I0)') 'NOC',nOC
   Write(u6,'(A)') 'FindOCIndex: configuration not found'
   Call Abend()
End If

If (nMatch > 1_iwp) Then
   Write(u6,'(A)') 'FindOCIndex: duplicate occupation configuration'
   Write(u6,'(A,2(1X,I0))') 'IDOC ISOC',iDOC,iSOC
   Write(u6,'(A,1X,I0)') 'NMATCH',nMatch
   Write(u6,'(A)') 'FindOCIndex: duplicate configuration'
   Call Abend()
End If

End Subroutine FindOCIndex

Subroutine FindSOCIndex(iDetOC,nElec,SOCTable,iSOC)

Implicit None

Integer(kind=iwp), Intent(In)  :: nElec
Integer(kind=iwp), Intent(In)  :: iDetOC(nElec)
Integer(kind=iwp), Intent(In)  :: SOCTable(:,:)
Integer(kind=iwp), Intent(Out) :: iSOC

Integer(kind=iwp) :: i
Integer(kind=iwp) :: j
Integer(kind=iwp) :: nSOC
Logical           :: Match

nSOC = Size(SOCTable,1)
iSOC = 0_iwp

Do j=1,Size(SOCTable,2)

   Match = .True.

   Do i=1,nSOC
      If (Count(iDetOC == SOCTable(i,j)) /= 1) Then
         Match = .False.
         Exit
      End If
   End Do

   If (Match) Then
      iSOC = j
      Exit
   End If

End Do

If (iSOC == 0_iwp) Then
   Error Stop 'FindSOCIndex: SOC list not found'
End If

End Subroutine FindSOCIndex

Subroutine FindDOCIndex(iDetOC,nElec,DOCTable,iDOC)

Use Definitions, Only: iwp

Implicit None

Integer(kind=iwp), Intent(In)  :: nElec
Integer(kind=iwp), Intent(In)  :: iDetOC(nElec)
Integer(kind=iwp), Intent(In)  :: DOCTable(:,:)
Integer(kind=iwp), Intent(Out) :: iDOC

Integer(kind=iwp) :: DOCOcc(Size(DOCTable,1))
Integer(kind=iwp) :: nDOC
Integer(kind=iwp) :: nDOCTable
Integer(kind=iwp) :: iElec
Integer(kind=iwp) :: jDOC

nDOC      = 0_iwp
nDOCTable = Size(DOCTable,2)
iElec     = 1_iwp
iDOC      = 0_iwp

Do While (iElec <= nElec)

   If (iElec < nElec) Then
      If (iDetOC(iElec) == iDetOC(iElec+1_iwp)) Then

         nDOC = nDOC + 1_iwp

         If (nDOC > Size(DOCOcc)) Then
            Error Stop 'FindDOCIndex: too many double occupations'
         End If

         DOCOcc(nDOC) = iDetOC(iElec)
         iElec = iElec + 2_iwp
         Cycle

      End If
   End If

   iElec = iElec + 1_iwp

End Do

If (nDOC == 0_iwp) Then

   If (Size(DOCTable,1) /= 0) Then
      Error Stop 'FindDOCIndex: expected zero-length DOC table'
   End If

   If (Size(DOCTable,2) /= 1) Then
      Error Stop 'FindDOCIndex: invalid empty DOC table'
   End If

   iDOC = 1_iwp
   Return

End If

If (nDOC /= Size(DOCTable,1)) Then
   Error Stop 'FindDOCIndex: inconsistent number of double occupations'
End If

Do jDOC = 1, nDOCTable
   If (All(DOCOcc(1:nDOC) == DOCTable(:,jDOC))) Then
      iDOC = jDOC
      Exit
   End If
End Do

If (iDOC == 0_iwp) Then
   Error Stop 'FindDOCIndex: DOC occupation not found'
End If

End Subroutine FindDOCIndex

Subroutine ExtractOccupationPattern(iDetOC,iDetSP,nElec, &
                                    DOCOcc,SOCOcc,OpenSP,  &
                                    nDouble,nOpen)

Use Definitions, Only: iwp

Implicit None

Integer(kind=iwp), Intent(In)  :: nElec
Integer(kind=iwp), Intent(In)  :: iDetOC(nElec)
Integer(kind=iwp), Intent(In)  :: iDetSP(nElec)

Integer(kind=iwp), Intent(Out) :: DOCOcc(nElec)
Integer(kind=iwp), Intent(Out) :: SOCOcc(nElec)
Integer(kind=iwp), Intent(Out) :: OpenSP(nElec)
Integer(kind=iwp), Intent(Out) :: nDouble
Integer(kind=iwp), Intent(Out) :: nOpen

Integer(kind=iwp) :: iElec

DOCOcc = 0_iwp
SOCOcc = 0_iwp
OpenSP = 0_iwp

nDouble = 0_iwp
nOpen   = 0_iwp

iElec = 1_iwp

Do While (iElec <= nElec)

   If (iElec < nElec) Then

      If (iDetOC(iElec) == iDetOC(iElec+1_iwp)) Then

         ! A repeated orbital number represents a doubly
         ! occupied orbital. BuildOrderedDeterminant places
         ! the alpha electron before the beta electron.

         If (iDetSP(iElec) /= 1_iwp .or. &
             iDetSP(iElec+1_iwp) /= -1_iwp) Then
            Error Stop &
               'ExtractOccupationPattern: invalid double occupation'
         End If

         nDouble = nDouble + 1_iwp
         DOCOcc(nDouble) = iDetOC(iElec)

         iElec = iElec + 2_iwp
         Cycle

      End If

   End If

   ! A non-repeated orbital number represents a singly
   ! occupied orbital.

   If (Abs(iDetSP(iElec)) /= 1_iwp) Then
      Error Stop 'ExtractOccupationPattern: invalid spin label'
   End If

   nOpen = nOpen + 1_iwp
   SOCOcc(nOpen) = iDetOC(iElec)
   OpenSP(nOpen) = iDetSP(iElec)

   iElec = iElec + 1_iwp

End Do

! Every electron must belong either to a doubly occupied
! orbital or to a singly occupied orbital.

If (2_iwp*nDouble+nOpen /= nElec) Then
   Error Stop 'ExtractOccupationPattern:*inconsistent electron count'
End If

End Subroutine ExtractOccupationpattern

Subroutine BuildOrderedDeterminant(AlphaString,BetaString, &
                                   nAlpha,nBeta,            &
                                   iDetOC,iDetSP,iSign)

Use Definitions, Only: iwp

Implicit None

Integer(kind=iwp), Intent(In)  :: nAlpha
Integer(kind=iwp), Intent(In)  :: nBeta
Integer(kind=iwp), Intent(In)  :: AlphaString(Max(1_iwp,nAlpha))
Integer(kind=iwp), Intent(In)  :: BetaString(Max(1_iwp,nBeta))

Integer(kind=iwp), Intent(Out) :: iDetOC(nAlpha+nBeta)
Integer(kind=iwp), Intent(Out) :: iDetSP(nAlpha+nBeta)
Integer(kind=iwp), Intent(Out) :: iSign

Integer(kind=iwp) :: iAlpha
Integer(kind=iwp) :: iBeta
Integer(kind=iwp) :: iElec
Integer(kind=iwp) :: nCross

iAlpha = 1_iwp
iBeta  = 1_iwp
iSign  = 1_iwp

Do iElec = 1, nAlpha+nBeta

   If (iAlpha <= nAlpha .and. iBeta <= nBeta) Then

      If (AlphaString(iAlpha) <= BetaString(iBeta)) Then

         ! Add the next alpha electron.  The use of <= ensures
         ! that alpha precedes beta in a doubly occupied orbital.

         iDetOC(iElec) = AlphaString(iAlpha)
         iDetSP(iElec) = 1_iwp

         iAlpha = iAlpha + 1_iwp

      Else

         ! Add the next beta electron.

         iDetOC(iElec) = BetaString(iBeta)
         iDetSP(iElec) = -1_iwp

         ! This beta electron passes all remaining alpha electrons.

         nCross = nAlpha - iAlpha + 1_iwp

         If (Mod(nCross,2_iwp) /= 0_iwp) Then
            iSign = -iSign
         End If

         iBeta = iBeta + 1_iwp

      End If

   Else If (iBeta > nBeta) Then

      ! No beta electrons remain.

      iDetOC(iElec) = AlphaString(iAlpha)
      iDetSP(iElec) = 1_iwp

      iAlpha = iAlpha + 1_iwp

   Else

      ! No alpha electrons remain.

      iDetOC(iElec) = BetaString(iBeta)
      iDetSP(iElec) = -1_iwp

      iBeta = iBeta + 1_iwp

   End If

End Do

! Basic consistency checks.

Do iElec = 2, nAlpha+nBeta
   If (iDetOC(iElec) < iDetOC(iElec-1_iwp)) Then
      Error Stop 'BuildOrderedDeterminant: unordered occupation list'
   End If

   If (iDetOC(iElec) == iDetOC(iElec-1_iwp)) Then
      If (iDetSP(iElec-1_iwp) /= 1_iwp .or. &
          iDetSP(iElec)       /= -1_iwp) Then
         Error Stop 'BuildOrderedDeterminant: invalid double occupation'
      End If
   End If
End Do

End Subroutine BuildOrderedDeterminant

  subroutine Count_Strings_By_Symmetry(SGS,nElec,nStringSym)

    implicit none

    type(SGStruct), intent(in) :: SGS

    integer(kind=iwp), intent(in)  :: nElec
    integer(kind=iwp), intent(out) :: nStringSym(SGS%nSym)

    integer(kind=iwp) :: iOrb
    integer(kind=iwp) :: iSym
    integer(kind=iwp) :: nOcc
    integer(kind=iwp) :: newSym

    integer(kind=iwp) :: nString(0:nElec,SGS%nSym)
    integer(kind=iwp) :: nStringNew(0:nElec,SGS%nSym)

    ! nString(nOcc,iSym) is the number of strings constructed
    ! from the orbitals processed so far that contain nOcc
    ! electrons and have spatial symmetry iSym.
    nString(:,:) = 0_iwp

    ! The empty string is totally symmetric.
    nString(0,1) = 1_iwp

    do iOrb = 1, SGS%nLev

       ! Contributions in which orbital iOrb is unoccupied.
       nStringNew(:,:) = nString(:,:)

       ! Contributions in which orbital iOrb is occupied.
       do nOcc = 0, min(nElec-1_iwp,iOrb-1_iwp)
          do iSym = 1, SGS%nSym

             if (nString(nOcc,iSym) == 0_iwp) cycle

             newSym = Mul(iSym,SGS%ISm(iOrb))

             nStringNew(nOcc+1,newSym) =               &
                nStringNew(nOcc+1,newSym) +            &
                nString(nOcc,iSym)

          end do
       end do

       nString(:,:) = nStringNew(:,:)

    end do

    nStringSym(:) = nString(nElec,:)

  end subroutine Count_Strings_By_Symmetry

subroutine GenerateStringSet(nElec,targetSym,nStrings,Strings)

  integer(kind=iwp), intent(in) :: nElec
  integer(kind=iwp), intent(in) :: targetSym
  integer(kind=iwp), intent(in) :: nStrings

  integer(kind=iwp), intent(out) :: &
      Strings(max(1_iwp,nElec),nStrings)

  integer(kind=iwp), allocatable :: String(:)

  integer(kind=iwp) :: iString
  integer(kind=iwp) :: nOrb

  nOrb = SGS%nLev

  Strings = 0_iwp
  iString = 0_iwp

  ! The empty string has totally symmetric spatial symmetry.
  if (nElec == 0_iwp) then

    if (targetSym == 1_iwp) then
      if (nStrings /= 1_iwp) then
        error stop 'GenerateStringSet: invalid empty-string count'
      end if
    else
      if (nStrings /= 0_iwp) then
        error stop 'GenerateStringSet: invalid empty-string symmetry'
      end if
    end if

    return

  end if

  if (nElec < 0_iwp .or. nElec > nOrb) then
    error stop 'GenerateStringSet: invalid electron count'
  end if

  Call mma_allocate(String,nElec,Label='String')
  String = 0_iwp

  ! Generate combinations in the same colexicographic order
  ! used by LUCIA.
  call GenerateLevel(nElec,nOrb,nElec,targetSym,nStrings,              &
                     String,Strings,iString)

  if (iString /= nStrings) then
    error stop 'GenerateStringSet: incorrect number of strings'
  end if

  Call mma_deallocate(String)

End subroutine GenerateStringSet

recursive subroutine GenerateLevel(nChoose,maxOrb,nElec,targetSym,nStrings,           &
                                   String,Strings,iString)

  integer(kind=iwp), intent(in) :: nChoose
  integer(kind=iwp), intent(in) :: maxOrb
  integer(kind=iwp), intent(in) :: nElec
  integer(kind=iwp), intent(in) :: targetSym
  integer(kind=iwp), intent(in) :: nStrings

  integer(kind=iwp), intent(inout) :: String(nElec)
  integer(kind=iwp), intent(out) :: Strings(nElec,nStrings)
  integer(kind=iwp), intent(inout) :: iString

  integer(kind=iwp) :: i
  integer(kind=iwp) :: iOrb
  integer(kind=iwp) :: iSym

  if (nChoose == 0_iwp) then

    iSym = 1_iwp

    do i = 1, nElec
      iSym = Mul(iSym,SGS%ISm(String(i)))
    end do

    if (iSym == targetSym) then

      iString = iString + 1_iwp

      if (iString > nStrings) then
        error stop 'GenerateLevel: too many strings'
      end if

      Strings(:,iString) = String(:)

    end if

    return

  end if

  ! Generate strings in colexicographic order. The highest
  ! occupied orbital is the outermost ordering criterion.
  do iOrb = nChoose, maxOrb

    String(nChoose) = iOrb

    call GenerateLevel(                               &
        nChoose-1_iwp,iOrb-1_iwp,                     &
        nElec,targetSym,nStrings,                     &
        String,Strings,iString)

  end do

end subroutine GenerateLevel

end subroutine Mk_SDREO_SGUGA

pure function SOC_Sym(ISm,nSOC,SOCOcc) result(iSym)
use Symmetry_Info, only: Mul
use Definitions,   only: iwp
implicit none

integer(iwp), intent(in) :: ISm(:)
integer(iwp), intent(in) :: nSOC
integer(iwp), intent(in) :: SOCOcc(:)

integer(iwp) :: iSym
integer(iwp) :: i

iSym = 1
do i=1,nSOC
   iSym = Mul(iSym,ISm(SOCOcc(i)))
end do

end function SOC_Sym


subroutine Count_OpenShell_OCPs(ISm,nLev,nOpen,nSym,nOpenOCP)

use Symmetry_Info, only: Mul
use Definitions,   only: iwp

implicit none

integer(iwp), intent(in)  :: nLev
integer(iwp), intent(in)  :: nOpen
integer(iwp), intent(in)  :: nSym
integer(iwp), intent(in)  :: ISm(nLev)

integer(iwp), intent(out) :: nOpenOCP(nSym)

nOpenOCP(:) = 0

if (nOpen == 0) then
   nOpenOCP(1) = 1
   return
endif

call recurse(1,1,1)

if (sum(nOpenOCP) /= binomial(nOpen,nLev)) then
   write(u6,*) 'Count_OpenShell_OCPs failed'
   write(u6,*) 'NOPEN = ',nOpen
   write(u6,*) 'SUM(NOPENOCP) = ',sum(nOpenOCP)
   write(u6,*) 'EXPECTED      = ',binomial(nOpen,nLev)
   call Abend()
endif

if (any(nOpenOCP < 0)) then
   write(u6,*) 'Count_OpenShell_OCPs failed'
   write(u6,*) 'Negative count encountered'
   call Abend()
endif


contains

   recursive subroutine recurse(iPicked,iStart,iSymProd)

      integer(iwp), intent(in) :: iPicked
      integer(iwp), intent(in) :: iStart
      integer(iwp), intent(in) :: iSymProd

      integer(iwp) :: i

      if (iPicked > nOpen) then
         nOpenOCP(iSymProd) = nOpenOCP(iSymProd) + 1
         return
      endif

      do i=iStart,nLev-(nOpen-iPicked)

         call recurse( &
            iPicked+1, &
            i+1, &
            Mul(iSymProd,ISm(i)) )

      enddo

   end subroutine recurse

end subroutine Count_OpenShell_OCPs


subroutine MKMAW(SGS)
! CONSTRUCT A MODIFIED DIRECT ARC WEIGHT TABLE

type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, ISUM, IU, IV

call mma_allocate(SGS%MAW,[1,SGS%nVert],[0,3],Label='SGS%MAW')
SGS%MAW(:,:) = 0

! COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO MAW:
! From the first vertex of the MIDLEV to the last vertex copy the
! Direct Arc Weight table (DAW)
SGS%MAW(SGS%MVSta:SGS%nVert,0:3) = SGS%DAW(SGS%MVSta:SGS%nVert,0:3)

! COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO MAW. HOWEVER,
!    NOTE THAT THE MAW TABLE IS ACCESSED BY THE UPPER VERTEX.
do IU=1,SGS%MVSta-1          ! Loop over the vertices before the first vertex of level MIDLEV
  do IC=0,3                  ! Loop over step vector
    ID = SGS%Down(IU,IC)     ! Given vertex IU get the vertex index ID from which the step IC originates
                             ! If non-zero put in the reverse direct arc weight table (RAW) values
    if (ID /= 0) SGS%MAW(IU,IC) = SGS%RAW(ID,IC)
  end do
end do
! FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
ISUM = 1
do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    IU = SGS%Up(IV,IC)
    if (IU == 0) cycle
    SGS%MAW(IU,IC) = SGS%MAW(IU,IC)+ISUM
  end do
  ISUM = ISUM+SGS%RAW(IV,4)
end do

do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    if (SGS%Down(IV,IC) == 0) cycle
    SGS%MAW(IV,IC) = SGS%MAW(IV,IC)+ISUM
  end do
  ISUM = ISUM+SGS%DAW(IV,4)
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' THE MODIFIED ARC WEIGHT TABLE IN MKMAW:'
do IV=1,SGS%nVert
  write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%MAW(IV,0:3)
end do
write(u6,*)
#endif

end subroutine MKMAW

subroutine MKSEG(SGS,CIS,EXS)
! PURPOSE: CONSTRUCT THE TABLES ISGM AND VSGM.
! ISGM(IVLT,ISGT) REFERS TO A SEGMENT OF THE SEGMENT TYPE
!    ISGT=1,..,nSeg, WHOSE TOP LEFT VERTEX IS IVLT. ISGM GIVES
!    ZERO IF THE SEGMENT IS IMPOSSIBLE IN THE GRAPH DEFINED BY
!    THE PALDUS TABLE DRT, ELSE IT IS THE BOTTOM LEFT VERTEX
!    NUMBER OF THE SEGMENT. THE SEGMENT VALUE IS THEN VSGM.

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp) :: I_N, IA, IAL, IB, IBL, IBR, INL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB, IVLT, IVRB, IVRT, LEV, MV, MVLL
real(kind=wp) :: V

call mma_allocate(CIS%IVR,SGS%nVert,4,Label='CIS%IVR')
call mma_allocate(CIS%ISGM,SGS%nVert,nSegCat,Label='CIS%ISGM')
call mma_allocate(CIS%VSGM,SGS%nVert,nSegCat,Label='CIS%VSGM')
call mma_allocate(EXS%MVL,CIS%nMidV,2,Label='EXS%MVL')
call mma_allocate(EXS%MVR,CIS%nMidV,2,Label='EXS%MVR')

! CONSTRUCT THE IVR TABLE.
! Make list of right vertices that connects to a left vertex on the same level
!  case 1: delta(b) = -1 , case 2: delta(b) = + 1
CIS%IVR(:,:) = 0   ! initiate
do LEV=1,SGS%nLev    ! Loop over all levels
  IV1 = SGS%LTV(LEV)
  IV2 = SGS%LTV(LEV-1)-1
  do IVL=IV1,IV2-1          ! loop over all vertices of level LEV, right to left, but the last vertex
    IAL = SGS%DRT(IVL,ATAB)! Pick up the a- and b-value of the left vertex
    IBL = SGS%DRT(IVL,BTAB)
    INL = 2*IAL+IBL
    do IV=IVL+1,IV2         ! loop over all vertices of level LEV, right to left, but the first vertex
      IA = SGS%DRT(IV,ATAB)! Pick up the a-value of the right vertex
      IB = SGS%DRT(IV,BTAB)
      I_N = 2*IA+IB

!     The number of particles in the node decrease left to right. If the difference between two vertices
!     is higher than 1 then there is no reason to explore further vertices to the right.

        if (INL-I_N > 1) cycle

!     Test if right vertex is consistent with a tail or an intermediate segment
!     Delta(b) = B(left vertex) - B(right vertex)
!     The intermediate segments are divided into five with delta(b)=-1, and five with delat(b)=+1
!     The tail segments are divided into two with delta(b)=-1, and two with delta(b)=+1

!     2a+b=N, where N is the number of electrons in the CSF described by the vertex
!     b=2S, where S is the total electonic spinf of the CSF described by the vertex

!     Tabulation of valid cases:
!     tail (u0,2d): delta(N)=-1, delta(b)=-1 => IA=IAL, and IB=IBL-1
!     tail (d0,2u): delta(N)=-1, delta(b)=+1 => IA=IAL-1, and IB=IBL-1
!     intermediate segment(00,uu,ud,dd,22): delta(N)=-1, delta(b)=-1
!     intermediate segment(00,uu,du,dd,22): delta(N)=-1, delta(b)=+1

      if (IA == IAL) then
        if (IB == (IBL-1)) CIS%IVR(IVL,1) = IV
      else if (IA == (IAL-1)) then
        if (IB == (IBL+1)) CIS%IVR(IVL,2) = IV
      end if

    end do
  end do
end do

! The partners of the triplet classes 5 and 6 (db = -2, +2): the vertex of
! the same level with the same number of electrons and b(right) = b(left)
! -+ 2, to the left or to the right.
do LEV=1,SGS%nLev
  IV1 = SGS%LTV(LEV)
  IV2 = SGS%LTV(LEV-1)-1
  do IVL=IV1,IV2
    IAL = SGS%DRT(IVL,ATAB)
    IBL = SGS%DRT(IVL,BTAB)
    do IV=IV1,IV2
      IA = SGS%DRT(IV,ATAB)
      IB = SGS%DRT(IV,BTAB)
      if ((IA == IAL+1) .and. (IB == IBL-2)) CIS%IVR(IVL,3) = IV
      if ((IA == IAL-1) .and. (IB == IBL+2)) CIS%IVR(IVL,4) = IV
    end do
  end do
end do

! CONSTRUCT THE MVL AND MVR TABLES:
EXS%MVR(:,:)=0
do IVL=SGS%MVSta,SGS%MVEnd   ! Loop over vertices of the MIDLEV (absolute index)
  MVLL = IVL-SGS%MVSta+1     ! Compute relative index
  ! If there is a valid node for delta(b)=-1, store relative index.
  if (CIS%IVR(IVL,1) /= 0) EXS%MVR(MVLL,1) = CIS%IVR(IVL,1)-SGS%MVSta+1
  if (CIS%IVR(IVL,2) /= 0) EXS%MVR(MVLL,2) = CIS%IVR(IVL,2)-SGS%MVSta+1 ! Dito delta(b)=+1
end do

! constructe the recipical version MVL
EXS%MVL(:,:)=0
do MV=1,CIS%nMidV
  if (EXS%MVR(MV,1) /= 0) EXS%MVL(EXS%MVR(MV,1),1) = MV
  if (EXS%MVR(MV,2) /= 0) EXS%MVL(EXS%MVR(MV,2),2) = MV
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' MIDVERT PAIR TABLES MVL,MVR IN MKSEG:'
write(u6,*) ' MVL TABLE:'
do MV=1,CIS%nMidV
write(u6,1234) MV,EXS%MVL(MV,1),EXS%MVL(MV,2)
enddo
write(u6,*) ' MVR TABLE:'
do MV=1,CIS%nMidV
write(u6,1234) MV,EXS%MVR(MV,1),EXS%MVR(MV,2)
enddo
1234 format('  MV=',I2,'    UPPER WALKS:',2I6)
write(u6,*)
write(u6,*) ' VERTEX PAIR TABLE IVR IN MKSEG:'
do IVL=1,SGS%nVert
write(u6,1235) IVL,CIS%IVR(IVL,1),CIS%IVR(IVL,2)
enddo
1235 format(' IVL=',I2,'    UPPER WALKS:',2I6)
#endif

! A genertic segemet is described as
!
!              (VLT)      (VRT)
!                |          |
!              (VLB)      (VRB)
!
!
!
! INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
CIS%ISGM(:,:) = 0
CIS%VSGM(:,:) = Zero

! FOR EACH VERTEX LOOP OVER ALL POSSIBLE SEGMENTS IT CAN BE A PART OF
do IVLT=1,SGS%nVert     ! Vertex Left Top
  do ISGT=1,nSegCat     ! Segment types (the singlet weights nSeg+1:nSegTot are not tabulated)
    if ((ISGT > nSeg) .and. (ISGT <= nSegTot)) cycle

    ! Find all four vertex indices of this segment, branch out if the segment is invalid
    ITT = ITVPT(ISGT)   ! 0-6, denotes class of the upper nodes of the segment
    IVRT = IVLT               ! Vertex Right Top is the same as the Vertex Top Left, default setting
    if (ITT == 1 .or. ITT == 2) IVRT = CIS%IVR(IVLT,ITT)  ! Pick up the associated upper right vertex, depends on delta(b)
    if (ITT == 5) IVRT = CIS%IVR(IVLT,3)  ! triplet, db = -2
    if (ITT == 6) IVRT = CIS%IVR(IVLT,4)  ! triplet, db = +2
    if (IVRT == 0) cycle          ! Branch out if there is no vertex that will contruct the top of the segment.
    ! Pick up the vertex index for the Vertex Left Bottom as a function of the segment type
    IVLB = SGS%Down(IVLT,IC1(ISGT))
    if (IVLB == 0) cycle          ! Branch out if there is no such vertex
    ! Pick up the vertex index for the Vertex Right Bottom v as a function of the case
    IVRB = SGS%Down(IVRT,IC2(ISGT))
    if (IVRB == 0) cycle          ! Branch out if there is no such vertex
    ! SEGMENT IS NOW ACCEPTED AS POSSIBLE.

    IB = SGS%DRT(IVLT,BTAB)    ! Pick up the b-value from the DRT
    if (ISGT >= iSegTHeadBeg) then
      ! triplet operator: the ASTA values, referred to the b of the right
      ! (ket) top vertex; zero segments (b = 0 with db = 0) are not valid
      IBR = SGS%DRT(IVRT,BTAB)
      if (ISGT <= iSegTHeadEnd) then
        V = ASTA_Head_T(IC1(ISGT),IC2(ISGT),IBR-IB,IBR)
      else if (ISGT <= iSegTTopEnd) then
        V = ASTA_Top_T(IC1(ISGT),IC2(ISGT),IBR-IB,IBR)
      else
        V = ASTA_Weight_T(IC1(ISGT),IC2(ISGT),IBR-IB,IBR)
      end if
      if (V == Zero) cycle
      CIS%ISGM(IVLT,ISGT) = IVLB
      CIS%VSGM(IVLT,ISGT) = V
      cycle
    end if
    CIS%ISGM(IVLT,ISGT) = IVLB  ! Mark that the segment is valid by changing the default value, 0,  to the index of
                                ! lower left vertex being a part of the segment.
!   Note that we use the segment values according to the ASTA method.
    select case (ISVC(ISGT))
      case (0)
        V = Zero
      case (1)
        V = One
      case (2)
        V = -One
      case (3)
        V = (-One)**IB
      case (4)
        V = (-One)**(IB+1)
      case (5)
        V =              sqrt(real(  IB,kind=wp)/real(1+IB,kind=wp))
      case (6)
        V = (-One)**IB * sqrt(real(1+IB,kind=wp)/real(2+IB,kind=wp))
      case (7)
        V =              sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
      case (8)
        V = (-One)**IB * sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
      case (9)
        V = (-One)**IB * One/sqrt(real(  IB,kind=wp)*real(1+IB,kind=wp))
      case (10)
        V = sqrt(real(  IB,kind=wp)*real(2+IB,kind=wp))/real(1+IB,kind=wp)
      case (11)
        V = (-One)**IB * One/sqrt(real(1+IB,kind=wp)*real(2+IB,kind=wp))
      case (12)
        V = sqrt(real(1+IB,kind=wp)*real(3+IB,kind=wp))/real(2+IB,kind=wp)
      case (13)
        V = Two
      case default
        V = Zero ! Dummy assignment
        call Abend()
    end select
    CIS%VSGM(IVLT,ISGT) = V   ! Store away the segment value.
  end do
end do

end subroutine MKSEG

subroutine MkNRCOUP(SGS,CIS,EXS,TRS)

type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
type(TRStruct), intent(in)    :: TRS

integer(kind=iwp) :: IBSYM, INDEO, INDEOB, INDEOT, IP, IPQ, IQ, ISGT, ISYDS1, ISYM, ISYUS1, ITSYM, IVLB, IVLT, LEV, LFTSYM, &
                     MV, MV1, MV2, MVSGM, MV4, MV5, MXDWN, MXUP, NDWNS1, NSGMX, NSGTMP, NT1TMP, NT2TMP, NT3TMP, NT4TMP, NT5TMP, &
                     NUPS1, INDEO0
integer(kind=iwp) :: ICLASS, IT0, NT, K, ITR
integer(kind=iwp) :: ITOP, IBOT
integer(kind=iwp), parameter :: nOpenBands = 4
integer(kind=iwp) :: NRL_OpenBlock
integer(kind=iwp) :: EXS_OpenBlock
logical(kind=iwp) :: ActiveBand(nOpenBands)
integer(kind=iwp) :: band, Memory, INDEO_NRL, INDEO_EXS
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IS, IST, NCP, NUW
#endif
  integer(kind=iwp), allocatable :: NRL(:,:,:)

ActiveBand = .false.
ActiveBand(1) = .true.
ActiveBand(2) = .true.

call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW',safe='*')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW',safe='*')

call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF',safe='*')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF',safe='*')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF',safe='*')

NRL_OpenBlock = nOpenBands * SGS%nLev
EXS_OpenBlock=0
  do band=1,nOpenBands
    if (ActiveBand(band)) EXS_OpenBlock = EXS_OpenBlock+SGS%nLev
  end do

EXS%MxEO = NRL_OpenBlock + (SGS%nLev*(SGS%nLev+1))/2
Memory= EXS_OpenBlock + (SGS%nLev*(SGS%nLev+1))/2

call mma_allocate(EXS%NOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
call mma_allocate(EXS%IOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
EXS%NOCP(:,:,:)=0
EXS%IOCP(:,:,:)=0

! NRL(sym,vertex,indeo): number of valid segment paths of a given symmetry and operators class arriving at a given vertex
call mma_allocate(NRL,[1,SGS%nSym],[1,SGS%nVert],[0,EXS%MxEO],Label='NRL')
! indeo=0 denotes a ordinary half-walk with no open- or closed-loop attached
INDEO0=0

! For upper walks
NRL(:,1:SGS%MVEnd,:) = 0
NRL(1,1,INDEO0) = 1           ! Initiate the head vertex for the INDEO0 to 1

! The NRLs are updated according to graph theory and network analysis.
! Segment originates from a vertex and targets another vertex. In doing so, the
! number of segment walks that access a particular vertex is the sum of all possible walks to
! the vertices of at the preceeding level. In the split graph approach we do this count, for
! the upper walks from the head vertex to the vertices before the mid level vertices, while for the
! lower walks we start at the tail vertex and strive upwards to all the mid vertices.

! IVLT=0 is the root node.
do IVLT = 1, SGS%MVSta-1    ! This formally denotes the upper left node of an segment
  LEV = SGS%DRT(IVLT,LTAB)

  ! loop over all possible transition classes: upper walk, delta(b)=-1/+1, or lower walk
  ! A segment is classified by its top and bottom transition class
  do ICLASS = 0, TRS%nClass-1

    IT0 = TRS%ITR0(IVLT,ICLASS)   ! Offset to list of valid segments
    NT  = TRS%NTR(IVLT,ICLASS)    ! # of valid segments for this vertex

    if (NT == 0) cycle

    do K = 1, NT                  ! loop over all segment that got us here!
      ITR  = IT0 + K              ! counter
      ISGT = TRS%ISGT(ITR)        ! segment index
      IVLB = TRS%IVLB(ITR)        ! vertex index of lower left node of the segments
      ITOP = TRS%ITOP(ITR)        ! index of top connection type (delta(b)=+1 or -1), ITOP=ICLASS
      IBOT = TRS%IBOT(ITR)        ! dito lower connection type
      ISYM = TRS%ISYM(ITR)        ! the symmetry index of the segment

      do ITSYM = 1, SGS%nSym      ! loop over the symmetry index of the top left  vertex
        IBSYM = Mul(ITSYM,ISYM)   ! get the symmetry index of left bottom vertex given the symmetry index of the segment

        select case (TRS%IFLAG(ITR))

        case (TR_WALK)
          ! ordinary upper walk
          NRL(IBSYM,IVLB,INDEO0) = NRL(IBSYM,IVLB,INDEO0) + NRL(ITSYM,IVLT,INDEO0)

        case (TR_OPEN)
          ! loop opening
          IP = LEV
          INDEO = IP + (OpenBand(IBOT)-1)*SGS%nLev
          NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO) + NRL(ITSYM,IVLT,INDEO0)

        case (TR_MID)
          ! intermediate open-loop propagation

          ! Include propagations from all active open loops IP>LEV
          do IP = LEV+1, SGS%nLev
            INDEOT = IP + (OpenBand(ITOP)-1)*SGS%nLev
            INDEOB = IP + (OpenBand(IBOT)-1)*SGS%nLev
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB) + NRL(ITSYM,IVLT,INDEOT)
          end do

        case (TR_CLOSE)
          ! loop closing
          IQ = LEV
          ! Include propagations from all active loops IP>LEV
          do IP = LEV+1, SGS%nLev
            INDEOT = IP + (OpenBand(ITOP)-1)*SGS%nLev
            IPQ = (IP*(IP-1))/2 + IQ
            INDEOB = NRL_OpenBlock + IPQ
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB) + NRL(ITSYM,IVLT,INDEOT)
          end do

        case (TR_TAIL)
          ! tail/downwalk propagation of closed loops
          do IPQ = 1, SGS%nLEV*(SGS%nLev+1)/2
            INDEO = NRL_OpenBlock + IPQ
            NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO) + NRL(ITSYM,IVLT,INDEO)
          end do

        case default
          write(u6,*) 'MkNRCOUP(upper): unexpected TRS%IFLAG = ', TRS%IFLAG(ITR)
          write(u6,*) '  ITR  = ', ITR
          write(u6,*) '  ISGT = ', ISGT
          call Abend()

        end select

      end do
    end do

  end do
end do

! store the accumulated number of upper walks
MXUP = 0
do MV=1,CIS%nMidV                  ! loop over midverticies
  IVLT = MV+SGS%MVSta-1            ! Get the absolute vertex index
  do LFTSYM=1,SGS%nSym             ! Loop over symmetries of the walks that end up at this vertex
    CIS%NOW(1,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO0)  ! Store away the number of walks to this mid vertex.
    MXUP = max(MXUP,CIS%NOW(1,LFTSYM,MV))           ! The max upper walks to any midvertex

    ! ---- open loops ----
    do band = 1, nOpenBands
      if (.not. ActiveBand(band)) cycle

      do IP = 1, SGS%nLev
        INDEO = IP + (band-1)*SGS%nLev
        EXS%NOCP(INDEO,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO)
      end do

    end do

    ! ---- closed loops ----
    do IP = 1, SGS%nLev
    do IQ = 1, IP
      IPQ = IP*(IP-1)/2 + IQ
      INDEO_EXS = IPQ + EXS_OpenBlock
      INDEO_NRL = IPQ + NRL_OpenBlock
      EXS%NOCP(INDEO_EXS,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO_NRL)
    end do
    end do

  end do
end do

! For lower walks
NRL(:,SGS%MVSta:SGS%nVert,:) = 0
NRL(1,SGS%nVert,0) = 1

do IVLT = SGS%nVert-1, SGS%MVSta, -1
  LEV = SGS%DRT(IVLT,LTAB)

  do ICLASS = 0, TRS%nClass-1

    IT0 = TRS%ITR0(IVLT,ICLASS)
    NT  = TRS%NTR(IVLT,ICLASS)

    if (NT == 0) cycle

    do K = 1, NT
      ITR  = IT0 + K
      ISGT = TRS%ISGT(ITR)
      IVLB = TRS%IVLB(ITR)
      ITOP = TRS%ITOP(ITR)
      IBOT = TRS%IBOT(ITR)
      ISYM = TRS%ISYM(ITR)

      do ITSYM = 1, SGS%nSym
        IBSYM = Mul(ITSYM,ISYM)

        select case (TRS%IFLAG(ITR))

        case (TR_TAIL)
          ! ordinary lower walk
          NRL(ITSYM,IVLT,INDEO0) = NRL(ITSYM,IVLT,INDEO0) + NRL(IBSYM,IVLB,INDEO0)

        case (TR_CLOSE)
          ! create open-loop class from below
          IQ = LEV
          INDEO = IQ + (OpenBand(ITOP)-1)*SGS%nLev
          NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO) + NRL(IBSYM,IVLB,INDEO0)

        case (TR_MID)
          ! propagate open-loop class upward
          do IQ = 1, LEV-1
            INDEOT = IQ + (OpenBand(ITOP)-1)*SGS%nLev
            INDEOB = IQ + (OpenBand(IBOT)-1)*SGS%nLev
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT) + NRL(IBSYM,IVLB,INDEOB)
          end do

        case (TR_OPEN)
          ! close open loop into closed-loop class
          IP = LEV
          do IQ = 1, LEV-1
            INDEOB = IQ + (OpenBand(IBOT)-1)*SGS%nLev
            IPQ    = (IP*(IP-1))/2 + IQ
            INDEOT = NRL_OpenBlock + IPQ
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT) + NRL(IBSYM,IVLB,INDEOB)
          end do

        case (TR_WALK)
          ! propagate closed loops upward
          do IPQ = 1, (LEV*(LEV-1))/2
            INDEO = NRL_OpenBlock + IPQ
            NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO) + NRL(IBSYM,IVLB,INDEO)
          end do

        case default
          write(u6,*) 'MkNRCOUP(lower): unexpected TRS%IFLAG = ', TRS%IFLAG(ITR)
          write(u6,*) '  ITR  = ', ITR
          write(u6,*) '  ISGT = ', ISGT
          call Abend()

        end select

      end do
    end do

  end do
end do

MXDWN = 0
do MV=1,CIS%nMidV
  IVLT = MV+SGS%MVSta-1
  do LFTSYM=1,SGS%nSym
    CIS%NOW(2,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO0)
    MXDWN = max(MXDWN,CIS%NOW(2,LFTSYM,MV))

    ! ---- open loops ----
    do band = 1, nOpenBands
      if (.not. ActiveBand(band)) cycle

      do IQ = 1, SGS%nLev
        INDEO = IQ + (band-1)*SGS%nLev
        EXS%NOCP(INDEO,LFTSYM,MV) = max(EXS%NOCP(INDEO,LFTSYM,MV),NRL(LFTSYM,IVLT,INDEO))
      end do

    end do

    ! ---- closed loops ----
    do IP = 1, SGS%nLev
    do IQ = 1, IP
      IPQ = IP*(IP-1)/2 + IQ
      INDEO_EXS = IPQ + EXS_OpenBlock
      INDEO_NRL = IPQ + NRL_OpenBlock
      EXS%NOCP(INDEO_EXS,LFTSYM,MV) = max(EXS%NOCP(INDEO_EXS,LFTSYM,MV),NRL(LFTSYM,IVLT,INDEO_NRL))
    end do
    end do

  end do
end do

EXS%nICOup = 0
do INDEO=1,size(EXS%IOCP,1)
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      EXS%IOCP(INDEO,LFTSYM,MV) = EXS%nICOup
      EXS%nICOup = EXS%nICOup +  EXS%NOCP(INDEO,LFTSYM,MV)
    end do
  end do
end do

call Mk_IOW(CIS,SGS)

call CSFCOUNT(CIS,SGS)

! INSERT FOR USE IN SIGMA ROUTINE

NSGMX  = 1
NDWNS1 = 0  ! Dummy initialization
NSGTMP = MXUP*MXDWN   ! Max size of intermediate sigma block

do MVSGM=1,CIS%nMidV
  MV1 = EXS%MVL(MVSGM,2)
  MV2 = EXS%MVL(MVSGM,1)

  MV4 = EXS%MVR(MVSGM,1)
  MV5 = EXS%MVR(MVSGM,2)

  do ISYUS1=1,SGS%nSym
    NUPS1 =    CIS%NOW(1,ISYUS1,MVSGM)
    do ISYDS1=1,SGS%nSym
      NDWNS1 = CIS%NOW(2,ISYDS1,MVSGM)

      NSGMX = max(NSGMX,CIS%NOCSF(ISYUS1,MVSGM,ISYDS1))

      if (MV1 /= 0) then
        NT4TMP = NUPS1 * CIS%NOW(2,ISYDS1,MV1)
        NT5TMP = CIS%NOW(1,ISYUS1,MV1)* NDWNS1
        NSGTMP = max(NSGTMP,NT4TMP,NT5TMP)
      end if

      if (MV2 /= 0) then
        NT3TMP = NUPS1 * CIS%NOW(2,ISYDS1,MV2)
        NT5TMP = CIS%NOW(1,ISYUS1,MV2) * NDWNS1
        NSGTMP = max(NSGTMP,NT3TMP,NT5TMP)
      end if

      if (MV4 /= 0) then
        NT1TMP = NUPS1 * CIS%NOW(2,ISYDS1,MV4)
        NT5TMP = CIS%NOW(1,ISYUS1,MV4) * NDWNS1
        NSGTMP = max(NSGTMP,NT1TMP,NT5TMP)
      end if

      if (MV5 /= 0) then
        NT2TMP = NUPS1 * CIS%NOW(2,ISYDS1,MV5)
        NT5TMP = CIS%NOW(1,ISYUS1,MV5) * NDWNS1
        NSGTMP = max(NSGTMP,NT2TMP,NT5TMP)
      end if

    end do
  end do
end do

! Set up the code to wheater or not intermediate sigma vectors are to be reused.
call mma_maxDBLE(NT1TMP)
if (NSGTMP*2*SGS%nSym*SGS%MidLev < NT1TMP/4) then
   call mma_allocate(EXS%SGTMP,NSGTMP*2*SGS%nSym*SGS%MidLev,Label='EXS%SGTMP')
   EXS%Reuse_SGTMP=.true.
!  EXS%Reuse_SGTMP=.false.
else
   call mma_allocate(EXS%SGTMP,NSGTMP,Label='EXS%SGTMP')
   EXS%Reuse_SGTMP=.false.
endif
EXS%SGTMP(:)=Zero
EXS%Save_Valid=.false.

#ifdef _DEBUGPRINT_
write(u6,600) MXUP,MXDWN,NSGMX,NSGMX,NSGTMP

NUW=CIS%NUW
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*) ' TOTAL NR OF COUPL COEFFS ',EXS%nICOup
INDEO = EXS_OpenBlock+1
write(u6,*) '         OF TYPE 1&2 ONLY:',EXS%IOCP(INDEO,1,1)
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)

write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN MkNRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,*)
  write(u6,1234) MV,(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,1235) (CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do IST=1,SGS%nSym
    write(u6,1236) IST,(CIS%NOCSF(IS,MV,IST),IS=1,SGS%nSym)
  end do
end do
write(u6,*)
write(u6,*) ' NR OF COUPLING COEFFICIENTS:'
write(u6,*) ' 1. OPEN LOOPS TYPE 1:'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      NCP = EXS%NOCP(IP,IS,MV)
      if (NCP /= 0) write(u6,2345) IP,MV,IS,NCP
    end do
  end do
end do
write(u6,*)
write(u6,*) ' 2. OPEN LOOPS TYPE 2:'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      NCP = EXS%NOCP(SGS%nLev+IP,IS,MV)
      if (NCP /= 0) write(u6,2345) IP,MV,IS,NCP
    end do
  end do
end do
write(u6,*)
write(u6,*) ' 3. CLOSED LOOPS:'
do IP=2,SGS%nLev
  do IQ=1,IP-1
    INDEO = EXS_OpenBlock+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do IS=1,SGS%nSym
        NCP = EXS%NOCP(INDEO,IS,MV)
        if (NCP /= 0) write(u6,2346) IP,IQ,MV,IS,NCP
      end do
    end do
  end do
end do
600 format(/,' MAXIMUM NUMBER OF WALKS', &
           /,' UPPER ',I6,' LOWER ',I6, &
           /,' LENGTH OF LARGEST WORK ARRAYS IN SUBROUTINE SIGMA', &
           /,' TEMPORARY SGM1 ',I7, &
           /,' TEMPORARY SGM2 ',I7, &
           /,' NSGTMP         ',I7)
1234 format('  MV=',I2,'    UPPER WALKS:',2I6)
1235 format('           LOWER WALKS:',8I6)
1236 format(' IST=',I2,'  CONFIGURATIONS:',8I6)
2345 format(' P=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
2346 format(' P=',I2,'  Q=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
#endif

call mma_deallocate(NRL)

end subroutine MkNRCOUP

subroutine MKCOUP(SGS,CIS,EXS,TRS)
! Purpose:
!   Transition-table version of MKCOUP.
!
! Scope of this first version:
!   - current class system only (0,1,2,3)
!   - uses TRS bucket traversal instead of scanning ISGT=1..nSeg
!   - keeps the old path/state logic as intact as possible
!   - keeps old MAW lookups
!   - keeps old VTab logic
!
! Important:
!   ISGPTH(ISEG,LEV) is now a bucket cursor K
!   ISGPTH(IRSEG,LEV) stores the raw segment number ISGT

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(inout) :: CIS
  type(EXStruct), intent(inout) :: EXS
  type(TRStruct), intent(in)    :: TRS

  integer(kind=iwp), allocatable :: VHashKey(:), VHashVal(:)
  integer(kind=iwp) :: hsize
  integer(kind=iwp) :: i, i1, i2
  integer(kind=iwp) :: IAWS, IC, ICL, ICOP, ICR
  integer(kind=iwp) :: IHALF, ILND, INDEO, IP, IPOS, IQ
  integer(kind=iwp) :: ISG, ISGT, ISYM, IT, ITYP, ITYPT
  integer(kind=iwp) :: IVLB, IVLT, IVRT, IVRTOP, IVTAB
  integer(kind=iwp) :: IVTOP, IVTSTA, IVTEND
  integer(kind=iwp) :: K, ITR, IT0, NT
  integer(kind=iwp) :: LEV, LEV1, LEV2, LL, L
  integer(kind=iwp) :: LFTSYM, MV, NCHECK, NVTAB_FINAL
  integer(kind=iwp) :: ITYPMX

  real(kind=wp) :: C

  integer(kind=iwp), allocatable :: ILNDW(:), ISGPTH(:,:)
  real(kind=wp),    allocatable :: val(:), VTab(:)

  ! Rows of ISGPTH
  integer(kind=iwp), parameter :: IVLFT = 1
  integer(kind=iwp), parameter :: ITYPE = 2
  integer(kind=iwp), parameter :: IAWSL = 3
  integer(kind=iwp), parameter :: IAWSR = 4
  integer(kind=iwp), parameter :: ILS   = 5
  integer(kind=iwp), parameter :: ICS   = 6
  integer(kind=iwp), parameter :: ISEG  = 7
  integer(kind=iwp), parameter :: IRSEG = 8

  integer(kind=iwp), parameter :: nOpenBands = 2
  integer(kind=iwp) :: MxEO_Block

  integer(kind=iwp), parameter :: nVTab = 5000
#ifdef _DEBUGPRINT_
  integer(kind=iwp) ICOP1, ICOP2
#endif

  call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
  ! a new coupling table invalidates its index
  call mma_deallocate(EXS%IOrd,safe='*')
  call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')

  MxEO_Block = nOpenBands * SGS%nLev
  ! Special case
  if (SGS%nLev == 1) then
    NVTAB_FINAL = 0
    call mma_allocate(EXS%VTab,NVTAB_FINAL,Label='EXS%VTab')
    return
  end if

  ! NOW is reused as a counter array and restored later by higher-level logic
  CIS%NOW(1:2,1:SGS%nSym,1:CIS%nMidV) = 0

  ! Same idea for NOCP
  EXS%NOCP(:,1:SGS%nSym,1:CIS%nMidV) = 0

  call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
  call mma_allocate(ISGPTH,[1,8],[0,SGS%nLev],Label='ISGPTH')
  call mma_allocate(val,[0,SGS%nLev],Label='val')
  call mma_allocate(VTab,nVTab,Label='VTab')
  hsize = 2*nVTab + 1
  call mma_allocate(VHashKey,[0,hsize-1],Label='VHashKey')
  call mma_allocate(VHashVal,[0,hsize-1],Label='VHashVal')
  VHashKey(:) = -huge(VHashKey)
  VHashVal(:) = 0

  ! Coupling coefficient value table
  NVTAB_FINAL = 2
  VTab(1) =  One
  VTab(2) = -One

  NCHECK = 0

  do IHALF = 1, 2  ! Upper and lower SG factors

    if (IHALF == 1) then
      IVTSTA = 1
      IVTEND = 1
      LEV1   = SGS%nLev
      LEV2   = SGS%MidLev
      ITYPMX = 0
    else
      IVTSTA = SGS%MVSta
      IVTEND = SGS%MVEnd
      LEV1   = SGS%MidLev
      LEV2   = 0
      ITYPMX = 2
    end if

    ! For HALF==1 the level is 1, for HALF==2 the level is SGS%MidLEV
    do IVTOP = IVTSTA, IVTEND ! Loop over all vertex on this level, left to right

      do ITYP = 0, ITYPMX

        IVRTOP = IVTOP   !
        if (ITYP > 0) IVRTOP = CIS%IVR(IVTOP,ITYP)
        if (IVRTOP == 0) cycle

        LEV = LEV1

        ! Initialize path table at starting level
        ISGPTH(IVLFT,LEV) = IVTOP
        ISGPTH(ITYPE,LEV) = ITYP
        ISGPTH(IAWSL,LEV) = 0
        ISGPTH(IAWSR,LEV) = 0
        ISGPTH(ILS,  LEV) = 1
        ISGPTH(ICS,  LEV) = 0
        ISGPTH(ISEG, LEV) = 0   ! bucket cursor
        ISGPTH(IRSEG,LEV) = 0   ! raw segment number

        val(LEV) = One

        ! Walk from top to midlevel, or from midlevel to root
        do while (LEV <= LEV1)

          ITYPT = ISGPTH(ITYPE,LEV)
          IVLT  = ISGPTH(IVLFT,LEV)

          IT0 = TRS%ITR0(IVLT,ITYPT)
          NT  = TRS%NTR(IVLT,ITYPT)

          K = ISGPTH(ISEG,LEV) + 1

          if (K > NT) then
            ! No more transitions left in this bucket: reset and backtrack
            ISGPTH(ISEG,LEV)  = 0
            ISGPTH(IRSEG,LEV) = 0
            LEV = LEV + 1
            cycle
          end if

          ITR  = IT0 + K
          ISGT = TRS%ISGT(ITR)
          IVLB = TRS%IVLB(ITR)

          ICL  = TRS%ICL(ITR)
          ICR  = TRS%ICR(ITR)
          ISYM = TRS%ISYM(ITR)

          ! Right upper vertex
          IVRT = IVLT
          if (TRS%IPRT(ITR) /= 0) IVRT = CIS%IVR(IVLT,TRS%IPRT(ITR))

          ! Store current bucket cursor and raw segment number
          ISGPTH(ISEG,LEV)  = K
          ISGPTH(IRSEG,LEV) = ISGT
          ISGPTH(ICS,LEV)   = ICL

          ! Descend one level
          LEV = LEV - 1

          ! Keep old MAW lookups in this first transition-based version
          ISGPTH(IAWSL,LEV) = ISGPTH(IAWSL,LEV+1) + SGS%MAW(IVLT,ICL)
          ISGPTH(IAWSR,LEV) = ISGPTH(IAWSR,LEV+1) + SGS%MAW(IVRT,ICR)

          val(LEV) = val(LEV+1) * TRS%VSEG(ITR)

          ISGPTH(ILS,LEV)   = Mul(ISYM,ISGPTH(ILS,LEV+1))
          ISGPTH(IVLFT,LEV) = IVLB
          ISGPTH(ITYPE,LEV) = TRS%IBOT(ITR)
          ISGPTH(ISEG,LEV)  = 0
          ISGPTH(IRSEG,LEV) = 0
          ISGPTH(ICS,LEV)   = 0

          if (LEV > LEV2) cycle

          ! ------------------------------------------------------
          ! Bottom of current half-path reached
          ! ------------------------------------------------------
          MV = ISGPTH(IVLFT,SGS%MidLev) + 1 - SGS%MVSta
          LFTSYM = ISGPTH(ILS,LEV2)

          IT = ISGPTH(ITYPE,SGS%MidLev)
          if (IT == 0) IT = 3
          if (ISGPTH(ITYPE,LEV2) == 0) IT = 0

          if (IT == 0) then

            ! Ordinary walk
            ILND = 1 + CIS%NOW(IHALF,LFTSYM,MV)
            IAWS = ISGPTH(IAWSL,LEV2)
            ILNDW(IAWS) = ILND
            CIS%NOW(IHALF,LFTSYM,MV) = ILND

            IPOS = CIS%IOW(IHALF,LFTSYM,MV) + (ILND-1)*CIS%nIpWlk

            do LL = LEV2+1, LEV1, nPack
              IC = 0
              do L = min(LL+nPack-1,LEV1), LL, -1
                IC = 4*IC + ISGPTH(ICS,L)
              end do
              IPOS = IPOS + 1
              CIS%ICase(IPOS) = IC
            end do

          else

            ! Open or closed loop
            IP = 0
            IQ = 0

            do L = LEV2+1, LEV1
              ISG = ISGPTH(IRSEG,L)

              if ((ISG >= 5) .and. (ISG <= 8))  IP = L
              if ((ISG >= 19) .and. (ISG <= 22)) IQ = L
            end do

            if (IP == 0) IP = IQ

            if (IT == 3) THEN
               INDEO = MxEO_Block + (IP*(IP-1))/2 + IQ
            else
               INDEO = SGS%nLev*(IT-1) + IP
            endif

            ! Increment the local index and store
            ICOP = 1 + EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
            ! Add the offset and get the global index
            ICOP = EXS%IOCP(INDEO,LFTSYM,MV) + ICOP

            NCHECK = NCHECK + 1

            if (ICOP > EXS%nICoup) then
              write(u6,*) 'ERROR in MKCOUP: ICOP > EXS%nICoup'
              write(u6,*) ' ICOP      = ', ICOP
              write(u6,*) ' nICoup    = ', EXS%nICoup
              write(u6,*) ' INDEO     = ', INDEO
              write(u6,*) ' MV        = ', MV
              write(u6,*) ' LFTSYM    = ', LFTSYM
              call Abend()
            end if

            ! ==========================================================
            ! M5 BEGIN
            ! Keep all of the following logic identical to the original
            ! MKCOUP until the transition-based path traversal is fully
            ! validated:
            !   - VTab lookup / insertion
            !   - ICoup write
            !   - ILNDW-based Lund renumbering
            !   - final EXS%VTab materialization
            ! ==========================================================

            ! ----------------------------------------------------------
            ! M5: keep VTab logic identical to the original MKCOUP
            ! Do not optimize or refactor this yet.
            ! ----------------------------------------------------------
            C = val(LEV2)

            ! Check the current table to see if this specific value is in the
            ! current table, VTAB
            do i = 1, NVTAB_FINAL
              IVTAB = i
              if (abs(C - VTab(i)) < 1.0e-10_wp) exit
            end do

            ! If not in the table add it at the end
            if (i > NVTAB_FINAL) then
              NVTAB_FINAL = NVTAB_FINAL + 1
              if (NVTAB_FINAL > nVTab) then
                write(u6,*) 'MKCOUP: NVTAB_FINAL exceeded nVTab'
                call Abend()
              end if
              VTab(NVTAB_FINAL) = C
              IVTAB = NVTAB_FINAL
            end if

            ! ----------------------------------------------------------
            ! M5: keep ICoup write logic identical to the original MKCOUP
            ! ----------------------------------------------------------

            EXS%ICoup(1,ICOP) = ISGPTH(IAWSL,LEV2)
            EXS%ICoup(2,ICOP) = ISGPTH(IAWSR,LEV2)
            EXS%ICoup(3,ICOP) = IVTAB

            if (ICOP > EXS%nICoup) then
              write(u6,*) 'MKCOUP: ICOP > EXS%nICoup after write'
              call Abend()
            end if

          end if

          ! Back up one level and continue exploring
          LEV = LEV + 1

        end do
      end do
    end do
  end do

! ------------------------------------------------------------
! M5: keep Lund renumbering identical to the original MKCOUP
! ------------------------------------------------------------
  ! ------------------------------------------------------------
  ! Renumber coupling coefficient indices by Lund scheme
  ! ------------------------------------------------------------
  do ICOP = 1, EXS%nICoup
    i1 = EXS%ICoup(1,ICOP)
    i2 = EXS%ICoup(2,ICOP)
    EXS%ICoup(1,ICOP) = ILNDW(i1)
    EXS%ICoup(2,ICOP) = ILNDW(i2)
  end do

! ============================================================
! M5 END
! ============================================================
  call mma_deallocate(val)
  call mma_deallocate(ISGPTH)
  call mma_deallocate(ILNDW)

#ifdef _DEBUGPRINT_
  ICOP1 = 0
  ICOP2 = 0

  write(u6,*) 'NR OF DIFFERENT COUPLING VALUES: ', NVTAB_FINAL
  ! The value  1.0 is stored at i1=1
  ! The value -1.0 is stored at i1=2
  do ICOP = 1, EXS%nICoup
    i1 = EXS%ICoup(3,ICOP)
    if (i1 == 1) ICOP1 = ICOP1 + 1
    if (i1 == 2) ICOP2 = ICOP2 + 1
  end do

  write(u6,*) 'NR OF COUPS WITH VALUE  1.0: ', ICOP1
  write(u6,*) 'NR OF COUPS WITH VALUE -1.0: ', ICOP2
#endif

  call mma_allocate(EXS%VTab,NVTAB_FINAL,Label='EXS%VTab')
  EXS%VTab(1:NVTAB_FINAL) = VTab(1:NVTAB_FINAL)
  call mma_deallocate(VTab)

  NCP_Max = max(1,maxval(EXS%NOCP(:,:,:)))
  !write(u6,*) 'NCP_MAX=',NCP_MAX
  !write(u6,*) 'SIZE(EXS%NOCP)=',SIZE(EXS%NOCP)
  !write(u6,*) 'EXS%NOCP=',EXS%NOCP
  call mma_allocate(EXS%I1List,NCP_Max,Label='I1List')
  call mma_allocate(EXS%I2List,NCP_Max,Label='I2List')
  call mma_allocate(EXS%XList,NCP_Max,Label='XList')
  call mma_deallocate(VHashKey)
  call mma_deallocate(VHashVal)

end subroutine MKCOUP

subroutine MKSGNUM(STSYM,SGS,CIS,EXS)
! PURPOSE: FOR ALL UPPER AND LOWER WALKS
!          COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
!          REVERSE ARC WEIGHT SUM, RESPECTIVELY.
!          STORE THE DATA IN THE TABLES USGN AND LSGN

integer(kind=iwp), intent(in) :: STSYM
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IC, ICODE, ICONF, IDAWSUM, ILOFF, ILW, IPOS, IRAWSUM, ISTEP, ISYM, IUOFF, IUW, JPOS, JSYM, LEV, LV, MIDV, &
                     NLW, NUW
integer(kind=iwp), allocatable :: ISTEPVEC(:)

EXS%Num_Sym = STSYM
call mma_allocate(EXS%USGN,SGS%MxUp,CIS%nMidV,Label='EXS%USGN',safe='*')
call mma_allocate(EXS%LSGN,SGS%MxDwn,CIS%nMidV,Label='EXS%LSGN',safe='*')
call mma_allocate(ISTEPVEC,SGS%nLev,Label='ISTEPVEC')

! INITIALIZE NUMBERING TABLES

EXS%USGN(:,:) = 0
EXS%LSGN(:,:) = 0

! MAIN LOOP RUNS OVER MIDVERTICES AND SYMMETRIES

ICONF = 0
do MIDV=1,CIS%nMidV
  do ISYM=1,SGS%nSym
    IUOFF = 1+CIS%IOW(1,ISYM,MIDV)
    NUW = CIS%NOW(1,ISYM,MIDV)
    JSYM = Mul(ISYM,STSYM)
    ILOFF = 1+CIS%IOW(2,JSYM,MIDV)
    NLW = CIS%NOW(2,JSYM,MIDV)
    if ((NUW == 0) .or. (NLW == 0)) cycle

    ! LOOP OVER ALL UPPER WALKS

    do IUW=1,NUW
      IPOS = IUOFF+CIS%nIpWlk*(IUW-1)
      ! UNPACK THE UPPER WALK STEP VECTOR
      ICODE = CIS%iCase(IPOS)
      JPOS = 0
      do LEV=SGS%MidLev+1,SGS%nLev
        JPOS = JPOS+1
        if (JPOS == nPack+1) then
          JPOS = 1
          IPOS = IPOS+1
          ICODE = CIS%iCase(IPOS)
        end if
        ISTEP = mod(ICODE,4)
        ISTEPVEC(LEV) = ISTEP
        ICODE = ICODE/4
      end do
      ! GET REVERSE ARC WEIGHT FOR UPPER WALK
      IRAWSUM = 1
      LV = 1
      do LEV=SGS%nLev,SGS%MidLev+1,-1
        IC = ISTEPVEC(LEV)
        LV = SGS%Down(LV,IC)
        IRAWSUM = IRAWSUM+SGS%RAW(LV,IC)
      end do
      EXS%USGN(IRAWSUM,MIDV) = IUW
    end do

    ! LOOP OVER ALL LOWER WALKS

    do ILW=1,NLW
      IPOS = ILOFF+CIS%nIpWlk*(ILW-1)
      ! UNPACK WALK STEP VECTOR
      ICODE = CIS%iCase(IPOS)
      JPOS = 0
      do LEV=1,SGS%MidLev
        JPOS = JPOS+1
        if (JPOS == nPack+1) then
          JPOS = 1
          IPOS = IPOS+1
          ICODE = CIS%iCase(IPOS)
        end if
        ISTEP = mod(ICODE,4)
        ISTEPVEC(LEV) = ISTEP
        ICODE = ICODE/4
      end do
      ! GET DIRECT ARC WEIGHT FOR THE LOWER WALK
      IDAWSUM = 1
      LV = SGS%nVert
      do LEV=1,SGS%MidLev
        IC = ISTEPVEC(LEV)
        LV = SGS%Up(LV,IC)
        IDAWSUM = IDAWSUM+SGS%DAW(LV,IC)
      end do
      EXS%LSGN(IDAWSUM,MIDV) = ICONF
      ICONF = ICONF+NUW
    end do

  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' LSGN IN SUBROUTINE MKSGNUM'
do MIDV=1,CIS%nMidV
  write(u6,'(1X,''MIDV='',I3,/,(20I6))') MIDV,EXS%LSGN(:,MIDV)
end do
write(u6,*)
write(u6,*) ' USGN IN SUBROUTINE MKSGNUM'
do MIDV=1,CIS%nMidV
  write(u6,'(1X,''MIDV='',I3,/,(20I6))') MIDV,EXS%USGN(:,MIDV)
end do
write(u6,*)
#endif

call mma_deallocate(ISTEPVEC)

end subroutine MKSGNUM



pure function NeedsRightPartner(ICLASS)

  logical(kind=iwp) :: NeedsRightPartner
  integer(kind=iwp), intent(in) :: ICLASS

  select case (ICLASS)
    case (1,2)
      NeedsRightPartner = .true.
    case default
  NeedsRightPartner = .false.
  end select

end function NeedsRightPartner

pure function PartnerSlot(ICLASS)

  integer(kind=iwp) :: PartnerSlot
  integer(kind=iwp), intent(in) :: ICLASS

  select case (ICLASS)
  case (1)
    PartnerSlot = 1
  case (2)
    PartnerSlot = 2
    case default
      PartnerSlot = 0
  end select

end function PartnerSlot

pure function OpenBand(ICLASS)

  integer(kind=iwp) :: OpenBand
  integer(kind=iwp), intent(in) :: ICLASS

  select case (ICLASS)
  case (1)
    OpenBand = 1
  case (2)
    OpenBand = 2
  case (3)
    OpenBand = 3
  case (4)
    OpenBand = 4
    case default
      OpenBand = 0
  end select

end function OpenBand

subroutine TRANS_Free(TRS)

  type(TRStruct), intent(inout) :: TRS

  call mma_deallocate(TRS%NTR,    safe='*')
  call mma_deallocate(TRS%ITR0,   safe='*')

  call mma_deallocate(TRS%ISGT,   safe='*')
  call mma_deallocate(TRS%IVLT,   safe='*')
  call mma_deallocate(TRS%IVLB,   safe='*')

  call mma_deallocate(TRS%ITOP,   safe='*')
  call mma_deallocate(TRS%IBOT,   safe='*')

  call mma_deallocate(TRS%ICL,    safe='*')
  call mma_deallocate(TRS%ICR,    safe='*')

  call mma_deallocate(TRS%IPRT,   safe='*')
  call mma_deallocate(TRS%ISYM,   safe='*')
  call mma_deallocate(TRS%IOBAND, safe='*')
  call mma_deallocate(TRS%IFLAG,  safe='*')

  call mma_deallocate(TRS%MAWL,   safe='*')
  call mma_deallocate(TRS%MAWR,   safe='*')

  call mma_deallocate(TRS%VSEG,   safe='*')

end subroutine TRANS_Free

subroutine MKTRANS(SGS,CIS,TRS)

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(TRStruct), intent(inout) :: TRS
  integer(kind=iwp) :: IBOT, ICL, ICLASS, ICR, IFLAG, IOBAND, IPRT, ISGT, ISYM, ITOP, ITR, IVLB, IVLT, LEV, N
  real(kind=wp) :: VSEG
  integer(kind=iwp), allocatable :: IPOS(:,:)

  ! Initialize metadata
  call TRANS_Free(TRS)

  TRS%nClass    = 4
  TRS%nOpenBand = 2

  ! Allocate bucket tables
  call mma_deallocate(TRS%NTR,safe='*')
  call mma_deallocate(TRS%ITR0,safe='*')
  call mma_allocate(TRS%NTR,[1,SGS%nVert],[0,TRS%nClass-1],Label='TRS%NTR')
  call mma_allocate(TRS%ITR0,[1,SGS%nVert],[0,TRS%nClass-1],Label='TRS%ITR0')

  TRS%NTR(:,:)  = 0
  TRS%ITR0(:,:) = 0

  ! First pass: couny transitions per bucket
  ! For every source vertex IVLT, count how many valid segments require each top class ITOP
  do IVLT = 1, SGS%nVert
     do ISGT = 1, nSeg
      IVLB = CIS%ISGM(IVLT,ISGT)
      if (IVLB == 0) cycle
      ITOP = ITVPT(ISGT)
      TRS%NTR(IVLT,ITOP) = TRS%NTR(IVLT,ITOP) + 1
    end do
  end do

  ! Build offsets
  N = 0
  do ICLASS = 0, TRS%nClass-1
    do IVLT = 1, SGS%nVert
      TRS%ITR0(IVLT,ICLASS) = N
      N = N + TRS%NTR(IVLT,ICLASS)
    end do
  end do

  TRS%nTrans = N
  ! Now the bucket layout is fixed.

  ! Allocate the flat transition arrays.

  call mma_allocate(TRS%ISGT,   TRS%nTrans, Label='TRS%ISGT')
  call mma_allocate(TRS%IVLT,   TRS%nTrans, Label='TRS%IVLT')
  call mma_allocate(TRS%IVLB,   TRS%nTrans, Label='TRS%IVLB')

  call mma_allocate(TRS%ITOP,   TRS%nTrans, Label='TRS%ITOP')
  call mma_allocate(TRS%IBOT,   TRS%nTrans, Label='TRS%IBOT')

  call mma_allocate(TRS%ICL,    TRS%nTrans, Label='TRS%ICL')
  call mma_allocate(TRS%ICR,    TRS%nTrans, Label='TRS%ICR')

  call mma_allocate(TRS%IPRT,   TRS%nTrans, Label='TRS%IPRT')
  call mma_allocate(TRS%ISYM,   TRS%nTrans, Label='TRS%ISYM')
  call mma_allocate(TRS%IOBAND, TRS%nTrans, Label='TRS%IOBAND')
  call mma_allocate(TRS%IFLAG,  TRS%nTrans, Label='TRS%IFLAG')

  call mma_allocate(TRS%MAWL,   TRS%nTrans, Label='TRS%MAWL')
  call mma_allocate(TRS%MAWR,   TRS%nTrans, Label='TRS%MAWR')

  call mma_allocate(TRS%VSEG,   TRS%nTrans, Label='TRS%VSEG')

  ! Allocate a fill-position array
  call mma_allocate(IPOS,[1,SGS%nVert],[0,TRS%nClass-1],Label='IPOS')
  IPOS(:,:) = TRS%ITR0(:,:)   ! IPOS is a write pointer per bucket.

  ! Second pass: fill transition arrays

  do IVLT = 1, SGS%nVert      ! Loop over all vertices
    LEV = SGS%DRT(IVLT,LTAB)  ! Pick up the level index of the vertex

    ! loop over all possible segment types allowed for this vertex.
    do ISGT = 1, nSeg

      IVLB = CIS%ISGM(IVLT,ISGT) ! Get the index of Vertex Left Bottom
      if (IVLB == 0) cycle

      ITOP = ITVPT(ISGT)   ! top connection class
      IBOT = IBVPT(ISGT)   ! bottom connection class

      ICL = IC1(ISGT)      ! left case
      ICR = IC2(ISGT)      ! right case

      ! symmetry index of the vertex pair (iVLT-iVLB)
      ISYM = 1
      if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)

      IPRT   = PartnerSlot(ITOP)
      IOBAND = OpenBand(ITOP)
      VSEG   = CIS%VSGM(IVLT,ISGT)

      ! Generate the flag unique for the segment
      select case (ISGT)
      case (1:4)
        IFLAG = TR_WALK
      case (5:8)
        IFLAG = TR_OPEN
      case (9:18)
        IFLAG = TR_MID
      case (19:22)
        IFLAG = TR_CLOSE
      case default
        IFLAG = TR_TAIL
      end select

      IPOS(IVLT,ITOP) = IPOS(IVLT,ITOP) + 1
      ITR = IPOS(IVLT,ITOP)

      TRS%ISGT(ITR) = ISGT
      TRS%IVLT(ITR) = IVLT
      TRS%IVLB(ITR) = IVLB

      TRS%ITOP(ITR) = ITOP
      TRS%IBOT(ITR) = IBOT

      TRS%ICL(ITR)  = ICL
      TRS%ICR(ITR)  = ICR

      TRS%IPRT(ITR)   = IPRT
      TRS%ISYM(ITR)   = ISYM
      TRS%IOBAND(ITR) = IOBAND
      TRS%IFLAG(ITR)  = IFLAG

      TRS%VSEG(ITR) = VSEG

      ! Leave these as placeholders in the first version
      TRS%MAWL(ITR) = 0
      TRS%MAWR(ITR) = 0

    end do
  end do

  call mma_deallocate(IPOS)

end subroutine MKTRANS

!> A new tag for a CI vector (see SG_Epq_Psi, CI_Tag).
function SG_New_Tag() result(iTag)
integer(kind=iwp) :: iTag
SG_Tag_Last = SG_Tag_Last+1
if (SG_Tag_Last >= huge(SG_Tag_Last)) SG_Tag_Last = 1
iTag = SG_Tag_Last
end function SG_New_Tag

subroutine SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM,CI_Tag)

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout), target :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
real(kind=wp), intent(in) :: CPQ, CI(*)
real(kind=wp), intent(_OUT_) :: SGM(*)
! CI_Tag: the tag of the vector CI (SG_New_Tag); the intermediates of a
! previous call with the same tag may be reused. SG_No_Reuse: no reuse.
! Without it the vector is recognized by its fingerprints (CI_Fingerprint,
! a pass over the vector).
integer(kind=iwp), intent(in), optional :: CI_Tag
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2, JQ, JP
real(kind=wp) :: X

! declarations to facilitate the reuse of sigma vectors if possible (the key
! is kept in EXS, next to SGTMP)
logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2, Moved, Tagged, Same_Key
integer(kind=iwp) :: iTag
real(kind=wp) :: CI_ID(2)
integer(kind=iwp) ::  iOff, jOff

!*********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!*********************************************************************

nCSFs=CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
ISYP = SGS%ISm(IP)
ISYQ = SGS%ISm(IQ)
ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)


if (IQ/=IP) then

! Unified code for the case P/=Q
Swap = IQ > IP
Test1 = MERGE(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
Test2 = MERGE(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

if (Test1) then

  ! EXCITING CASE, IQ<IP<=MIDLEV.
  INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
  do MVSGM=1,CIS%nMidV
    do ISYUSG=1,SGS%nSym

      NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
      if (NS1 == 0) cycle
      ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
      ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
      NDWNC = CIS%NOW(2,ISYDC,MVSGM)
      if (NDWNC == 0) cycle
      ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
      NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
      IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
      NCP = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVSGM),Swap)
      if (NCP == 0) cycle
      NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)

      ! get the off-set to the block of compressed coupling coefficients.
      LICP = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVSGM),Swap)
      ! CASE IS: LOWER HALF, EXCITE:
      call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,Swap,Moved)
      ! the walk index (EXS%IOrd) of a reordered block must follow
      if (Moved .and. allocated(EXS%IOrd)) call Index_Coup_Block(EXS,LICP,NCP)
      call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)

    end do
  end do
else if (Test2) then
  ! EXCITING CASE, MIDLEV<IQ<IP
  INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
  do MVSGM=1,CIS%nMidV
    do ISYUSG=1,SGS%nSym
      NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
      if (NS1 == 0) cycle
      ISYUC = Mul(ISYPQ,ISYUSG)
      NUPC = CIS%NOW(1,ISYUC,MVSGM)
      if (NUPC == 0) cycle
      ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
      NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
      ISYDSG = Mul(ISYUSG,ISYSGM)
      NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
      IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
      NCP = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVSGM),Swap)
      if (NCP == 0) cycle
      LICP = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVSGM),Swap)
      ! CASE IS: UPPER HALF, EXCITE:
      call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
    end do
  end do
else
  ! EXCITING CASE, IQ<=MIDLEV<IP
  !
  !   The canonical case, p>q
  !   sigma(Jl(m),Ju(m)) = Sum(Il(m'),Iu(m')) C(Il(m'),Iu(m'))
  !                       x A^p(Ju(m),Iu(m')) x A^q(Ju(m),Iu(m'))
  iOff=1
  JP=MERGE(IQ,IP,SWAP)
  JQ=MERGE(IP,IQ,SWAP)
  ! The intermediates depend on the CI vector, on the upper-half orbital JP
  ! and on the symmetry of the lower-half orbital JQ.
  ! A tagged vector is recognized by its tag (the fingerprints only as the
  ! sanity check of SG_Checks), an untagged one by its fingerprints.
  Reuse_Sigma = .false.
  CI_ID(:) = Zero
  iTag = 0
  if (present(CI_Tag)) iTag = CI_Tag
  Tagged = iTag > 0
  if (EXS%Reuse_SGTMP .and. (iTag >= 0)) then
    Same_Key = EXS%Save_Valid .and. (EXS%Save_Owner == 1) .and. (EXS%Save_Swap .eqv. Swap) .and. (EXS%Save_U == JP) .and. &
               (EXS%Save_LSym == SGS%ISm(JQ)) .and. (EXS%Save_CISym == ISYCI) .and. (EXS%Save_nCSF == nCSFs) .and. &
               (EXS%Save_Tag == iTag)
    if (Tagged) then
      Reuse_Sigma = Same_Key
      if (SG_Checks) then
        CI_ID(:) = CI_Fingerprint(CI,nCSFs)
        if (Reuse_Sigma .and. any(EXS%Save_ID(:) /= CI_ID(:))) then
          write(u6,*) 'SG_Epq_Psi: the same tag for different CI vectors, tag ',iTag
          call Abend()
        end if
      end if
    else
      CI_ID(:) = CI_Fingerprint(CI,nCSFs)
      Reuse_Sigma = Same_Key .and. all(EXS%Save_ID(:) == CI_ID(:))
    end if
  end if
  ! Loop over all mid vertex
  do MVSGM=1,CIS%nMidV  ! m
    ! loop over delta b = -1 and +1
    do MV = 1, 2
      ! Find the connecting vertex, MVX=m' , in the case it exist.
      !
      ! In the the canonical case, p>q, we have the matrix element
      !
      !   A^p(Ju(m),Iu(m')), m'>m
      !
      !   Accessed INDEO,Gamma(Iu),m'
      MVX = MERGE(EXS%MVR(MVSGM,MV),EXS%MVL(MVSGM,MV),SWAP)
      if (MVX == 0) cycle
      if ((Swap .and. (MVX < MVSGM)) .or. ((.not. Swap) .and. (MVX > MVSGM))) then
        write(u6,*) 'SG_Epq_Psi: midvertex partner in the wrong order (MV, MVSGM, MVX, Swap):',MV,MVSGM,MVX,Swap
        call Abend()
      end if

      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        ISYUC = MERGE(Mul(ISYQ,ISYUSG),Mul(ISYP,ISYUSG),Swap)
        ISYDC = MERGE(Mul(ISYP,ISYDSG),Mul(ISYQ,ISYDSG),Swap)

        NUPC = CIS%NOW(1,ISYUC,MVX)
        if (NUPC == 0) cycle
        NDWNC = CIS%NOW(2,ISYDC,MVX)
        if (NDWNC == 0) cycle

!       INDEO =merge(JP, SGS%nLev+JP, MV==1)
        INDEO =(MV-1)*SGS%nLev+JP
        NCP1 = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVX),Swap)
        if (NCP1 == 0) cycle

        ! CASE IS: UPPER HALF, EXCITE:
        LICP1 = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVX),Swap)
        IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

        ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
        ! Put it at the offset jOff
        jOff=iOff
        if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) then
           call Build_Crossing_Temporary( &
                  NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff), &
                  NCP1,EXS%ICOUP(1,LICP1+1),Swap)
        end if

        ! Update the pointer in SGTMP to the next slot.
        if (EXS%Reuse_SGTMP) then
          iOff = iOff + NUPSG*NDWNC
          if (iOff-1 > size(EXS%SGTMP)) then
            write(u6,*) 'iOff-1>SIZE(EXS%SGTMP) -- 1'
            write(u6,*) iOff-1,size(EXS%SGTMP)
            call Abend()
          end if
        end if

!       INDEO = merge(JQ, SGS%nLev+JQ, MV==1)
        INDEO = (MV-1)*SGS%nLev+JQ
        NCP2 = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVX),Swap)
        if (NCP2 == 0) cycle

        ! In case of no reuse compute the sigma vector here.
        if (.not. EXS%Reuse_SGTMP) then
           jOff = 1
           call Build_Crossing_Temporary( &
                  NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff), &
                  NCP1,EXS%ICOUP(1,LICP1+1),Swap)
        end if

        ! CASE IS: LOWER HALF, EXCITE:
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP2 = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVX),Swap)
        call Contract_Crossing_Temporary( &
               CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG, &
               SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP2+1),Swap)

      end do
    end do
  end do

  if (EXS%Reuse_SGTMP) then
    if (.not. Reuse_Sigma) then
      EXS%Save_Valid = .true.
      EXS%Save_Owner = 1
      EXS%Save_Swap = Swap
      EXS%Save_U = JP
      EXS%Save_LSym = SGS%ISm(JQ)
      EXS%Save_CISym = ISYCI
      EXS%Save_nCSF = nCSFs
      EXS%Save_Tag = iTag
      EXS%Save_ID(:) = CI_ID(:)
      ! a vector used only once leaves nothing to reuse
      if (iTag < 0) EXS%Save_Valid = .false.
    end if
  else
    ! SGTMP was overwritten
    EXS%Save_Valid = .false.
  end if

end if

else

  ! THEN THIS IS A SPECIAL CASE.
  if (IP > SGS%MidLev) then

    ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
    ! IP=IQ>MIDLEV
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOUW = CIS%IOW(1,ISYUSG,MVSGM)
        IPSHFT = 2*(IP-1-SGS%MidLev)
        LUW = 1+IOUW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
        IPPOW = 2**IPSHFT
        do I=1,NUPSG
          IC = CIS%ICase(LUW+I*CIS%nIpWlk)
          ICS = mod(IC/IPPOW,4)
          if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          ISTA = ISGSTA+I
          call DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
        end do
      end do
    end do

  else

    ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
    ! IP=IQ < MIDLEV.
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOLW = CIS%IOW(2,ISYDSG,MVSGM)
        IPSHFT = 2*(IP-1)
        LLW = 1+IOLW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
        IPPOW = 2**IPSHFT
        do J=1,NDWNSG
          JC = CIS%ICase(LLW+J*CIS%nIpWlk)
          ICS = mod(JC/IPPOW,4)
          if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          JSTA = ISGSTA + NUPSG*(J-1) + 1
          SGM(JSTA:JSTA+NUPSG-1) = SGM(JSTA:JSTA+NUPSG-1)+X*CI(JSTA:JSTA+NUPSG-1)
        end do
      end do
    end do

  end if
end if

contains

  ! Two weighted checksums of the CI vector, with distinct, non-repeating
  ! weights 1+frac(i*w), so that vectors with equal component sums (e.g.
  ! E_pq e_I = c e_J) are told apart.
  function CI_Fingerprint(C,n) result(F)
    integer(kind=iwp), intent(in) :: n
    real(kind=wp), intent(in) :: C(n)
    real(kind=wp) :: F(2)
    real(kind=wp), parameter :: W1 = 0.6180339887498949_wp, W2 = 0.4142135623730951_wp
    real(kind=wp) :: X
    integer(kind=iwp) :: k
    F(:) = Zero
    do k=1,n
      X = real(k,kind=wp)
      F(1) = F(1)+C(k)*(One+(X*W1-aint(X*W1)))
      F(2) = F(2)+C(k)*(One+(X*W2-aint(X*W2)))
    end do
  end function CI_Fingerprint

  subroutine sort_icoup_block(ICOUP,NCP,swap,Moved)

  integer(kind=iwp), intent(in) :: NCP
  integer(kind=iwp), intent(inout) :: ICOUP(3,NCP)
  logical(kind=iwp), intent(in) :: swap
  logical(kind=iwp), intent(out) :: Moved
    integer(kind=iwp) :: i, j, temp(3)

  Moved = .false.
  do i = 2, NCP
     temp = ICOUP(:,i)
     j = i-1
     do while (j >= 1)
        if (swap) then
           if (ICOUP(1,j) <= temp(1)) exit
        else
           if (ICOUP(2,j) <= temp(2)) exit
        end if
        ICOUP(:,j+1) = ICOUP(:,j)
        j = j - 1
        Moved = .true.
     end do
     ICOUP(:,j+1) = temp
  end do
end subroutine sort_icoup_block

! Reusable construction stage for a crossing excitation operator.
subroutine Build_Crossing_Temporary( &
              NDWN,NUPC,CI,NUPSG,TMP,NCP,ICOUP,Swap)
  integer(kind=iwp), intent(in) :: NDWN,NUPC,NUPSG,NCP
  integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CI(NUPC,NDWN)
  real(kind=wp), intent(out) :: TMP(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: Swap
  TMP=Zero
  call Apply_row(One,NDWN,NUPC,CI,NUPSG,TMP,NCP,ICOUP,swap=Swap)
end subroutine Build_Crossing_Temporary

! Reusable contraction stage for a crossing excitation operator.
subroutine Contract_Crossing_Temporary( &
              Scale,NUP,NDWNC,TMP,NDWNSG,SIGMA,NCP,ICOUP,Swap)
  integer(kind=iwp), intent(in) :: NUP,NDWNC,NDWNSG,NCP
  integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
  real(kind=wp), intent(in) :: Scale,TMP(NUP,NDWNC)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: Swap
  call Apply_col(Scale,NUP,NDWNC,TMP,NDWNSG,SIGMA,NCP,ICOUP,swap=Swap)
end subroutine Contract_Crossing_Temporary

subroutine apply_col(CPQ, NUP, NDWNC, CI, NDWNSG, SIGMA, NCP, ICOUP, swap)

  ! SIGMA(:,I2) = SIGMA(:,I2) + CPQ*V*CI(:,I1) for every coupling (I1,I2,V)
  ! of the lower half (columns of the CSF block). One contiguous column update
  ! per coupling: the blocks are small, so gathering the columns for a GEMM
  ! cost more (copies, temporaries, call overhead) than it gained.
  integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUP,NDWNC)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp) :: I1, I2, ICP, iFrom, iTo
  real(kind=wp) :: X

  ! with swap the roles of the walk indices are exchanged
  iTo = merge(1,2,swap)
  iFrom = merge(2,1,swap)
  do ICP=1,NCP
    I2 = ICOUP(iTo,ICP)
    I1 = ICOUP(iFrom,ICP)
    X = CPQ*EXS%VTab(ICOUP(3,ICP))
    SIGMA(:,I2) = SIGMA(:,I2)+X*CI(:,I1)
  end do

end subroutine apply_col

  subroutine apply_row(CPQ,NDWN,NUPC,CI,NUPSG,SIGMA,NCP,ICOUP,swap)

  ! SIGMA(I2,:) = SIGMA(I2,:) + CPQ*V*CI(I1,:) for every coupling (I1,I2,V)
  ! of the upper half (rows of the CSF block). Columns outside, couplings
  ! inside: every access stays within one column of CI and of SIGMA (no
  ! strided rows, no tile copies).
  integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUPC,NDWN)
  real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp) :: ICP, IDWN
  integer(kind=iwp), pointer :: I1LIST(:), I2LIST(:)
  real(kind=wp), pointer :: XLIST(:)

  XLIST => EXS%XLIST
  I1LIST => EXS%I1LIST
  I2LIST => EXS%I2LIST

  if (swap) then
    do ICP=1,NCP
      I1list(ICP) = ICOUP(2,ICP)
      I2list(ICP) = ICOUP(1,ICP)
      Xlist(ICP) = CPQ*EXS%VTab(ICOUP(3,ICP))
    end do
  else
    do ICP=1,NCP
      I1list(ICP) = ICOUP(1,ICP)
      I2list(ICP) = ICOUP(2,ICP)
      Xlist(ICP) = CPQ*EXS%VTab(ICOUP(3,ICP))
    end do
  end if

  do IDWN=1,NDWN
    do ICP=1,NCP
      SIGMA(I2list(ICP),IDWN) = SIGMA(I2list(ICP),IDWN)+Xlist(ICP)*CI(I1list(ICP),IDWN)
    end do
  end do

  nullify(XLIST,I1LIST,I2LIST)

end subroutine apply_row

end subroutine SG_Epq_Psi

SubRoutine SG_Apply_Ak_Psi(iState,LSM,T,k,Scale,CI,Sigma)

  use Constants, only: Zero, One
  use stdalloc,  only: mma_allocate, mma_deallocate

  Implicit None

  Integer(kind=iwp), Intent(In) :: iState
  Integer(kind=iwp), Intent(In) :: LSM
  Integer(kind=iwp), Intent(In) :: k

  Real(kind=wp), Intent(In) :: T(:,:)
  Real(kind=wp), Intent(In) :: Scale
  Real(kind=wp), Intent(In) :: CI(:)

  Real(kind=wp), Intent(InOut) :: Sigma(:)

  Integer(kind=iwp) :: p
  Integer(kind=iwp) :: nLev

  Real(kind=wp) :: CPK
  Real(kind=wp), Parameter :: Thr = 1.0e-14_wp

  Real(kind=wp), Allocatable :: Tmp(:)
  Logical(kind=iwp) :: Save_Reuse

  Save_reuse=EXS(istate)%ReUse_SGTMP
  EXS(istate)%ReUse_SGTMP=.False.

  nLev = size(T,1)

  Call mma_allocate(Tmp,size(CI),Label='Tmp')

  Do p=1,nLev

!
!   Symmetry constraint
!

     If (SGS(iState)%ISM(p) /= SGS(iState)%ISM(k)) Cycle

!
!   CPK = Tpk - delta_pk
!

     CPK = T(p,k)

     If (p == k) CPK = CPK-One

     CPK = Scale*CPK

     If (Abs(CPK) <= Thr) Cycle

     Tmp(:) = Zero

     Call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState), &
                     p,k,CPK,LSM,CI,Tmp)

     Sigma(:) = Sigma(:) + Tmp(:)

  End Do

  Call mma_deallocate(Tmp)

  EXS(istate)%ReUse_SGTMP=Save_reuse
End SubRoutine SG_Apply_Ak_Psi

#ifdef _NOT_USED_
subroutine SG_Sum_Hpq_Epq_Psi(iState,LSM,Hpq,nLev,CI,Sigma)

  use Constants, only: Zero, One
  use stdalloc,  only: mma_allocate,mma_deallocate

  implicit none

  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp), intent(in) :: LSM
  integer(kind=iwp), intent(in) :: nLev

  real(kind=wp), intent(in)  :: Hpq(nLev,nLev)
  real(kind=wp), intent(in)  :: CI(:)
  real(kind=wp), intent(out) :: Sigma(:)

  integer(kind=iwp) :: ip,iq
  real(kind=wp)  :: Factor
  real(kind=wp), Parameter :: Thr = 1.0d-14

  real(kind=wp), allocatable :: Tmp(:)

  call mma_allocate(Tmp,size(CI),Label='Tmp')

  Sigma(:)=Zero

  do iq=1,nLev
    do ip=1,nLev

      ! Epq may only connect orbitals of the same irrep
      If (SGS(iState)%ISM(ip) /= SGS(iState)%ISM(iq)) Cycle

      Factor = Hpq(ip,iq)

      if(abs(Factor) < Thr) cycle

      call SG_Epq_Psi( SGS(iState), &
                       CIS(iState), &
                       EXS(iState), &
                       ip, iq, One, LSM, &
                       CI, Tmp )

      Sigma(:)=Sigma(:)+Factor*Tmp(:)

    end do
  end do

  call mma_deallocate(Tmp)

end subroutine SG_Sum_Hpq_Epq_Psi
#endif

subroutine SG_EpqEqp_Diag(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,DIAG)

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout), target :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
real(kind=wp), intent(in) :: CPQ
real(kind=wp), intent(inout) :: CI(*)
real(kind=wp), intent(_OUT_) :: DIAG(*)
! todo: cleanup
!integer(kind=iwp) :: INDEO, IOC, ISGSTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
integer(kind=iwp) :: INDEO, IOC, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, LICP, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2

! declarations to facilitate the reuse of the intermediates in SGTMP (the key
! is kept in EXS; here they depend only on the couplings, not on CI)
logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2
integer(kind=iwp) ::  iOff, jOff, JP, JQ

!*********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!*********************************************************************

nCSFs=CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
ISYP = SGS%ISm(IP)
ISYQ = SGS%ISm(IQ)
ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)

Swap = IQ > IP
Test1=MERGE(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
Test2=MERGE(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

if (IQ /= IP) Then
  Swap = IQ > IP
  Test1=MERGE(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
  Test2=MERGE(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)
  If (Test1) Then
    ! EXCITING CASE, IQ<IP<=MIDLEV.
    INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym

        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
        ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)
        if (NDWNC == 0) cycle
!       ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
        NCP = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVSGM),SWAP)
        if (NCP == 0) cycle
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)

        ! off-set to the block of compressed coupling coefficients.
        LICP = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVSGM),SWAP)
        ! CASE IS: LOWER HALF, EXCITE:

!       call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,Swap)
!       call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
        call Apply_colX(CPQ,NUPSG,NDWNC,DIAG(IOC+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)

      end do
    end do
  Else If (Test2) Then
    ! EXCITING CASE, MIDLEV<IQ<IP
    INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,SWAP)
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISYUC = Mul(ISYPQ,ISYUSG)
        NUPC = CIS%NOW(1,ISYUC,MVSGM)
        if (NUPC == 0) cycle
!       ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
        NCP = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVSGM),SWAP)
        if (NCP == 0) cycle
        LICP = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVSGM),SWAP)
        ! CASE IS: UPPER HALF, EXCITE:
!       call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
        call apply_rowX(CPQ,NDWNSG,NUPC,DIAG(IOC+1),NCP,EXS%ICOUP(1,LICP+1),swap=Swap)
      end do
    end do
  Else

    ! EXCITING CASE, IQ<=MIDLEV<IP

    iOff=1
    JP=MERGE(IQ,IP,SWAP)
    JQ=MERGE(IP,IQ,SWAP)
    Reuse_Sigma = EXS%Reuse_SGTMP .and. EXS%Save_Valid .and. (EXS%Save_Owner == 2) .and. (EXS%Save_Swap .eqv. Swap) .and. &
                  (EXS%Save_U == JP) .and. (EXS%Save_LSym == SGS%ISm(JQ)) .and. (EXS%Save_CISym == ISYCI)

    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
        MVX = MERGE(EXS%MVR(MVSGM,MV),EXS%MVL(MVSGM,MV),SWAP)
        if (MVX == 0) cycle
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
          if (NS1 == 0) cycle
!         ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          ISYUC = MERGE(Mul(ISYQ,ISYUSG),Mul(ISYP,ISYUSG),SWAP)
          ISYDC = MERGE(Mul(ISYP,ISYDSG),Mul(ISYQ,ISYDSG),SWAP)

          NUPC = CIS%NOW(1,ISYUC,MVX)
          if (NUPC == 0) cycle
          NDWNC = CIS%NOW(2,ISYDC,MVX)
          if (NDWNC == 0) cycle

          INDEO = MERGE(merge(IQ, SGS%nLev+IQ, MV==1),merge(IP, SGS%nLev+IP, MV==1),SWAP)
          NCP1 = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVX),SWAP)
          if (NCP1 == 0) cycle

          ! CASE IS: UPPER HALF, EXCITE:
          LICP1 = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVX),Swap)
          IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

          ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
          if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) then
             NTMP = NUPC
             EXS%SGTMP(iOff:iOff+NTMP-1) = Zero
!            call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
             call Apply_rowY(One,NUPC,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
          end if

          jOff=iOff
          if (EXS%Reuse_SGTMP) then
            iOff = iOff + NUPC
            if (iOff-1 > size(EXS%SGTMP)) then
              write(u6,*) 'iOff-1>SIZE(EXS%SGTMP) -- 1'
              write(u6,*) iOff-1,size(EXS%SGTMP)
              call Abend()
            end if
          end if

          INDEO = MERGE(merge(IP, SGS%nLev+IP, MV==1),merge(IQ, SGS%nLev+IQ, MV==1),SWAP)
          NCP2 = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVX),SWAP)
          if (NCP2 == 0) cycle

          if (.not. EXS%Reuse_SGTMP) then
             NTMP = NUPC
             EXS%SGTMP(jOff:jOff+NTMP-1) = Zero
!            call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
             call Apply_rowY(One,NUPC,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
          end if

          ! CASE IS: LOWER HALF, EXCITE:
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          LICP2 = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVX),SWAP)
!         call Apply_col(CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP2+1),swap=Swap)
          call Apply_colY(CPQ,NUPC,EXS%SGTMP(jOff),NDWNC,DIAG(IOC+1),NCP2,EXS%ICOUP(1,LICP2+1),swap=Swap)

        end do
      end do
    end do
    if (EXS%Reuse_SGTMP) then
      if (.not. Reuse_Sigma) then
        EXS%Save_Valid = .true.
        EXS%Save_Owner = 2
        EXS%Save_Swap = Swap
        EXS%Save_U = JP
        EXS%Save_LSym = SGS%ISm(JQ)
        EXS%Save_CISym = ISYCI
        EXS%Save_nCSF = nCSFs
      end if
    else
      ! SGTMP was overwritten
      EXS%Save_Valid = .false.
    end if
  End If
Else
  CI(1:nCSFs)=One
  Call SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,DIAG)
  DIAG(1:nCSFs) = DIAG(1:nCSFs)**2
End If

contains

subroutine apply_colX(CPQ, NUP, NDWNSG, SIGMA, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NUP, NDWNSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC1), =.true. (DEX1)

  integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1, idx
  real(kind=wp) :: X
  real(kind=wp), pointer :: VTAB(:)

vTab => EXS%VTab
do ICP=1,NCP
   JLFT = ICOUP(1,ICP)
   JRGT = ICOUP(2,ICP)

   I1 = MERGE(JRGT,JLFT,Swap)

   idx=ICOUP(3,ICP)
   X = (CPQ * VTab(idx))**2
   !$OMP SIMD
   do IUP=1,NUP
     SIGMA(IUP,I1) = SIGMA(IUP,I1) + X
   end do
 end do
vTab => Null()
end subroutine apply_colX

subroutine apply_rowX(CPQ, NDWN, NUPSG, SIGMA, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NDWN, NUPSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ
  real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC2), =.true. (DEX2)

  integer(kind=iwp) :: ICP, ILFT, IRGT, IDWN, I1, idx
  real(kind=wp) :: X
  real(kind=wp), pointer :: VTAB(:)

vTab => EXS%VTab
do ICP=1,NCP
  ILFT = ICOUP(1,ICP)
  IRGT = ICOUP(2,ICP)

  I1=MERGE(IRGT,ILFT,Swap)

  idx=ICOUP(3,ICP)
  X = (CPQ * VTab(idx))**2
  !$OMP SIMD
  do IDWN=1,NDWN
    SIGMA(I1,IDWN) = SIGMA(I1,IDWN) + X
  end do
end do
vTab => Null()
end subroutine apply_rowX

subroutine apply_rowY(CPQ, NUPSG, UP, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NUPSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ
  real(kind=wp), intent(inout) :: UP(NUPSG)
  logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC2), =.true. (DEX2)

! todo: cleanup
! integer(kind=iwp) :: ICP, ILFT, IRGT, I1, I2, idx
  integer(kind=iwp) :: ICP, ILFT, IRGT, I1,     idx
  real(kind=wp) :: X
  real(kind=wp), pointer :: VTAB(:)

vTab => EXS%VTab
do ICP=1,NCP
  ILFT = ICOUP(1,ICP)
  IRGT = ICOUP(2,ICP)

  I1=MERGE(IRGT,ILFT,SWAP)
! I2=MERGE(ILFT,IRGT,SWAP)

  idx=ICOUP(3,ICP)
  X = (CPQ * VTab(idx))**2
  UP(I1) = UP(I1) + X
end do
vTab => Null()
end subroutine apply_rowY

subroutine apply_colY(CPQ, NUP, Up, NDWNSG, SIGMA, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NUP, NDWNSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, Up(NUP)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap   ! =.false. (EXC1), =.true. (DEX1)

! todo: cleanup
! integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1, I2, idx
  integer(kind=iwp) :: ICP, JLFT, JRGT, IUP, I1,     idx
  real(kind=wp) :: X
  real(kind=wp), pointer :: VTAB(:)

vTab => EXS%VTab
do ICP=1,NCP
  JLFT = ICOUP(1,ICP)
  JRGT = ICOUP(2,ICP)

  I1=MERGE(JRGT,JLFT,Swap)
! I2=MERGE(JLFT,JRGT,Swap)

  idx=ICOUP(3,ICP)
  X = (CPQ * VTab(idx))**2
  !$OMP SIMD
  do IUP=1,NUP
    SIGMA(IUP,I1) = SIGMA(IUP,I1) + X*Up(IUP)
  end do
end do
vTab => Null()
end subroutine apply_colY

end subroutine SG_EpqEqp_Diag

#ifdef _NOT_USED_
Subroutine SG_Multiply_By_Nk(iState,LSM,k,tkk,CI)

  use Constants, only: Zero, One, Half, Three, Four
  use stdalloc,  only: mma_allocate,mma_deallocate

  Integer(kind=iwp), Intent(In)    :: iState
  Integer(kind=iwp), Intent(In)    :: LSM
  Integer(kind=iwp), Intent(In)    :: k
  Real(kind=wp) , Intent(In)    :: tkk
  Real(kind=wp) , Intent(InOut) :: CI(:)

  Real(kind=wp) :: b,c

  Real(kind=wp), Allocatable :: N1(:),N2(:)
  Logical(kind=iwp) :: ReUse_Save

  ReUse_Save=EXS(istate)%ReUse_SGTMP
  EXS(istate)%ReUse_SGTMP=.False.

  Call mma_allocate(N1,size(CI),Label='N1')
  Call mma_allocate(N2,size(CI),Label='N2')

  b = Half*(-tkk*tkk + Four*tkk - Three)
  c = Half*(tkk-One)*(tkk-One)

  ! N1 = Ekk Psi

  Call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState),k,k,One,LSM,CI,N1)

  ! N2 = Ekk^2 Psi

  Call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState),k,k,One,LSM,N1,N2)

  CI(:)=CI(:)+b*N1(:)+c*N2(:)

  Call mma_deallocate(N1)
  Call mma_deallocate(N2)

  EXS(istate)%ReUse_SGTMP=ReUse_Save

End Subroutine SG_Multiply_By_Nk
#endif

SubRoutine SG_TRACI(iState,LSM,T,nT,CI,nCSF)

  use molcas, only: MxAct
  use Constants, only: Zero, Half, Three
  use stdalloc,  only: mma_allocate, mma_deallocate

  Implicit None

  Integer(kind=iwp), Intent(In)    :: iState
  Integer(kind=iwp), intent(In) :: LSM
  Integer(kind=iwp), Intent(In)    :: nT
  Integer(kind=iwp), Intent(In)    :: nCSF

  Real(kind=wp), Intent(In), target   :: T(nT)
  Real(kind=wp), Intent(InOut):: CI(nCSF)

  Integer(kind=iwp) :: k, iOff, iSym, nAct, p, kLocal, pLocal
  Integer(kind=iwp) :: nLev, Level(MxAct)

  Real(kind=wp) :: tkk, X

  real(kind=wp), allocatable :: Tmp(:)
  Real(kind=wp), Allocatable :: T2(:,:)


  nLev = SGS(iState)%NLEV

! Call RecPrt('SG_TRACI; CI',' ',CI,1,nCSF)
  Call mma_allocate(T2,nLev,nLev,Label='T2')
! Call RecPrt('T',' ',T,1,nT)

  T2(:,:) = Zero
  iOff = 0

  Do iSym = 1,SGS(iState)%nSym

     nAct = 0

     Do p = 1,nLev
        If (SGS(iState)%ISM(p) == iSym) Then
           nAct = nAct + 1
           Level(nAct) = p
        End If
     End Do

     Do kLocal = 1,nAct
        k = Level(kLocal)

        Do pLocal = 1,nAct
           p = Level(pLocal)
           T2(p,k) = T(iOff + pLocal + (kLocal-1)*nAct)
        End Do
     End Do
     iOff = iOff + nAct*nAct
  End Do
! Call RecPrt('T2',' ',T2,nLev,nLev)
!
! Work array
!
  call mma_allocate(Tmp,nCSF,Label='tmp')
!
! Apply row operators sequentially.
!

  Do k=1,nLev

!    Write (6,*) 'ik=',k
     TMP(:)=Zero

     ! tmp = 1/2 A_k CI
     Call SG_Apply_Ak_Psi(iState,LSM,T2,k,Half,CI,TMP)
!    Call RecPrt('TMP',' ',TMP,1,nCSF)

     tkk = T2(k,k)
     X = Three-tkk
!    Write (6,*) 'Ckk, X=',TKK, X
     CI(:)=CI(:)+X*TMP(:)
!    Call RecPrt('CI',' ',CI,1,nCSF)

     Call SG_Apply_Ak_Psi(iState,LSM,T2,k,One,TMP,CI)
!    Call RecPrt('CI',' ',CI,1,nCSF)

  End Do

!
! Cleanup
!
  call mma_deallocate(Tmp)
  call mma_deallocate(T2)

End SubRoutine SG_TRACI

#ifdef _NOT_USED_
subroutine build_trow_operator(T,k,tkk,H1)

  use constants, only: Zero
  implicit none

  integer(kind=iwp),intent(in) :: k
  real(kind=wp) ,intent(in) :: T(:,:)
  real(kind=wp) ,intent(out):: tkk
  real(kind=wp) ,intent(out):: H1(:,:)

  integer(kind=iwp) :: p, nLev

  nLev = size(T,1)
  H1(:,:)=Zero

  tkk = T(k,k)

  do p=1,nLev
     If (p == k) Cycle
     H1(p,k) = T(p,k)/tkk
  end do

end subroutine build_trow_operator
#endif



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

#ifdef _NOT_USED_
subroutine Encode_CSF( &
           CIS,ISYCI, &
           MV,ISYUP,IUP,IDWN, &
           iCSF )

use Definitions, only: iwp

type(CIStruct), intent(in) :: CIS

integer(kind=iwp), intent(in)  :: ISYCI
integer(kind=iwp), intent(in)  :: MV
integer(kind=iwp), intent(in)  :: ISYUP
integer(kind=iwp), intent(in)  :: IUP
integer(kind=iwp), intent(in)  :: IDWN

integer(kind=iwp), intent(out) :: iCSF

integer(kind=iwp) :: NUP
integer(kind=iwp) :: IOFF

NUP  = CIS%NOW(1,ISYUP,MV)
IOFF = CIS%IOCSF(ISYUP,MV,ISYCI)

iCSF = IOFF + IUP + (IDWN-1)*NUP

end subroutine Encode_CSF
#endif

!> Decode_CSF for CIS(ME_iState), remembering the last two CSFs (least
!> recently used slot replaced): the matrix-element routines decode the same
!> bra and ket for every operator string of a matrix element. A new CSF is
!> located by a binary search over the sorted block offsets ME_Blk.
subroutine Decode_ME(ISYCI,iCSF,MV,ISYUP,IUP,IDWN)

use Definitions, only: iwp

integer(kind=iwp), intent(in)  :: ISYCI, iCSF
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
           MV,ISYUP,IUP,IDWN )

use Definitions, only: iwp,u6

type(CIStruct), intent(in) :: CIS

integer(kind=iwp), intent(in)  :: ISYCI
integer(kind=iwp), intent(in)  :: iCSF

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
      NCI  = CIS%NOCSF(ISYUP,MV,ISYCI)

      ! CSF iCSF lies in this (MV,ISYUP) block
      if ((NCI > 0) .and. (iCSF > IOFF) .and. (iCSF <= IOFF+NCI)) then

         ILOC = iCSF-IOFF

         NUP = CIS%NOW(1,ISYUP,MV)

         IDWN = (ILOC-1)/NUP + 1
         IUP  = ILOC - (IDWN-1)*NUP

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
Return
Do iSym = 1, SGS%nSym
Do p = 1, SGS%nLev
Do q = 1, SGS%nLev

End Do
End Do
End Do

end subroutine MK_FLAT_Epq

subroutine sg_csfdet( nopen,     &
                         ndet,      &
                         detspin,      &
                         ncsf,      &
                         twos,      &
                         cdc )

use Constants, only: Half, One, Two, Zero
use Definitions, only: wp, iwp

implicit none

integer(iwp), intent(in)  :: nopen
integer(iwp), intent(in)  :: ndet
integer(iwp), intent(in)  :: ncsf
integer(iwp), intent(in)  :: detspin(nopen,ndet)
integer(iwp), intent(in)  :: twos(0:nopen,ncsf)

real(wp), intent(out) :: cdc(ndet,ncsf)

integer(iwp) :: idet
integer(iwp) :: icsf
integer(iwp) :: i

real(wp) :: coef, norm
real(wp) :: sign
real(wp) :: s
real(wp) :: m(Max(1,nopen))

logical :: up

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) '========================================'
write(u6,*) ' SG_CSFDET'
write(u6,*) '========================================'
write(u6,*) ' NOPEN, nDet, nCSF = ',nopen,ndet,ncsf
write(u6,*)
#endif

do idet=1,ndet
   call det_mchain(nopen,detspin(:,idet),m)

#ifdef _DEBUGPRINT_
   write(u6,*)
   write(u6,*) '----------------------------------------'
   write(u6,*) ' Determinant ',idet
   write(u6,*) '----------------------------------------'
   write(u6,'(A,100I3)') ' DETSPIN = ', (detspin(i,idet),i=1,nopen)
   write(u6,'(A,100F8.3)') ' MCHAIN  = ', (m(i),i=1,nopen)
#endif

   do icsf=1,ncsf

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) ' CSF = ',icsf
      write(u6,'(A,100I3)') ' TWOS = ', (twos(i,icsf),i=0,nopen)
      write(u6,'(A)',advance='no')
      write(u6,*) ' SCHAIN:'
      do i=1,nopen
         write(u6,'(I4,F10.4)') i,0.5_wp*real(twos(i,icsf),wp)
      end do
#endif

      coef = One
      sign = One

      do i=1,nopen
         s = Half*real(twos(i,icsf),wp)
         ! walk step:
         ! twos(i)=twos(i-1)+1  => up
         ! twos(i)=twos(i-1)-1  => down

         up = (twos(i,icsf) > twos(i-1,icsf))
         if (up) then
            if (detspin(i,idet) == 1) then
               coef = coef*(s+m(i)) / (Two*s)
            else
               coef = coef*(s-m(i)) / (Two*s)
            end if

         else

            if (detspin(i,idet) == 1) then
               coef = coef*(s-m(i)+One) / (Two*s+Two)
               sign = -sign
            else
               coef = coef*(s+m(i)+One) / (Two*s+Two)
            end if

         end if

#ifdef _DEBUGPRINT_
         write(u6,*)
         write(u6,*) ' Step i = ',i
         write(u6,*) ' twos(i-1) = ',twos(i-1,icsf)
         write(u6,*) ' twos(i)   = ',twos(i,icsf)
         write(u6,*) ' S = ',s
         write(u6,*) ' M = ',m(i)

         if (up) then
            write(u6,*) ' walk step = UP'
         else
            write(u6,*) ' walk step = DOWN'
         endif

         if (detspin(i,idet) == 1) then
            write(u6,*) ' determinant spin = ALPHA'
         else
            write(u6,*) ' determinant spin = BETA'
         endif
#endif

      end do

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) ' FINAL PRODUCT = ',coef
      write(u6,*) ' FINAL SIGN    = ',sign
#endif

      if (coef < Zero) then
         cdc(idet,icsf)=Zero
      else
         cdc(idet,icsf)=sign*sqrt(coef)
      end if

#ifdef _DEBUGPRINT_
      write(u6,*) ' CDC(',idet,',',icsf,') = ',cdc(idet,icsf)
#endif

   end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('Final CDC matrix',' ',cdc,ndet,ncsf)
#endif

do icsf=1,ncsf
   norm=Zero

   do idet=1,ndet
      norm=norm+cdc(idet,icsf)**2
   end do

   if (Abs(One-Norm)>1.0e-12_wp) Then
      write(u6,*)
      write(u6,*) ' CSF ',icsf,' norm = ',norm
      Call Abend()
   End If

end do

contains
subroutine det_mchain(nopen,det,m)

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none

integer(iwp), intent(in)  :: nopen
integer(iwp), intent(in)  :: det(Max(1,nopen))
real(wp),     intent(out) :: m(Max(1,nopen))

integer(iwp) :: i

m(1)=real(det(1),wp)-Half

do i=2,nopen
   m(i)=m(i-1)+real(det(i),wp)-Half
end do

end subroutine det_mchain

end subroutine sg_csfdet

#ifdef _TEST_CSFDET_
subroutine test_sg_csfdet()

use Definitions, only: wp, iwp

implicit none

integer(iwp), parameter :: maxopen=8
integer(iwp) :: nopen
integer(iwp) :: mult
integer(iwp) :: ms2

integer(iwp), allocatable :: detspin(:,:)
integer(iwp), allocatable :: twos(:,:)

real(wp), allocatable :: cdc(:,:)
integer(iwp) :: ndet, ncsf

#ifdef _DEBUGPRINT_
integer(iwp) :: idet, icsf

character(len=1) :: walk(Maxopen)
character(len=1) :: str(Maxopen)
integer(iwp) :: i

write(u6,*)
write(u6,*) '===================================='
write(u6,*) ' TEST SG_CSFDET'
write(u6,*) '===================================='
write(u6,*)
#endif

do nopen = 1,maxopen

#ifdef _DEBUGPRINT_
   write(u6,*)
   write(u6,*) '========================================='
   write(u6,*) 'NOPEN = ',nopen
   write(u6,*) '========================================='
#endif

   do mult = mod(nopen,2)+1, nopen+1, 2

      ms2 = mult-1

#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'MULT = ',mult
      write(u6,*) 'MS2  = ',ms2
#endif

      call generate_detspin(nopen,ms2,detspin,ndet)
#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'Determinant strings'
      do idet=1,ndet
         do i=1,nopen
            if (detspin(i,idet) == 1) then
               str(i)='a'
            else
               str(i)='b'
            end if
         end do
         write(u6,'(I4,2X,100A1)') idet,(str(i),i=1,nopen)
      end do
#endif

      call generate_twos(nopen,mult,twos,ncsf)
#ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'CSF walks'
      do icsf=1,ncsf
         do i=1,nopen
            if (twos(i,icsf) > twos(i-1,icsf)) then
               walk(i)='u'
            else
               walk(i)='d'
            end if
         end do
         write(u6,'(I4,2X,100A1)') icsf,(walk(i),i=1,nopen)
      end do
#endif

#ifdef _TEST_CODE_
      if (ms2 == mult-1) then
         call check_twos( nopen, ndet, detspin, ncsf, twos )
      endif
#endif

      call mma_allocate(cdc,nDet,nCSF,Label='cdc')
      call sg_csfdet(nopen,ndet,detspin,ncsf,twos,cdc)

      call mma_deallocate(detspin)
      call mma_deallocate(twos)
      call mma_deallocate(cdc)

   end do

end do

end subroutine test_sg_csfdet
#endif


subroutine generate_detspin(nopen,ms2,detspin,ndet)

use Definitions, only: iwp
use stdalloc,   only: mma_allocate

implicit none

integer(iwp), intent(in) :: nopen
integer(iwp), intent(in) :: ms2

integer(iwp), intent(out) :: ndet
integer(iwp), allocatable, intent(out) :: detspin(:,:)

integer(iwp) :: ibit
integer(iwp) :: idet
integer(iwp) :: i
integer(iwp) :: mx

!
! Count determinants
!

ndet = 0
mx   = 2**nopen

do ibit=0,mx-1

   if (2*popcnt(ibit)-nopen == ms2) then
      ndet = ndet + 1
   endif

enddo

call mma_allocate(detspin,[1,Max(1,nopen)],[1,ndet],Label='detspin')

!
! Build determinant list
!

idet = 0

do ibit=0,mx-1

   if (2*popcnt(ibit)-nopen /= ms2) cycle

   idet = idet + 1

   do i=1,nopen
      detspin(i,idet) = ibits(ibit,i-1,1)
   enddo

enddo

end subroutine generate_detspin


subroutine generate_twos(nopen,mult,twos,ncsf)

use Definitions, only: iwp
use stdalloc,   only: mma_allocate

implicit none

integer(iwp), intent(in) :: nopen
integer(iwp), intent(in) :: mult

integer(iwp), intent(out) :: ncsf
integer(iwp), allocatable, intent(out) :: twos(:,:)

integer(iwp) :: ibit
integer(iwp) :: i
integer(iwp) :: ms2l
integer(iwp) :: mx
integer(iwp) :: icsf

logical :: upper

!
! Count walks
!

ncsf = 0
mx   = 2**nopen

do ibit=0,mx-1
   ms2l  = 0
   upper = .true.
   do i=1,nopen
      if (ibits(ibit,nopen-i,1) == 1) then
         ms2l = ms2l + 1
      else
         ms2l = ms2l - 1
      endif

      if (ms2l < 0) then
         upper = .false.
         exit
      endif
   enddo

   if (.not.upper) cycle
   if (ms2l /= mult-1) cycle
   ncsf = ncsf + 1

enddo

call mma_allocate(twos,[0,nopen],[1,ncsf],Label='twos')

!
! Build walks
!

icsf = 0
do ibit=0,mx-1

   ms2l  = 0
   upper = .true.
   do i=1,nopen
      if (ibits(ibit,i-1,1) == 1) then
         ms2l = ms2l + 1
      else
         ms2l = ms2l - 1
      endif

      if (ms2l < 0) then
         upper = .false.
         exit
      endif
   enddo

   if (.not.upper) cycle
   if (ms2l /= mult-1) cycle

   icsf = icsf + 1
   twos(0,icsf) = 0
   ms2l = 0

   do i=1,nopen
      if (ibits(ibit,i-1,1) == 1) then
         ms2l = ms2l + 1
      else
         ms2l = ms2l - 1
      endif

      twos(i,icsf) = ms2l
   enddo
enddo

end subroutine generate_twos

#ifdef _TEST_CODE_
subroutine check_twos(nopen, ndet, detspin, ncsf, twos)

use Definitions, only: iwp, u6

implicit none

integer(iwp), intent(in) :: nopen
integer(iwp), intent(in) :: ndet
integer(iwp), intent(in) :: ncsf

integer(iwp), intent(in) :: detspin(nopen, ndet)
integer(iwp), intent(in) :: twos(0:nopen, ncsf)

integer(iwp) :: idet
integer(iwp) :: iup
integer(iwp) :: i
integer(iwp) :: ms2l
integer(iwp) :: step

logical :: upper

iup = 0

do idet = 1, ndet

   upper = .true.
   ms2l  = 0

   do i = 1, nopen

      if (detspin(i, idet) == 1) then
         ms2l = ms2l + 1
      else
         ms2l = ms2l - 1
      end if

      if (ms2l < 0) then
         upper = .false.
         exit
      end if

   end do

   if (.not. upper) cycle

   iup = iup + 1

   if (iup > ncsf) then
      write(u6,*) 'CHECK_TWOS FAILED'
      write(u6,*) 'Too many upper determinants'
      write(u6,*) 'IUP  = ', iup
      write(u6,*) 'NCSF = ', ncsf
      call Abend()
   end if

   do i = 1, nopen

      step = twos(i,iup) - twos(i-1,iup)

      if (detspin(i,idet) == 1) then

         if (step /= 1) then
            write(u6,*) 'CHECK_TWOS FAILED'
            write(u6,*) 'IDET = ', idet
            write(u6,*) 'IUP  = ', iup
            write(u6,*) 'POSITION = ', i
            call Abend()
         end if

      else

         if (step /= -1) then
            write(u6,*) 'CHECK_TWOS FAILED'
            write(u6,*) 'IDET = ', idet
            write(u6,*) 'IUP  = ', iup
            write(u6,*) 'POSITION = ', i
            call Abend()
         end if

      end if

   end do

end do

if (iup /= ncsf) then
   write(u6,*) 'CHECK_TWOS FAILED'
   write(u6,*) 'Number of upper determinants /= number of CSFs'
   write(u6,*) 'IUP  = ', iup
   write(u6,*) 'NCSF = ', ncsf
   call Abend()
end if

end subroutine check_twos
#endif

#ifdef _GOSPER_
!Todo: implement in generate_detspin and generate_twos
subroutine next_combination(x)
   implicit none
   integer, intent(inout) :: x

   integer :: c, r

   c = iand(x, -x)
   r = x + c
   x = ior(r, ishft(ieor(r,x)/c, -2))

end subroutine next_combination
integer :: x, m, n

m = 5
n = 3

x = 2**n - 1

do while (x < 2**m)

   write(u6,'(I4,2X,B5.5)') x, x

   call next_combination(x)

end do
#endif

!***********************************************************************
! Generate the next orbital occupation configuration in Lucia order.
!
! Occupations are represented by the occupied orbital index of each
! electron:
!
!    OC = (1,1,2,2,3)
!
! corresponds to
!
!    DOC = (1,2)
!    SOC = (3)
!
! The occupation list is stored in nondecreasing order and no orbital
! may occur more than twice. With INI=1 the first configuration is
! generated. With INI=0 the next configuration in lexical Lucia order
! is generated.
!
! Input:
!    nElec  - number of electrons
!    nLev   - number of active orbitals
!    Ini    - 1: initialize, 0: increment
!
! Output:
!    OC     - orbital occupation configuration
!    NoNew  - 0: new configuration generated
!             1: no further configuration exists
!
!***********************************************************************
Subroutine Nxt_OC( OC, nElec, nLev, Ini, NoNew )

Implicit None

Integer, Intent(In)    :: nElec
Integer, Intent(In)    :: nLev
Integer, Intent(In)    :: Ini
Integer, Intent(InOut) :: OC(nElec)
Integer, Intent(Out)   :: NoNew

Integer :: I
Integer :: IAdd
Integer :: IEL
Integer :: Increase
Integer :: JOrb
Integer :: nDouble

If (Ini == 1) Then
   If (nElec <= 2*nLev) Then
      NoNew = 0
      nDouble = nElec/2
      Do I=1,nDouble
         OC(2*I-1) = I
         OC(2*I  ) = I
      End Do
      If (2*nDouble /= nElec) &
         OC(2*nDouble+1) = nDouble + 1
   Else
      NoNew = 1
   End If
Else
   IAdd = 1
   IEL  = 0
   Do While (IAdd == 1)
      IEL = IEL + 1
      Increase = 0
      If (IEL < nElec) Then
         If (OC(IEL) < OC(IEL+1)-1) Then
            Increase = 1
         Else If (OC(IEL) == OC(IEL+1)-1) Then
            If (IEL == nElec-1) Then
               Increase = 1
            Else If (OC(IEL+1) /= OC(IEL+2)) Then
               Increase = 1
            End If
         End If
      Else If (IEL == nElec) Then
         If (OC(IEL) < nLev) Then
            Increase = 1
         Else
            NoNew = 1
            Exit
         End If
      Else
         NoNew = 1
         Exit
      End If

      If (Increase == 1) Then
         NoNew  = 0
         OC(IEL)= OC(IEL) + 1
         nDouble = (IEL-1)/2
         Do JOrb=1,nDouble
            OC(2*JOrb-1) = JOrb
            OC(2*JOrb  ) = JOrb
         End Do
         If (2*nDouble < IEL-1) &
            OC(IEL-1) = nDouble + 1
         IAdd = 0
      End If
   End Do
End If

End Subroutine Nxt_OC

Subroutine Decode_OC(OC,nElec,nDOC,DOCOcc,nSOC,SOCOcc)

Implicit None

Integer, Intent(In)  :: nElec
Integer, Intent(In)  :: OC(nElec)

Integer, Intent(Out) :: nDOC
Integer, Intent(Out) :: nSOC
Integer, Intent(Out) :: DOCOcc(:)
Integer, Intent(Out) :: SOCOcc(:)

Integer :: i

nDOC = 0
nSOC = 0
i    = 1

Do While (i <= nElec)

   If (i < nElec) Then
      If (OC(i) == OC(i+1)) Then
         nDOC = nDOC + 1
         DOCOcc(nDOC) = OC(i)
         i = i + 2
      Else
         nSOC = nSOC + 1
         SOCOcc(nSOC) = OC(i)
         i = i + 1
      End If
   Else
      nSOC = nSOC + 1
      SOCOcc(nSOC) = OC(i)
      i = i + 1
   End If

End Do

End Subroutine Decode_OC


! PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
!          AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
!          RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
!          VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
!          CORRESPONDING STEP VECTOR.
subroutine STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)

use Definitions, only: iwp

implicit none

integer(iwp), intent(in)  :: NCL, NOP, NORB
integer(iwp), intent(in)  :: ICL(NCL), IOP(NOP)
integer(iwp), intent(in)  :: ISPIN(NOP)
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
enddo
!
!  Mark singly occupied orbitals
!
do iOpen=1,NOP
   if (ISPIN(iOpen) == 1) then
      IWALK(IOP(iOpen)) = 1
   else
      IWALK(IOP(iOpen)) = 2
   endif
enddo

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
! level (step d = 0,1,2,3 gives occupation 0,1,1,2), nPack levels per
! integer. Same layout as read in the diagonal case of MatEl_Epq.
!
subroutine CSF_Occupation(ISY,I,OccVec)
  integer(kind=iwp), intent(in)  :: ISY, I
  integer(kind=iwp), intent(out) :: OccVec(:)
  integer(kind=iwp) :: MV, ISYU, ISYD, IUP, IDN, iLev, IOW0, IPSHFT, IC, ICS
  call Decode_ME(ISY,I,MV,ISYU,IUP,IDN)
  ISYD = Mul(ISY,ISYU)
  do iLev = 1, SGS(ME_iState)%nLev
    if (iLev > SGS(ME_iState)%MidLev) then
      IOW0   = CIS(ME_iState)%IOW(1,ISYU,MV)
      IPSHFT = 2*(iLev-1-SGS(ME_iState)%MidLev)
      IC     = CIS(ME_iState)%ICase(1+IOW0-CIS(ME_iState)%nIpWlk+IPSHFT/(2*nPack)+IUP*CIS(ME_iState)%nIpWlk)
    else
      IOW0   = CIS(ME_iState)%IOW(2,ISYD,MV)
      IPSHFT = 2*(iLev-1)
      IC     = CIS(ME_iState)%ICase(1+IOW0-CIS(ME_iState)%nIpWlk+IPSHFT/(2*nPack)+IDN*CIS(ME_iState)%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*nPack)
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
      IC = CIS%ICase(1+IOW0-CIS%nIpWlk+IPSHFT/(2*nPack)+IUP*CIS%nIpWlk)
    else
      IOW0 = CIS%IOW(2,ISYD,MV)
      IPSHFT = 2*(iLev-1)
      IC = CIS%ICase(1+IOW0-CIS%nIpWlk+IPSHFT/(2*nPack)+IDN*CIS%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*nPack)
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
  if (.not. allocated(CIS(iA)%ICase)) call MkCOT(SGS(iA),CIS(iA))
  if (.not. allocated(CIS(iB)%ICase)) call MkCOT(SGS(iB),CIS(iB))
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

!> The walk tables of state iState for the CSFs of symmetry LSM: ICase
!> (MkCOT) and the walk numbering USGN/LSGN (MkSgNum), which depends on the
!> symmetry and is rebuilt when another symmetry is asked for (a state can
!> serve several symmetries, e.g. MCLR's state for the whole run).
subroutine SG_Num_Setup(iState,LSM)
  integer(kind=iwp), intent(in) :: iState, LSM
  if (.not. allocated(CIS(iState)%ICase)) call MkCOT(SGS(iState),CIS(iState))
  if ((.not. allocated(EXS(iState)%USGN)) .or. (EXS(iState)%Num_Sym /= LSM)) &
    call MkSgNum(LSM,SGS(iState),CIS(iState),EXS(iState))
end subroutine SG_Num_Setup

!> Set up the single-CSF matrix elements for the CSFs of symmetry LSM of
!> state iState (SG_CSF_Occ, SG_MatEl_H).
subroutine SG_MatEl_Init(iState,LSM)
  integer(kind=iwp), intent(in) :: iState, LSM
  ME_iState = iState
  ME_LSM = LSM
  call Reset_Decode_ME()
  call SG_Num_Setup(iState,LSM)
  if (.not. allocated(EXS(iState)%IOrd)) call Mk_Coup_Index(EXS(iState))
end subroutine SG_MatEl_Init

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
