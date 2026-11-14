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
module SGUGA

use Molcas, only: MxGas, MxLev, MxSym
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6
use CI_Verify_Options, only: CI_Verify

#include "intent.fh"

implicit none
private

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
type SGStruct
  integer(kind=iwp) :: NSym = 0, nActEl = 0, IFRAS = 0
  integer(kind=iwp) :: IA0, IB0, IC0, iSpin, nLev, nVert, nVert0, MidLev, MVSta, MVEnd, MXUP, MXDWN
  integer(kind=iwp), allocatable :: ISm(:), DRT(:,:), DRT0(:,:), Down(:,:), Down0(:,:), Up(:,:), Ver(:), MAW(:,:), LTV(:), &
                                    DAW(:,:), RAW(:,:)
  integer(kind=iwp), pointer :: DRTP(:,:), DOWNP(:,:)
  integer(kind=iwp) :: L2ACT(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: LEVEL(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: NRAS(MxSym,MxGAS), NRASEL(MxGAS), nRsPrt=0
  type(CDCInfo), allocatable :: CDCInfo(:,:)
  integer(kind=iwp), Allocatable :: OC(:)
  type(SDReoInfo), allocatable :: SDREO(:)
end type SGStruct

! index for DRT
integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2, ATAB = 3, BTAB = 4, CTAB = 5

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
type CIStruct
  integer(kind=iwp) :: nMidV, nIpWlk, nWalk, NUW
  integer(kind=iwp), allocatable :: NOW(:,:,:), IOW(:,:,:), NCSF(:), NOCSF(:,:,:), IOCSF(:,:,:), ICase(:), IVR(:,:), ISGM(:,:)
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
type EXStruct
  integer(kind=iwp) :: MxEO, nICoup
  integer(kind=iwp), allocatable :: NOCP(:,:,:), IOCP(:,:,:), ICoup(:,:), MVL(:,:), MVR(:,:), USGN(:,:), LSGN(:,:), IOrd(:,:)
  real(kind=wp), allocatable :: VTab(:), SGTMP(:)
  logical(kind=iwp) :: Reuse_SGTMP
  integer(kind=iwp), allocatable:: I1list(:), I2list(:)
  real(kind=wp), allocatable:: XList(:)
end type EXStruct

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

integer(kind=iwp), parameter ::                                                                                                    &
                                ITVPT(nSegTot) =  &
[ 0, 0, 0, 0,  0, 0, 0, 0,  1, 1, 1, 1, 1,  2, 2, 2, 2, 2,  1, 1, 2, 2,  3, 3, 3, 3,  0, 0, 0, 0,  3, 3, 3, 3],&
                                IBVPT(nSegTot) =  &
[ 0, 0, 0, 0,  1, 1, 2, 2,  1, 1, 2, 1, 1,  2, 2, 1, 2, 2,  3, 3, 3, 3,  3, 3, 3, 3,  0, 0, 0, 0,  3, 3, 3, 3],&
                                IC1(nSegTot)   =  &
[ 0, 1, 2, 3,  0, 2, 0, 1,  0, 1, 1, 2, 3,  0, 1, 2, 2, 3,  1, 3, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3],&
                                IC2(nSegTot)   =  &
[ 0, 1, 2, 3,  1, 3, 2, 3,  0, 1, 2, 2, 3,  0, 1, 1, 2, 3,  0, 2, 0, 1,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3],&
                                ISVC(nSegTot)  =  &
[ 1, 1, 1, 1,  1, 7, 8, 4,  1, 2, 9,10, 2,  1, 2,11,12, 2,  1, 5, 6, 3,  1, 1, 1, 1,  0, 1, 1,13,  0, 1, 1,13]

integer(kind=iwp), parameter :: TR_WALK  = 1
integer(kind=iwp), parameter :: TR_OPEN  = 2
integer(kind=iwp), parameter :: TR_MID   = 4
integer(kind=iwp), parameter :: TR_CLOSE = 8
integer(kind=iwp), parameter :: TR_TAIL  = 16
integer(kind=iwp), parameter :: TR_WEIGHT = 32

public :: SGStruct, CIStruct, EXStruct, SG_Free, SG_Init, SG_Init_Light, SG_Epq_Psi, SG_ReOrd, MkCOT, MkSgNum, SG_NUM, SG_PHASE
public :: SG_Lucia_Walk
public :: SG_EpqEqp_Diag, sg_csfdet, generate_detspin, generate_twos, SG_ReOrd_Orig, StepVec, SG_TRACI
public :: SG_Test_EpqErs_Ref
public :: SG_Mk_SDREO
!#define _TEST_CSFDET_
#ifdef _TEST_CSFDET_
public :: test_sg_csfdet
#endif

! Set nPack to the number of cases (2 bit per case) that can be packed in one integer.
#ifdef SIZE_INITIALIZATION
integer(kind=iwp), parameter :: nPack = storage_size(1_iwp)/2-1
#elif defined (_I8_)
integer(kind=iwp), parameter:: nPack=32-1
#else
integer(kind=iwp), parameter:: nPack=16-1
#endif

public :: nPack

Type (SGStruct), target :: SGS(3)
Type (CIStruct) :: CIS(3)
Type (EXStruct) :: EXS(3)
Logical(kind=iwp) :: State_is_used(3)=[.False.,.False.,.False.]

Public:: SGS, CIS, EXS, State_is_used

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

public :: SG_CSF_Occ, SG_MatEl_H, SG_MatEl_Init


contains

subroutine MKSGUGA(SGS,CIS)
! PURPOSE: MAKE THE GUGA TABLES
! NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!          THE START ADDRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!          THREE USER DEFINED TYPES. Consult the sguga module for the details.

  type(SGStruct), target, intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS

  ! COMPUTE TOP ROW OF THE GUGA TABLE

  call mknVert0(SGS)

  ! SET UP A FULL PALDUS DRT TABLE:
  ! (INITIALLY NO RESTRICTIONS ARE PUT UP)

  SGS%nVert = SGS%nVert0
  if (SGS%IFRAS /= 0) then
    call mma_allocate(SGS%DRT0,SGS%nVert0,5,Label='DRT0')
    call mma_allocate(SGS%DOWN0,[1,SGS%nVert0],[0,3],Label='DOWN0')
    SGS%DRTP => SGS%DRT0
    SGS%DOWNP => SGS%DOWN0
  else
    call mma_allocate(SGS%DRT,SGS%nVert,5,Label='SGS%DRT')
    call mma_allocate(SGS%DOWN,[1,SGS%nVert],[0,3],Label='SGS%DOWN')
    SGS%DRTP => SGS%DRT
    SGS%DOWNP => SGS%DOWN
  end if

  call mkDRT0(SGS)

  ! IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
  ! VERTICES WHICH VIOLATE THE FORMER.

  if (SGS%IFRAS /= 0) then
    call RmVert(SGS) ! REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
    call mkDRT(SGS) ! IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED DRT TABLE
  end if

  nullify(SGS%DOWNP,SGS%DRTP)

  call MKDAW(SGS) ! COMPUTE DOWNCHAIN TABLE AND ARC WEIGHT
  call MKRAW(SGS) ! COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
  call MKLTV(SGS) ! COMPUTE LTV TABLES.
  call MKMID(SGS) ! COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
  call MKCOT(SGS,CIS) ! CONSTRUCT THE CASE LIST
contains

  subroutine mknVert0(SGS)

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IAC

    SGS%IB0 = SGS%iSpin-1
    SGS%IA0 = (SGS%nActEl-SGS%IB0)/2
    SGS%IC0 = SGS%nLev-SGS%IA0-SGS%IB0

    if ((2*SGS%IA0+SGS%IB0 /= SGS%nActEl) .or. (SGS%IA0 < 0) .or. (SGS%IB0 < 0) .or. (SGS%IC0 < 0)) then
      write(u6,*) 'mknVert0 Error: Impossible specifications.'
      write(u6,'(1x,a,3I8)') 'NACTEL,NLEV,ISPIN:',SGS%nActEl,SGS%nLev,SGS%iSpin
      write(u6,'(1x,a,3I8)') 'IA0,IB0,IC0:      ',SGS%IA0,SGS%IB0,SGS%IC0
      write(u6,*) ' This is a severe internal error, or possibly'
      write(u6,*) ' indicates a strange input which should have been'
      write(u6,*) ' diagnosed earlier. Please submit a bug report.'
      call Abend()
    end if

!   See Eqs. (2) of 10.1002/jcc.26080

    IAC = min(SGS%IA0,SGS%IC0)
    SGS%nVert0 = ((SGS%IA0+1)*(SGS%IC0+1)*(2*SGS%IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

  end subroutine mknVert0

subroutine mkDRT0(SGS)
  type(SGStruct), target, intent(inout) :: SGS

  integer(kind=iwp) :: ADWN, BDWN, CDWN
  integer(kind=iwp) :: AUP, BUP, CUP
  integer(kind=iwp) :: LEV, ISTEP
  integer(kind=iwp) :: VUP, VSTA, VEND, VNEW, VDWN
  integer(kind=iwp) :: key, maxKey

  integer(kind=iwp), parameter :: Steps(4,0:3)=reshape( &
       [0,0,1,0, 0,1,0,1, 1,-1,1,1, 1,0,0,2],[4,4])

  ! --- bounds for encoding ---
  integer(kind=iwp) :: Amax, Bmax, Cmax
  integer(kind=iwp), allocatable :: vmap(:)

  SGS%DRTP(:,:) = 0
  SGS%DownP(:,:) = 0

  SGS%NACTEL = 2*SGS%IA0 + SGS%IB0
  SGS%nLev   = SGS%IA0 + SGS%IB0 + SGS%IC0

  ! maximum ranges (safe upper bounds)
  Amax = SGS%IA0
  Bmax = SGS%IB0 + SGS%IA0   ! loose but safe
  Cmax = SGS%IC0

  maxKey = (Amax+1)*(Bmax+1)*(Cmax+1)

  call mma_allocate(vmap,[0,maxKey-1],Label='vmap')
  vmap(:) = 0

  ! ------------------------------------------
  ! helper: inline key encoding
  ! ------------------------------------------
#   define _KEY_(A,B,C) ((A)*(Bmax+1)*(Cmax+1)+(B)*(Cmax+1)+(C))

  ! TOP vertex
  SGS%DRTP(1,LTAB) = SGS%nLev
  SGS%DRTP(1,NTAB) = SGS%NACTEL
  SGS%DRTP(1,ATAB) = SGS%IA0
  SGS%DRTP(1,BTAB) = SGS%IB0
  SGS%DRTP(1,CTAB) = SGS%IC0

  key = _KEY_(SGS%IA0,SGS%IB0,SGS%IC0)
  vmap(key) = 1

  VSTA = 1
  VEND = 1

  ! ------------------------------------------
  ! main loop
  ! ------------------------------------------
  do LEV = SGS%nLev, 1, -1
    VNEW = 0

    do VUP = VSTA, VEND
      AUP = SGS%DRTP(VUP,ATAB)
      BUP = SGS%DRTP(VUP,BTAB)
      CUP = SGS%DRTP(VUP,CTAB)

      do ISTEP = 0, 3

        ADWN = AUP - Steps(1,ISTEP)
        if (ADWN < 0) cycle

        BDWN = BUP - Steps(2,ISTEP)
        if (BDWN < 0) cycle

        CDWN = CUP - Steps(3,ISTEP)
        if (CDWN < 0) cycle

        key = _KEY_(ADWN,BDWN,CDWN)

        VDWN = vmap(key)

        if (VDWN /= 0) then
          SGS%DownP(VUP,ISTEP) = VDWN
        else
          VNEW = VNEW + 1
          VDWN = VEND + VNEW

          SGS%DownP(VUP,ISTEP) = VDWN

          SGS%DRTP(VDWN,ATAB) = ADWN
          SGS%DRTP(VDWN,BTAB) = BDWN
          SGS%DRTP(VDWN,CTAB) = CDWN

          vmap(key) = VDWN
        end if

      end do
    end do

    VSTA = VEND + 1
    VEND = VEND + VNEW

  end do

  ! ------------------------------------------
  ! finalize LTAB / NTAB
  ! ------------------------------------------
  SGS%DRTP(1:VEND,LTAB) = SGS%DRTP(1:VEND,ATAB) + &
                          SGS%DRTP(1:VEND,BTAB) + &
                          SGS%DRTP(1:VEND,CTAB)

  SGS%DRTP(1:VEND,NTAB) = 2*SGS%DRTP(1:VEND,ATAB) + &
                            SGS%DRTP(1:VEND,BTAB)

    call mma_deallocate(vmap)

#   undef _KEY_

end subroutine mkDRT0

  subroutine mkDRT(SGS)
  ! PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
  !          THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
  !          VERTICES FROM THE OLD TO THE NEW DRT TABLE

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, ID, IDNEW, IV, IVNEW

    call mma_allocate(SGS%DRT,SGS%nVert,5,Label='DRT')
    call mma_allocate(SGS%Down,[1,SGS%nVert],[0,3],Label='SGS%DOWN')

    ! Loop over all original vertices
    do IV=1,SGS%nVert0
      IVNEW = SGS%Ver(IV) ! Pick up the new index, if iVNEW==0
      if (IVNEW == 0) cycle
      ! Move over all relevant elements from DRT0
      SGS%DRT(IVNEW,:) = SGS%DRT0(IV,1:5)
      ! Loop over step vector
      do IC=0,3
        ID = SGS%Down0(IV,IC) ! Pick up the old index of the vertex which connects to vertex IV
        IDNEW = 0
        ! If a viable arc get the new index of the vertex ID
        if (ID /= 0) IDNEW = SGS%Ver(ID)
        ! Update the downward Chain
        SGS%Down(IVNEW,IC) = IDNEW
      end do
    end do
#   ifdef _DEBUGPRINT_
    do IV=1,SGS%nVert
      write(u6,*) 'DRT(i,:)=',SGS%DRT(IV,:)
    end do
#   endif

    call mma_deallocate(SGS%Ver)
    call mma_deallocate(SGS%DRT0)
    call mma_deallocate(SGS%Down0)

  end subroutine mkDRT

  subroutine MKDAW(SGS)
  ! PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, IDWN, IV

    call mma_allocate(SGS%DAW,[0,SGS%nVert],[0,4],Label='SGS%DAW')

    ! Note that this is done without symmetry.
    ! BEGIN TO CONSTRUCT DOWN CHAIN TABLE

    SGS%DAW(:,:) = 0
    SGS%DAW(SGS%nVert,4) = 1

    ! Loop over all vertices
    do IV=SGS%nVert-1,1,-1
      do IC=0,3
        IDWN = SGS%Down(IV,IC)
        SGS%DAW(IV,IC+1) = SGS%DAW(IV,IC) + SGS%DAW(IDWN,4)
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' DIRECT ARC WEIGHTS:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%DAW(IV,0:4)
    end do
    write(u6,*)
#   endif

  end subroutine MKDAW

  subroutine MKRAW(SGS)
  ! PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, IDWN, IU, IV

    call mma_allocate(SGS%UP,[1,SGS%nVert],[0,3],Label='SGS%UP')
    call mma_allocate(SGS%RAW,[0,SGS%nVert],[0,4],Label='SGS%RAW')

    ! BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:

    SGS%Up(:,:) = 0
    do IU=1,SGS%nVert-1
      do IC=0,3
        IDWN = SGS%Down(IU,IC)
        if (IDWN == 0) cycle
        SGS%Up(IDWN,IC) = IU
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' THE UPCHAIN TABLE IN MKRAW:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,4(1X,I6))') IV,SGS%Up(IV,0:3)
    end do
    write(u6,*)
#   endif

    ! USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:

    SGS%RAW(:,:) = 0
    SGS%RAW(1,4) = 1
    do IV=2,SGS%nVert
      do IC=0,3
        IU = SGS%Up(IV,IC)
        SGS%RAW(IV,IC+1) = SGS%RAW(IV,IC) + SGS%RAW(IU,4)
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' THE REVERSE ARC WEIGHT TABLE IN MKRAW:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%RAW(IV,0:4)
    end do
    write(u6,*)
#   endif

  end subroutine MKRAW

  subroutine MKLTV(SGS)

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IV, LEV

    !Note that the range is [-1,nLev]! An imaginary level, -1, is added.
    call mma_allocate(SGS%LTV,[-1,SGS%nLev],Label='LTV')

    ! SET UP A LEVEL-TO-VERTEX TABLE, LTV

    SGS%LTV(:) = 0

    ! Loop over all vertices
    do IV=1,SGS%nVert
      LEV = SGS%DRT(IV,LTAB) ! Pick up the level index of the vertex
      SGS%LTV(LEV) = SGS%LTV(LEV)+1 ! Increment the number of vertex for a particular level
    end do

    ! Loop over all levels, from the level corresponding to the head vertex to the level corresponding
    ! to the root vertex of the graph
    do LEV=SGS%nLev,0,-1
      SGS%LTV(LEV-1) = SGS%LTV(LEV-1)+SGS%LTV(LEV) ! Update such that  SGS%LEV(Lev) is the number of vertex at level
                                                   ! Lev and above.
    end do

    do LEV=-1,SGS%nLev-1
      SGS%LTV(LEV) = 1+SGS%LTV(LEV+1)
    end do

    ! SGS%LTV(LEV) now tabulates the index of the first vertex at level LEV if the vertices are numbered 1,..,nVert
    ! starting from the head vertex.

  end subroutine MKLTV

  subroutine MKMID(SGS)
  ! PURPOSE: FIND THE MIDLEVEL

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IL, IV, MINW, MV, NW

    ! USE DAW,RAW TABLES TO DETERMINE MIDLEV.
    ! THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
    ! IS THE BEST CHOICE.

    !hrl 980529 fix for nLev=0 (no orbitals in any active space)
    ! Since LTV(-1:nLev)  and the statement after the loop
    ! MVSta=LTV(MidLev) we have the condition MidLev>=nLev
    ! Hence MidLev=1 is inappropriate for nLev=0
    ! MIDLEV=1

    SGS%MidLev = min(SGS%nLev,1)

    MINW = 1000000 ! Initiate to a large number
    do IL=1,SGS%nLev-1 ! Loop from level 1 (one above the root vertex) to the
                       ! level below the head vertex.

      ! SGS%RAW(IV,4) is the number of possible paths from the root vertex up to vertex IV
      ! SGS%DAW(IV,4) is the number of possible paths from the head vertex down to vertex IV
      NW = 0
      do IV=SGS%LTV(IL),SGS%LTV(IL-1)-1  ! Loop over all vertex at this level
        NW = NW+abs(SGS%RAW(IV,4)-SGS%DAW(IV,4))
      end do
      NW = abs(NW)
      if (NW >= MINW) cycle
      SGS%MidLev = IL        ! Update level index indicating the midlevel
      MINW = NW              ! Update MINW
    end do

    SGS%MVSta = SGS%LTV(SGS%MidLev)      ! Staring vertex index at level MIDLEV
    SGS%MVEnd = SGS%LTV(SGS%MidLev-1)-1  ! Ending vertex index at level MIDLEV
    CIS%nMidV = SGS%MVEnd-SGS%MVSta+1    ! Number of vertices at level MIDLEV

    ! NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
    ! (DISREGARDING SYMMETRY)

    SGS%MxUp = 0
    SGS%MxDwn = 0
    do MV=SGS%MVSta,SGS%MVEnd ! Lopp over all vertices of level MIDLEV
      SGS%MxUp  = max(SGS%MxUp,SGS%RAW(MV,4))
      SGS%MxDwn = max(SGS%MxDwn,SGS%DAW(MV,4))
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,'(A,I3)') ' MIDLEVEL =             ',SGS%MidLev
    write(u6,'(A,I3)') ' NUMBER OF MIDVERTICES =',CIS%NMIDV
    write(u6,'(A,I3)') ' FIRST MIDVERTEX =      ',SGS%MVSta
    write(u6,'(A,I3)') ' LAST MIDVERTEX =       ',SGS%MVEnd
    write(u6,'(A,I3)') ' MAX. NO UPPER WALKS=   ',SGS%MxUp
    write(u6,'(A,I3)') ' MAX. NO LOWER WALKS=   ',SGS%MxDwn
    write(u6,*)
#   endif

  end subroutine MKMID

subroutine RMVERT(SGS)

  type(SGStruct), intent(inout) :: SGS

  integer(kind=iwp) :: IV, IC, ID
  integer(kind=iwp) :: L, N
  integer(kind=iwp) :: NCHANGES, NV, iRO, iSy, Lev
  integer(kind=iwp), allocatable :: CONN(:), Lim(:)

  !-----------------------------------------
  ! Build occupation limits
  !-----------------------------------------
  call mma_allocate(Lim,SGS%nLev,Label='Lim')
  Lim(:) = 0

  Lev = 0
  do iRO = 1, SGS%nRsPrt
    do iSy = 1, SGS%nSym
      Lev = Lev + SGS%nRas(iSy,iRO)
    end do
    if (Lev > 0) Lim(Lev) = SGS%nRasEl(iRO)
  end do

  !-----------------------------------------
  ! Initialize vertex mask
  !-----------------------------------------
  call mma_allocate(SGS%Ver,SGS%nVert0,Label='SGS%Ver')
  SGS%Ver(:) = 1

  call mma_allocate(CONN,SGS%nVert,Label='CONN')

  !-----------------------------------------
  ! Initial pruning based on occupation
  !-----------------------------------------
  do IV = 1, SGS%nVert-1
    L = SGS%DRT0(IV,LTAB)
    N = SGS%DRT0(IV,NTAB)

    if (N < Lim(L)) SGS%Ver(IV) = 0
  end do

  !-----------------------------------------
  ! Iterative pruning
  !-----------------------------------------
  do
    NCHANGES = 0

    !-------------------------------------
    ! Step 1: remove invalid downward arcs
    !-------------------------------------
    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) then
        SGS%Down0(IV,0:3) = 0
        cycle
      end if

      do IC = 0,3
          ID = SGS%Down0(IV,IC)
          if (ID < 1) cycle

        if (SGS%Ver(ID) == 0) then
           SGS%Down0(IV,IC) = 0
           NCHANGES = NCHANGES + 1
        end if
      end do
    end do

    !-------------------------------------
    ! Step 2: remove vertices with no children
    !-------------------------------------
    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) cycle

      if (all(SGS%Down0(IV,0:3) == 0)) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES + 1
      end if
    end do

    !-------------------------------------
    ! Step 3: rebuild upward connectivity
    !-------------------------------------
    CONN(:) = 0
    CONN(1) = SGS%Ver(1)

    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) cycle

      do IC = 0,3
          ID = SGS%Down0(IV,IC)
          if (ID < 1) cycle
          if (SGS%Ver(ID) == 1) CONN(ID) = 1
      end do
    end do

    ! remove vertices not connected from above
    do IV = 1, SGS%nVert
        if ((SGS%Ver(IV) == 1) .and. (CONN(IV) == 0)) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES + 1
      end if
    end do

    if (NCHANGES == 0) exit

  end do

  !-----------------------------------------
  ! Safety check
  !-----------------------------------------
  if (SGS%Ver(1) == 0) then
    write(u6,*) 'RASSI/RMVERT: Too severe restrictions.'
    call Abend()
  end if

  !-----------------------------------------
  ! Renumber vertices
  !-----------------------------------------
  NV = 0
  do IV = 1, SGS%nVert
    if (SGS%Ver(IV) == 1) then
      NV = NV + 1
      SGS%Ver(IV) = NV
    end if
  end do

  SGS%nVert = NV

  call mma_deallocate(CONN)
  call mma_deallocate(Lim)

end subroutine RMVERT

end subroutine MKSGUGA

subroutine SG_Init(iState,nSym,nActEl,iSpin,                      &
                   nRas,nRasEl,nRsPrt,                            &
                   xLevel,xL2Act,xNLEV,xNSM)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt),nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xnLev, xNSM(MxLev)

  type(TRStruct) :: TRS

  call SG_Init_Light(iState,nSym,nActEl,iSpin,        &
                      nRas,nRasEl,nRsPrt,             &
                      xLevel,xL2Act,xnLev,xNSM)

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
                          xLevel,xL2Act,xNLEV,xNSM,Do_MkSGUGA)

  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xNLEV, xNSM(MxLev)
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
 call mma_deallocate(EXS%MVL,safe='*')
 call mma_deallocate(EXS%MVR,safe='*')
 call mma_deallocate(EXS%USGN,safe='*')
 call mma_deallocate(EXS%LSGN,safe='*')
 call mma_deallocate(EXS%I1list,safe='*')
 call mma_deallocate(EXS%I2list,safe='*')
 call mma_deallocate(EXS%Xlist,safe='*')

End Associate

State_is_used(iState)=.False.

end subroutine SG_Free

subroutine Mk_CDC(SGS)

use second_quantization, only: binomial=>binom_coef

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
      if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

      if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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


   if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

   if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

   if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

      if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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


              if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
              If (SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC) Then
                 Write (u6,'(A)') 'SGS%CDCInfo(nOpen,Mult)%OCReverse(iDOC,iSOC,iRefSm) /= iOC'
                 Call Abend()
              End If
              end if

              Call FindSpinIndex(OpenSP,                                               &
                                 iOpen,                                                &
                                 SGS%CDCInfo(iOpen,mult)%DetSpin,                      &
                                 iSpinPattern)

              if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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
              if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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
     if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
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

use second_quantization, only: binomial=>binom_coef
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

subroutine MKCOT(SGS,CIS)
! For INIT==0
! PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
! NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
!          THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
!          WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
!          SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
!          THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.

! For INIT==1
! PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
!          STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
!          IN THE ARRAY ICASE. GROUPS OF nPack CASES ARE PACKED
!          INTO ONE INTEGER WORD.

type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp) :: IC, IHALF, ILND, INIT, IPOS, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, L, LEV, LEV1, LEV2, LL, MV
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: IS, NUW
# endif
integer(kind=iwp), allocatable :: SCR(:,:)
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

do INIT=0,1
  if (INIT == 0) then
    ! Compute the number of integers needed to store the packed case vectors.
    CIS%nIpWlk = 1+(SGS%MidLev-1)/nPack
    CIS%nIpWlk = max(CIS%nIpWlk,1+(SGS%nLev-SGS%MidLev-1)/nPack)
    call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW')
    call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW')
    call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF')
    call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF')
    call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF')

    call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr')  ! First index referenced by IVERT, SYM, and ISTEP

    ! CLEAR ARRAYS IOW AND NOW

    CIS%NOW(:,:,:) = 0
    CIS%IOW(:,:,:) = 0

    ! CLEAR ARRAYS IOCSF AND NOCSF

    CIS%IOCSF(:,:,:) = 0
    CIS%NOCSF(:,:,:) = 0
  else
    !call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr',safe='*')  ! First index referenced by IVERT, SYM, and ISTEP
    call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')
    ! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY
    CIS%NOW(:,:,:) = 0
  end if

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
  ! set the loop ranges:
  ! IVSTA-IVTEND: the top vertex, or loop over midlevel vertices
  ! LEV2-LEV1: loop over levels
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = SGS%nLev
    LEV2 = SGS%MidLev
  else
    IVTSTA = SGS%MVSta
    IVTEND = SGS%MVEnd
    LEV1 = SGS%MidLev
    LEV2 = 0
  end if

  ! LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH
  ! This is either the top vertex, or one of the midlevel vertices

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
    LEV = LEV1

    ! Store away vertex index, and initiate symmetry index, and the step vector index
    Scr(IVERT,LEV) = IVTOP
    Scr(ISYM,LEV)  =  1
    Scr(ISTEP,LEV) = -1

    do while (LEV <= LEV1)

      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
      IVT = Scr(IVERT,LEV)  ! Pickup the current vertex index for level LEV
      ! Continue scanning the step vectors, SCR(ISTEP,LEV) is the index of the next vector to explot
      do ISTP=Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) exit      ! Exits if step vector leads to a valid vertex below.
      end do

      ! IF NO SUCH ARC WAS POSSIBLE. GO UP ONE LEVEL AND TRY AGAIN.
      if (ISTP > 3) then
        Scr(ISTEP,LEV) = -1
        LEV = LEV+1
        cycle
      end if

      ! SUCH AN ARC WAS FOUND. WALK DOWN:
      Scr(ISTEP,LEV) = ISTP ! Store the current step vector index of level Level
      ! doubly occupied or empty orbital case are total symmetric. Singly occupied orbitals
      ! carry the symmetry of the orbital in the level.
      ISML = 1
          if (ISTP == 1 .or. ISTP == 2) ISML = SGS%ISm(LEV)

      LEV = LEV-1     ! Walk down one level

      ! Store away the accumulated symmetry, the new vertex, and initiate the step vector counter.
      Scr(ISYM,LEV) = Mul(ISML,Scr(ISYM,LEV+1))
      Scr(IVERT,LEV) = IVB
      Scr(ISTEP,LEV) = -1

      if (LEV > LEV2) cycle   ! Repeat

      ! WE HAVE NOW REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
      ! FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
      MV = Scr(IVERT,SGS%MidLev)+1-SGS%MVSta ! Pick up the relative index of the midlev vertex
      IWSYM = Scr(ISYM,LEV2)                 ! Pick up the symmetry of the walk

      ILND = CIS%NOW(IHALF,IWSYM,MV) + 1
      CIS%NOW(IHALF,IWSYM,MV) = ILND   ! Increment counter for how many walks there are given the
      ! midvertex index and the total symmetry of the subwalk.
      if (INIT == 1) then
        ! CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
        IPOS = CIS%IOW(IHALF,IWSYM,MV)+(ILND-1)*CIS%nIpWlk

        ! PACK THE STEPS IN GROUPS OF nPack LEVELS PER INTEGER:
        do LL=LEV2+1,LEV1,nPack

          IC = 0
          do L=min(LL+nPack-1,LEV1),LL,-1
            IC = 4*IC+Scr(ISTEP,L)
          end do

          IPOS = IPOS+1
          CIS%ICase(IPOS) = IC
        end do
      end if

      ! BACK UP ONE LEVEL AND EXPLORE NEW WALKS:
      LEV = LEV+1
    end do
  end do
end do

if (INIT == 0) then
! Generate an off-set array for the CIS%NOW array
call Mk_IOW(CIS,SGS)
call CSFCOUNT(CIS,SGS)

#ifdef _DEBUGPRINT_
NUW=CIS%NUW
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)
write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN MkNRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:',(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,'(A,8I6)') '           LOWER WALKS:',(CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do ISTP=1,SGS%nSym
    write(u6,'(A,I2,A,8I6)') ' ISTP=',ISTP,'  CONFIGURATIONS:',(CIS%NOCSF(IS,MV,ISTP),IS=1,SGS%nSym)
  end do
end do
#endif
end if

end do ! INIT

call mma_deallocate(Scr,safe='*')

end subroutine MKCOT

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
  integer(kind=iwp) :: I_N, IA, IAL, IB, IBL, INL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB, IVLT, IVRB, IVRT, LEV, MV, MVLL
real(kind=wp) :: V

call mma_allocate(CIS%IVR,SGS%nVert,2,Label='CIS%IVR')
call mma_allocate(CIS%ISGM,SGS%nVert,nSeg,Label='CIS%ISGM')
call mma_allocate(CIS%VSGM,SGS%nVert,nSeg,Label='CIS%VSGM')
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
  do ISGT=1,nSeg        ! Setment types

    ! Find all four vertex indices of this segment, branch out if the segment is invalid
    ITT = ITVPT(ISGT)   ! 0-3, denotes class of the upper nodes of the segment
    IVRT = IVLT               ! Vertex Right Top is the same as the Vertex Top Left, default setting
    if (ITT == 1 .or. ITT == 2) IVRT = CIS%IVR(IVLT,ITT)  ! Pick up the associated upper right vertex, depends on delta(b)
    if (IVRT == 0) cycle          ! Branch out if there is no vertex that will contruct the top of the segment.
    ! Pick up the vertex index for the Vertex Left Bottom as a function of the segment type
    IVLB = SGS%Down(IVLT,IC1(ISGT))
    if (IVLB == 0) cycle          ! Branch out if there is no such vertex
    ! Pick up the vertex index for the Vertex Right Bottom v as a function of the case
    IVRB = SGS%Down(IVRT,IC2(ISGT))
    if (IVRB == 0) cycle          ! Branch out if there is no such vertex
    ! SEGMENT IS NOW ACCEPTED AS POSSIBLE.

    CIS%ISGM(IVLT,ISGT) = IVLB  ! Mark that the segment is valid by changing the default value, 0,  to the index of
                                ! lower left vertex being a part of the segment.
    IB = SGS%DRT(IVLT,BTAB)    ! Pick up the b-value from the DRT
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
               If (IP==3 .and. IHALF==1) Write (6,*) 'IT,IP,INDEO=',IT,IP,INDEO
            endif

            ! Increment the local index and store
            ICOP = 1 + EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
            ! Add the offset and get the global index
            ICOP = EXS%IOCP(INDEO,LFTSYM,MV) + ICOP
            If (IP==3 .and. IT/=3 .and. IHALF==1) THEN
               Write (u6,*) 'INDEO,LFTSYM,MV=',INDEO,LFTSYM,MV
            End If

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
            If (IP==3 .and. IT/=3 .and. IHALF==1) Then
               Write (u6,*) 'ICOP,IVTAB=',ICOP,IVTAB
               Write (u6,*) 'Iu=',ISGPTH(IAWSL,LEV2)
               Write (u6,*) 'Ju=',ISGPTH(IAWSR,LEV2)
               Write (u6,*) 'A_r=',VTAB(IVTAB)
            End If

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

subroutine CSFCOUNT(CIS,SGS)
! CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
! SEPARATED BY MIDVERTICES AND SYMMETRY.
! FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.

type(CIStruct), intent(inout) :: CIS
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: ISYDWN, ISYTOT, ISYUP, MV, N

CIS%NCSF(:) = 0   ! Number of CSFs in each irrep
do ISYTOT=1,SGS%NSYM
  do MV=1,CIS%nMidV
    do ISYUP=1,SGS%NSYM
      ISYDWN = Mul(ISYTOT,ISYUP)
      N = CIS%NOW(1,ISYUP,MV)*CIS%NOW(2,ISYDWN,MV)
      CIS%NOCSF(ISYUP,MV,ISYTOT) = N
      CIS%IOCSF(ISYUP,MV,ISYTOT) = CIS%NCSF(ISYTOT)
      CIS%NCSF(ISYTOT) = CIS%NCSF(ISYTOT)+N
#     include "compiler_features.h"
#     ifdef _BUGGY_INTEL_LLVM_
      ! dummy statement to work around compiler bug, will never be executed
      if (ISYUP > 99) CIS%NCSF(ISYTOT) = -1
#     endif
    end do
  end do
end do

end subroutine CSFCOUNT

subroutine Mk_IOW(CIS,SGS)
! CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
! SEPARATED FOR EACH MIDVERTEX AND SYMMETRY
!
! CIS%IOW(1/2,ISYM,MV) is an off-set vector associated with CIS%NOW(1/2,ISYM,MV)
! CIS%NOW(1/2,ISYM,MV) is the number of upper/lower walks of symmetry ISYM that ends/starts in mid vertex MV (relative indexation).
! 1/2 indicate if the off-set is for upper or lower walks.

  type(CIStruct), intent(inout) :: CIS
  type(SGStruct), intent(inout) :: SGS
  integer(kind=iwp) :: ISYM, MV

CIS%NUW = 0   ! At completion the total number of upper walks.
do MV=1,CIS%nMidV
  do ISYM=1,SGS%NSYM
    CIS%IOW(1,ISYM,MV) = CIS%NUW*CIS%nIpWlk
    CIS%NUW = CIS%NUW+CIS%NOW(1,ISYM,MV)
  end do
end do
CIS%nWalk = CIS%NUW ! At completion the total number pf upper and lower walks

do MV=1,CIS%nMidV
  do ISYM=1,SGS%NSYM
    CIS%IOW(2,ISYM,MV) = CIS%nWalk*CIS%nIpWlk
    CIS%nWalk = CIS%nWalk+CIS%NOW(2,ISYM,MV)
  end do
end do
end subroutine Mk_IOW

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

subroutine SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM)

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
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2, JQ, JP
real(kind=wp) :: X

! declarations to facilitate the reuse of sigma vectors if possible
integer(kind=iwp), save:: i_save_p=0, i_save_q=0
integer(kind=iwp), save:: i_save_p_sym=-1, i_save_q_sym=-1
logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2, Test3, Moved
real(kind=wp) :: CI_ID=Zero
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
  Reuse_Sigma = .false.
  Test3=MERGE((i_save_q == IQ) .and. (i_save_p_sym == ISYP), &
              (i_save_p == IP) .and. (i_save_q_sym == ISYQ), Swap)
  if (EXS%Reuse_SGTMP .and.  Test3) Reuse_Sigma = CI_ID == sum(CI(1:nCSFs))+real(nCSFs,kind=wp)

  JP=MERGE(IQ,IP,SWAP)
  JQ=MERGE(IP,IQ,SWAP)
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
      If (Swap .and. MVX<MVSGM) Then
         Write (6,*) 'MV(1)=',MV
         Write (6,*) 'MVSGM=',MVSGM
         Write (6,*) 'MVX=',MVX
         Write (6,*) 'Urk'
         Call abend()
      End If
      If (.not.Swap .and. MVX>MVSGM) Then
         Write (6,*) 'MV(2)=',MV
         Write (6,*) 'MVSGM=',MVSGM
         Write (6,*) 'MVX=',MVX
         Write (6,*) 'Urk'
         Call abend()
      End If

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

  If (SWAP) Then
    i_save_q=IQ
    i_save_p_sym=ISYP
    i_save_p=0
    i_save_q_sym=-1
  else
    i_save_p=IP
    i_save_q_sym=ISYQ
    i_save_q=0
    i_save_p_sym=-1
  end if
  if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) CI_ID = sum(CI(1:nCSFs))+real(nCSFs,kind=wp)

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
  i_save_q=0
  i_save_p_sym=-1
  i_save_p=0
  i_save_q_sym=-1
  CI_ID=Zero
end if

contains

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

  integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUP,NDWNC)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp), parameter :: KBLOCK = 16
  integer(kind=iwp) :: blck, finish, finish2, i, I1a(KBLOCK), I1b(KBLOCK), I2, I2b, ICP, nk, nk2, offset, start, start2
  real(kind=wp) :: Xa(KBLOCK), Xb(KBLOCK), ABLOCK(NUP,KBLOCK), W(KBLOCK,2), TEMP(NUP,2)
  real(kind=wp), pointer :: VTAB(:)

  VTAB => EXS%VTab

  ICP = 1
  start2=0
  i2b=0
  finish2=0

  if (swap) then

  do while (ICP <= NCP)

    I2 = ICOUP(1,ICP)
    start = ICP
    finish = ICP
    do
      if (finish > NCP) exit
      if (ICOUP(1,finish)/=I2) exit
      finish = finish + 1
    end do

    nk = finish-start

    if (finish <= NCP) then
      I2b = ICOUP(1,finish)
      start2 = finish
      finish2 = start2
      do
        if (finish2 > NCP) exit
        if (ICOUP(1,finish2)/=I2b) exit
        finish2 = finish2 + 1
      end do
      nk2 = finish2-start2
    else
      nk2 = -1
    end if

    ! structural GEMM only if perfectly safe AND fits buffers
    if ((nk == nk2) .and. (nk > 0) .and. (nk <= KBLOCK)) then
      do i=1,nk
        I1a(i)=ICOUP(2,start+i-1)
        I1b(i)=ICOUP(2,start2+i-1)
      end do
      if (all(I1a(1:nk)==I1b(1:nk))) then
        do i=1,nk
          Xa(i)=CPQ*VTAB(ICOUP(3,start+i-1))
          Xb(i)=CPQ*VTAB(ICOUP(3,start2+i-1))
          ABLOCK(:,i)=CI(:,I1a(i))
          W(i,1)=Xa(i)
          W(i,2)=Xb(i)
        end do

        call DGEMM_('N','N',NUP,2,nk,One,ABLOCK,NUP,W,KBLOCK,Zero,TEMP,NUP)

        SIGMA(:,I2)  = SIGMA(:,I2)  + TEMP(:,1)
        SIGMA(:,I2b) = SIGMA(:,I2b) + TEMP(:,2)

        ICP = finish2
        cycle
      end if
    end if

    ! fallback blocked GEMV
    do offset = 1, nk, KBLOCK
      blck = min(KBLOCK,nk-offset+1)

      do i=1,blck
        I1a(i)=ICOUP(2,start+offset+i-2)
        Xa(i)=CPQ*VTAB(ICOUP(3,start+offset+i-2))
        ABLOCK(:,i)=CI(:,I1a(i))
      end do

      call DGEMM_('N','N', NUP, 1, blck, 1.0_wp, ABLOCK, NUP, Xa, KBLOCK, 0.0_wp, TEMP, NUP)

      SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
    end do

    ICP = finish

  end do

  else

  do while (ICP <= NCP)

    I2 = ICOUP(2,ICP)
    start = ICP
    finish = ICP
    do
      if (finish > NCP) exit
      if (ICOUP(2,finish)/=I2) exit
      finish = finish + 1
    end do

    nk = finish-start

    if (finish <= NCP) then
      I2b = ICOUP(2,finish)
      start2 = finish
      finish2 = start2
      do
        if (finish2 > NCP) exit
        if (ICOUP(2,finish2)/=I2b) exit
        finish2 = finish2 + 1
      end do
      nk2 = finish2-start2
    else
      nk2 = -1
    end if

    if ((nk == nk2) .and. (nk > 0) .and. (nk <= KBLOCK)) then
      do i=1,nk
        I1a(i)=ICOUP(1,start+i-1)
        I1b(i)=ICOUP(1,start2+i-1)
      end do
      if (all(I1a(1:nk)==I1b(1:nk))) then
        do i=1,nk
          Xa(i)=CPQ*VTAB(ICOUP(3,start+i-1))
          Xb(i)=CPQ*VTAB(ICOUP(3,start2+i-1))
          ABLOCK(:,i)=CI(:,I1a(i))
          W(i,1)=Xa(i)
          W(i,2)=Xb(i)
        end do

        call DGEMM_('N','N',NUP,2,nk,One,ABLOCK,NUP,W,KBLOCK,Zero,TEMP,NUP)

        SIGMA(:,I2)  = SIGMA(:,I2)  + TEMP(:,1)
        SIGMA(:,I2b) = SIGMA(:,I2b) + TEMP(:,2)

        ICP = finish2
        cycle
      end if
    end if

    do offset = 1, nk, KBLOCK
      blck = min(KBLOCK,nk-offset+1)

      do i=1,blck
        I1a(i)=ICOUP(1,start+offset+i-2)
        Xa(i)=CPQ*VTAB(ICOUP(3,start+offset+i-2))
        ABLOCK(:,i)=CI(:,I1a(i))
      end do

      call DGEMM_('N','N',NUP,1,blck,One,ABLOCK,NUP,Xa,KBLOCK,Zero,TEMP,NUP)

      SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
    end do

    ICP = finish

  end do

  end if

  nullify(VTAB)

end subroutine apply_col

  subroutine apply_row(CPQ,NDWN,NUPC,CI,NUPSG,SIGMA,NCP,ICOUP,swap)

    integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP, ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUPC,NDWN)
  real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp), parameter :: KBLOCK=16
    integer(kind=iwp) :: blk, i, ICP, IDWN, j, nblk
    real(kind=wp) :: CI_blk(KBLOCK,KBLOCK), X
    integer(kind=iwp), pointer :: I1LIST(:), I2LIST(:)
    real(kind=wp), pointer :: VTAB(:), XLIST(:)
  logical, parameter :: USE_OPT=.true.

  VTAB => EXS%VTab
  XLIST => EXS%XLIST
  I1LIST => EXS%I1LIST
  I2LIST => EXS%I2LIST

  if (.not. USE_OPT) then

    if (swap) then
      do ICP=1,NCP
        I1list(ICP)=ICOUP(2,ICP)
        I2list(ICP)=ICOUP(1,ICP)
        Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
      end do
    else
      do ICP=1,NCP
        I1list(ICP)=ICOUP(1,ICP)
        I2list(ICP)=ICOUP(2,ICP)
        Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
      end do
    end if


    do ICP=1,NCP
!$OMP SIMD
      do IDWN=1,NDWN
        SIGMA(I2list(ICP),IDWN)=SIGMA(I2list(ICP),IDWN)+Xlist(ICP)*CI(I1list(ICP),IDWN)
      end do
    end do

  else

    ! sorted input improves locality (optional external sort)
    if (swap) then
      do ICP=1,NCP
        I1list(ICP)=ICOUP(2,ICP)
        I2list(ICP)=ICOUP(1,ICP)
        Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
      end do
    else
      do ICP=1,NCP
        I1list(ICP)=ICOUP(1,ICP)
        I2list(ICP)=ICOUP(2,ICP)
        Xlist(ICP)=CPQ*VTAB(ICOUP(3,ICP))
      end do
    end if

    ! tiled processing
    do j=1,NDWN,KBLOCK
      nblk = min(KBLOCK, NDWN-j+1)

      do i=1,NCP,KBLOCK
        blk = min(KBLOCK, NCP-i+1)

        ! pack tile (transpose-like)
        do ICP=1,blk
          do IDWN=1,nblk
            CI_blk(IDWN,ICP)=CI(I1list(i+ICP-1), j+IDWN-1)
          end do
        end do

        ! compute
        do ICP=1,blk
          X = Xlist(i+ICP-1)
!$OMP SIMD
          do IDWN=1,nblk
            SIGMA(I2list(i+ICP-1), j+IDWN-1) =               SIGMA(I2list(i+ICP-1), j+IDWN-1) + X*CI_blk(IDWN,ICP)
          end do
        end do

      end do
    end do

  end if

    nullify(VTAB)

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

! declarations to facilitate the reuse of sigma vectors if possible
integer(kind=iwp), save:: i_save_p=0, i_save_q=0
integer(kind=iwp), save:: i_save_p_sym=-1, i_save_q_sym=-1
logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2, Test3
real(kind=wp) :: CI_ID=Zero
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
    Reuse_Sigma = .false.
    Test3=MERGE((i_save_q == IQ) .and. (i_save_p_sym == ISYP), &
                (i_save_p == IP) .and. (i_save_q_sym == ISYQ), Swap)
    if (EXS%Reuse_SGTMP .and. Test3) Reuse_Sigma = CI_ID == sum(CI(1:nCSFs))+real(nCSFs,kind=wp)

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
    If (SWAP) Then
      i_save_q=IQ
      i_save_p_sym=ISYP
      i_save_p=0
      i_save_q_sym=-1
    Else
      i_save_p=IP
      i_save_q_sym=ISYQ
      i_save_q=0
      i_save_p_sym=-1
    End If
    if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) CI_ID = sum(CI(1:nCSFs))+real(nCSFs,kind=wp)
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

! This is temporarily for the MCLR code
subroutine SG_ReOrd_orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
!***********************************************************************
!                                                                      *
!     Rearrange CI-vectors                                             *
!     iMode=0 --> from SGA to split graph GUGA order                   *
!     iMode=1 --> from split graph GUGA to SGA order                   *
!                                                                      *
!     calling arguments:                                               *
!     iRefSm  : integer                                                *
!               state symmetry                                         *
!     iMode   : integer                                                *
!               switch selecting reordering mode (see above)           *
!     nSm     : array of integer                                       *
!               symmetry per active orbital                            *
!     CIold   : array of real                                          *
!               incoming CI vector                                     *
!     CInew   : array of real                                          *
!               outgoing CI vector                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and J. Olsen                                      *
!     University of Lund, Sweden, 1990                                 *
!***********************************************************************

use spinfo, only: MINOP, NCNFTP, NCSFTP, NTYP
use Lucia_data, only: CFTP, CONF_Occ
use Molcas, only: MxAct
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iState, IREFSM, IMODE, nConf
real(kind=wp), intent(in) :: CIOLD(nConf)
real(kind=wp), intent(out) :: CINEW(nConf)
integer(kind=iwp) :: IC, ICL, ICNBS, ICNBS0, ICSBAS, ICSFJP, IIBCL, IIBOP, IICSF, IOPEN, IP, IPBAS, ISG, ITYP, &
                     IWALK(mxAct), JOCC, KCNF(MxAct), KOCC, KORB
real(kind=wp) :: Fact
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

If (.NOT.Allocated(CIS(iState)%ICASE)) Call MkCOT(SGS(istate),CIS(istate))
If (.NOT.Allocated(EXS(iState)%USGN)) Call MkSgNum(IREFSM,SGS(istate),CIS(istate),EXS(istate))


ICSFJP = 0
ICNBS0 = 0 ! dummy initialize
IPBAS = 0 ! dummy initialize
! LOOP OVER CONFIGURATIONS TYPES
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1
  ICL = (SGS(iState)%nActEl-IOPEN)/2
  ! BASE ADDRESS FOR CONFIGURATION OF THIS TYPE
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNFTP(ITYP-1,IREFSM)*(SGS(iState)%nActEl+IOPEN-1)/2
  end if
  ! BASE ADDRESS FOR PROTOTYPE SPIN COUPLINGS
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NCSFTP(ITYP-1)*(IOPEN-1)
  end if

  ! LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
  ! SPIN COUPLINGS

  do IC=1,NCNFTP(ITYP,IREFSM)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    do IICSF=1,NCSFTP(ITYP)
      ICSFJP = ICSFJP+1
      ICSBAS = IPBAS+(IICSF-1)*IOPEN
      KCNF(:)=0
      ! Obtain configuration in standard RASSCF form
      IIBOP = 1
      IIBCL = 1
      JOCC = ICL+IOPEN
      do KOCC=0,JOCC-1
        KORB = Conf_Occ(IREFSM)%A(ICNBS+KOCC)
        if (KORB < 0) then
          ! Doubly occupied orbital
          KCNF(IIBCL) = abs(KORB)
          IIBCL = IIBCL+1
        else
          ! Singly occupied orbital
          KCNF(ICL+IIBOP) = KORB
          IIBOP = IIBOP+1
        end if
      end do

      ! COMPUTE STEP VECTOR
      call STEPVEC(KCNF(1:ICL),KCNF(ICL+1),ICL,IOPEN,CFTP(ICSBAS),SGS(iState)%nLev,IWALK)

      ! GET SPLIT GRAPH ORDERING NUMBER
      ISG = SG_NUM(SGS(iState),EXS(istate),IWALK)
      ! GET PHASE PHASE FACTOR
      IP = SG_PHASE(SGS(istate),IWALK)
      Fact = merge(-One,One,IP < 0)
      if (IMODE == 0) then
        CINEW(ISG) = Fact * CIOLD(ICSFJP)
      else
        CINEW(ICSFJP) = Fact * CIOLD(ISG)
      end if
    end do
  end do
end do

#ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CIOLD(I),I=1,min(200,ICSFJP))
  write(u6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CINEW(I),I=1,min(200,ICSFJP))
  write(u6,*)
#endif

end subroutine SG_Reord_Orig

subroutine SG_ReOrd(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)

use Molcas, only: MxAct
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none

integer(kind=iwp), intent(in)  :: iState
integer(kind=iwp), intent(in)  :: IREFSM
integer(kind=iwp), intent(in)  :: IMODE
integer(kind=iwp), intent(in)  :: nConf

real(kind=wp), intent(in)  :: CIOLD(nConf)
real(kind=wp), intent(out) :: CINEW(nConf)

integer(kind=iwp) :: IC
integer(kind=iwp) :: ICL
integer(kind=iwp) :: ICSFJP
integer(kind=iwp) :: IICSF
integer(kind=iwp) :: IOPEN
integer(kind=iwp) :: IP
integer(kind=iwp) :: ISG
integer(kind=iwp) :: MinOpen, nOC
integer(kind=iwp) :: MaxOpen

integer(kind=iwp) :: IWALK(MxAct)
integer(kind=iwp) :: KCNF(MxAct)

real(kind=wp) :: Fact

#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif
If (SGS(iState)%nRsPrt /= 1) Then
   Call SG_ReOrd_orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
   Return
endif

if (.not. allocated(CIS(iState)%ICASE)) call MkCOT(SGS(iState),CIS(iState))
call MkSgNum(IREFSM,SGS(iState),CIS(iState),EXS(iState))

ICSFJP = 0

MinOpen = SGS(iState)%iSpin - 1
MaxOpen = min(SGS(iState)%nActEl, &
              2*SGS(iState)%nLev-SGS(iState)%nActEl)
if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen - 1

! loop over open-shell sectors. That is, run over all orbital configurations
! that have the same number of open-shells
do iOpen = MinOpen, MaxOpen, 2

   ICL = (SGS(iState)%nActEl-iOpen)/2

   nOC = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nOCSym(IREFSM)

!  Run over all prototype configurations of the sector
   do IC = 1, nOC

      call GetOccupation_SGUGA(iState,iOpen,SGS(iState)%iSpin,IREFSM,IC,KCNF)

!     Loop over all orbital combinations for this particular sector.
      do IICSF = 1, SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nCSF

         ICSFJP = ICSFJP + 1

         call BuildWalk(iState,ICL,iOpen,IICSF,KCNF,IWALK)

#ifdef _SANTITY_CHECK_
         if (any(IWALK(1:SGS(iState)%nLev) < 0) .or. &
             any(IWALK(1:SGS(iState)%nLev) > 3)) then

            write(u6,*) 'INVALID IWALK'
            write(u6,'(20I4)') IWALK(1:SGS(iState)%nLev)
            call Abend()
         endif
#endif

         ISG = SG_NUM(SGS(iState),EXS(iState),IWALK)
         IP = SG_PHASE(SGS(iState),IWALK)
         Fact = One
         if (IP < 0) Fact = -One
         if (IMODE == 0) then
            CINEW(ISG) = Fact*CIOLD(ICSFJP)
         else
            CINEW(ICSFJP) = Fact*CIOLD(ISG)
         endif
      enddo
   enddo
enddo

#ifdef _DEBUGPRINT_

write(u6,*)
write(u6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CIOLD(I),I=1,min(200,ICSFJP))

write(u6,*)
write(u6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CINEW(I),I=1,min(200,ICSFJP))
write(u6,*)

#endif

end subroutine SG_ReOrd

!> Step vector IWALK of the CSF ICSF of the Lucia (SGA) CSF basis of
!> symmetry IREFSM and the phase IP relating it to the GUGA CSF, as in
!> SG_ReOrd: CINEW(SG_NUM(IWALK)) = IP*CIOLD(ICSF) for IMODE = 0.
!> Only for a single RAS partition (nRsPrt = 1, the case of SG_ReOrd proper).
subroutine SG_Lucia_Walk(iState,IREFSM,ICSF,IWALK,IP)

use Molcas, only: MxAct

integer(kind=iwp), intent(in) :: iState, IREFSM, ICSF
integer(kind=iwp), intent(out) :: IWALK(MxAct), IP
integer(kind=iwp) :: iOpen, ICL, MinOpen, MaxOpen, nOC, nCSF, iOff, IC, IICSF
integer(kind=iwp) :: KCNF(MxAct)

if (SGS(iState)%nRsPrt /= 1) then
  write(u6,*) 'SG_Lucia_Walk: only for nRsPrt = 1'
  call Abend()
end if

MinOpen = SGS(iState)%iSpin-1
MaxOpen = min(SGS(iState)%nActEl,2*SGS(iState)%nLev-SGS(iState)%nActEl)
if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen-1

iOff = 0
do iOpen=MinOpen,MaxOpen,2
  nOC = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nOCSym(IREFSM)
  nCSF = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nCSF
  if (ICSF <= iOff+nOC*nCSF) then
    ICL = (SGS(iState)%nActEl-iOpen)/2
    IC = (ICSF-iOff-1)/nCSF+1
    IICSF = ICSF-iOff-(IC-1)*nCSF
    call GetOccupation_SGUGA(iState,iOpen,SGS(iState)%iSpin,IREFSM,IC,KCNF)
    call BuildWalk(iState,ICL,iOpen,IICSF,KCNF,IWALK)
    IP = SG_PHASE(SGS(iState),IWALK)
    return
  end if
  iOff = iOff+nOC*nCSF
end do

write(u6,*) 'SG_Lucia_Walk: CSF index out of range:',ICSF
call Abend()

end subroutine SG_Lucia_Walk

subroutine GetOccupation_SGUGA(iState,iOpen,mult,iSym,iOC,KCNF)
   use Molcas, only: MxAct
   integer(kind=iwp), intent(in)  :: iState
   integer(kind=iwp), intent(in)  :: iOpen
   integer(kind=iwp), intent(in)  :: mult
   integer(kind=iwp), intent(in)  :: iSym
   integer(kind=iwp), intent(in)  :: iOC
   integer(kind=iwp), intent(out) :: KCNF(MxAct)

   integer(kind=iwp) :: iDOC
   integer(kind=iwp) :: iSOC
   integer(kind=iwp) :: iCL
   integer(kind=iwp) :: i

   KCNF(:) = 0

   iCL = SGS(iState)%CDCInfo(iOpen,mult)%nClosed

   iDOC = SGS(iState)%CDCInfo(iOpen,mult)%OCDOC(iOC,iSym)
   iSOC = SGS(iState)%CDCInfo(iOpen,mult)%OCSOC(iOC,iSym)

   if (iCL > 0) then
      do i=1,iCL
         KCNF(i) = SGS(iState)%CDCInfo(iOpen,mult)%DOC(i,iDOC)
      end do
   end if

   if (iOpen > 0) then
      do i=1,iOpen
         KCNF(iCL+i) = SGS(iState)%CDCInfo(iOpen,mult)%SOC(i,iSOC)
      end do
   end if

end subroutine GetOccupation_SGUGA

subroutine BuildWalk(iState,ICL,IOPEN,IICSF,KCNF,IWALK)
   use Molcas, only: MxAct
   integer(kind=iwp), intent(in)  :: istate
   integer(kind=iwp), intent(in)  :: ICL
   integer(kind=iwp), intent(in)  :: IOPEN
   integer(kind=iwp), intent(in)  :: IICSF
   integer(kind=iwp), intent(in)  :: KCNF(MxAct)
   integer(kind=iwp), intent(out) :: IWALK(MxAct)

   integer(kind=iwp) :: i
   integer(kind=iwp) :: Twos(0:IOPEN)

   if (iOpen == 0) then
      IWALK(:) = 0
      do i=1,ICL
         IWALK(KCNF(i)) = 3
      end do
      return
   endif

   Twos(0:IOPEN) = SGS(iState)%CDCInfo(IOPEN,SGS(iState)%iSpin)%Twos(0:IOPEN,IICSF)

   IWALK(:) = 0

   ! Doubly occupied orbitals
   do i=1,ICL
      IWALK(KCNF(i)) = 3
   end do

   ! Singly occupied orbitals
   do i=1,IOPEN
      if (Twos(i) > Twos(i-1)) then
         IWALK(KCNF(ICL+i)) = 1
      else
         IWALK(KCNF(ICL+i)) = 2
      end if
   end do

End subroutine BuildWalk

!===========================================================
! Build walk from Lucia occupation class and CFTP data
!===========================================================

#ifdef _NOT_USED_ANY_LONGER_
subroutine BuildWalk_Lucia(ICL,IOPEN,ICSBAS,KCNF,IWALK)

   use Lucia_data, only: CFTP
   use Molcas, only: MxAct
   integer(kind=iwp), intent(in)  :: ICL
   integer(kind=iwp), intent(in)  :: IOPEN
   integer(kind=iwp), intent(in)  :: ICSBAS

   integer(kind=iwp), intent(in)  :: KCNF(MxAct)

   integer(kind=iwp), intent(out) :: IWALK(MxAct)

   call STEPVEC( &
      KCNF(1:ICL), &
      KCNF(ICL+1:ICL+IOPEN), &
      ICL, &
      IOPEN, &
      CFTP(ICSBAS:ICSBAS+IOPEN-1), &
      SGS(iState)%nLev, &
      IWALK )

end subroutine BuildWalk_Lucia
#endif


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

subroutine SG_Test_EpqErs_Ref(iState,LSM,Epq_Only)

use Constants,   only: Zero, One
use Definitions, only: wp, iwp, u6
use stdalloc,    only: mma_allocate, mma_deallocate
use, intrinsic :: iso_fortran_env, only: int64

implicit none

! If Epq_Only is present and true, only <I|Epq|J> is verified and the
! <I|EpqErs|J> tests are skipped.
logical(kind=iwp), intent(in), optional :: Epq_Only

integer(kind=iwp), intent(in) :: iState
integer(kind=iwp), intent(in) :: LSM

!
! Operator-product topologies.
!
integer(kind=iwp), parameter :: iUU = 1
integer(kind=iwp), parameter :: iLL = 2
integer(kind=iwp), parameter :: iUL = 3
integer(kind=iwp), parameter :: iLU = 4
integer(kind=iwp), parameter :: iUO = 5
integer(kind=iwp), parameter :: iOU = 6
integer(kind=iwp), parameter :: iLO = 7
integer(kind=iwp), parameter :: iOL = 8
integer(kind=iwp), parameter :: iOO = 9

integer(kind=iwp), parameter :: nTopology = 9

character(len=2), parameter :: TopologyName(nTopology) = &
   [ character(len=2) :: 'UU', 'LL', 'UL', 'LU', 'UO', 'OU', 'LO', 'OL', 'OO' ]

!
! Diagonal-status classes.
!
!     CC: closed   x closed
!     DC: diagonal x closed
!     CD: closed   x diagonal
!     DD: diagonal x diagonal
!
integer(kind=iwp), parameter :: iCC = 1
integer(kind=iwp), parameter :: iDC = 2
integer(kind=iwp), parameter :: iCD = 3
integer(kind=iwp), parameter :: iDD = 4

integer(kind=iwp), parameter :: nDiagonalCase = 4

character(len=2), parameter :: DiagonalName(nDiagonalCase) = &
   [ character(len=2) :: 'CC', 'DC', 'CD', 'DD' ]

!
! Orientation classes for closed x closed products.
!
!     CC: canonical x canonical
!     RC: reverse   x canonical
!     CR: canonical x reverse
!     RR: reverse   x reverse
!
integer(kind=iwp), parameter :: iCanCan = 1
integer(kind=iwp), parameter :: iRevCan = 2
integer(kind=iwp), parameter :: iCanRev = 3
integer(kind=iwp), parameter :: iRevRev = 4

integer(kind=iwp), parameter :: nOrientationCase = 4

character(len=2), parameter :: &
   OrientationName(nOrientationCase) = &
   [ character(len=2) :: 'CC', 'RC', 'CR', 'RR' ]

integer(kind=iwp) :: iTopology
integer(kind=iwp) :: iDiagonalCase
integer(kind=iwp) :: iOrientationCase

integer(kind=iwp) :: p
integer(kind=iwp) :: q
integer(kind=iwp) :: r
integer(kind=iwp) :: s

integer(kind=iwp) :: pMin
integer(kind=iwp) :: pMax
integer(kind=iwp) :: qMin
integer(kind=iwp) :: qMax
integer(kind=iwp) :: rMin
integer(kind=iwp) :: rMax
integer(kind=iwp) :: sMin
integer(kind=iwp) :: sMax

integer(kind=iwp) :: I, J, K

integer(kind=iwp) :: nCSF

!
! Counters indexed by operator-product topology.
!
integer(kind=iwp) :: nOperatorPairs(nTopology)
integer(kind=iwp) :: nMatrixElements(nTopology)
integer(kind=iwp) :: nDiscrepancies(nTopology)

!
! Diagonal-status counters:
!
!     first index  = topology
!     second index = CC, DC, CD, or DD
!
integer(kind=iwp) :: nPairDiagonal( &
                           nTopology,nDiagonalCase )

integer(kind=iwp) :: nErrDiagonal( &
                           nTopology,nDiagonalCase )

!
! Orientation counters for closed x closed products:
!
!     first index  = topology
!     second index = CC, RC, CR, or RR
!
integer(kind=iwp) :: nPairOrientation( &
                           nTopology,nOrientationCase )

integer(kind=iwp) :: nErrOrientation( &
                           nTopology,nOrientationCase )

real(kind=wp) :: XRef
real(kind=wp) :: XTest
real(kind=wp) :: Difference

real(kind=wp), parameter :: Threshold = 1.0e-12_wp

logical :: FirstDiagonal
logical :: SecondDiagonal

logical :: FirstCanonical
logical :: SecondCanonical

integer(kind=iwp) :: nTest, nFail
real(kind=wp)     :: MaxDiff
Integer(kind=iwp) :: iSYCI
!
! Quick zero test: orbital occupations of the CSFs of symmetry LSM,
! Occ(t,I) = occupation of level t in CSF I, and the number of matrix
! elements per topology that were screened as exactly zero.
!
integer(kind=iwp), allocatable :: Occ(:,:)
integer(kind=iwp) :: nScreened(nTopology)
!
! Sampling of the (I,J) pairs. For a given operator (pair), all pairs are
! tested if there are at most MaxPairs of them. Otherwise random pairs are
! drawn (at most MaxDraw, and at most twice the number of pairs) and
! accepted until QuotaNZ pairs that pass the
! quick zero test (compared with the reference) and QuotaZ pairs that are
! zero by occupation (fast routine must return exactly zero) are taken.
! The seed is taken from the clock, so each run tests a different subset,
! unless MOLCAS_RANDOM_SEED is set. The seed is printed before the first
! sampled pair, so a failed run can be repeated with
! >export MOLCAS_RANDOM_SEED=<seed>.
!
integer(kind=iwp), parameter :: MaxPairs = 40000, QuotaNZ = 10000, &
                                QuotaZ = 10000, MaxDraw = 2000000
integer(kind=int64) :: iRandState, iRandSeed
logical :: SeedFromEnv, SeedReported
integer(kind=iwp) :: nPairsIJ, nDrawMax, iDraw, nTakeZ, nTakeNZ
integer(kind=iwp) :: nSampledOps(nTopology)
logical :: Sampled
logical :: Screened
#ifdef _CHECK_QUICK_ZERO_
real(kind=wp) :: XChk
#endif

#ifdef _DEBUGPRINT_

Do ISYCI=1, SGS(istate)%nSym
  Write (u6,*) 'ISYCI=',ISYCI
  Call Print_CSF_Walks(SGS(iState),CIS(iState),ISYCI)
End Do
#endif

! State and symmetry of the module-level matrix-element routines
ME_iState = iState
ME_LSM = LSM
call Reset_Decode_ME()
if (.not. allocated(EXS(iState)%IOrd)) call Mk_Coup_Index(EXS(iState))
! Verify that MatEl_Epq is correct

call Init_Rand_Seed()
Do ISYCI=1, SGS(istate)%nSym
   call Check_MatEl_Epq(iState,ISYCI,1.0e-12_wp,nTest,nFail,MaxDiff)
   if (nFail /= 0) call Abend()
End Do

if (present(Epq_Only)) then
   if (Epq_Only) return
end if

nCSF = CIS(iState)%NCSF(LSM)

!
! Orbital occupations of all CSFs of symmetry LSM, for the quick zero test,
! unpacked from the step vectors (O(nCSF*nLev)).
!
call mma_allocate(Occ,SGS(iState)%nLev,nCSF,Label='Occ')
do I=1,nCSF
   call CSF_Occupation(LSM,I,Occ(:,I))
end do
nScreened = 0
nSampledOps = 0

nOperatorPairs  = 0
nMatrixElements = 0
nDiscrepancies  = 0

nPairDiagonal = 0
nErrDiagonal  = 0

nPairOrientation = 0
nErrOrientation  = 0

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' GENERAL TWO-OPERATOR TEST'
! List the topologies from TopologyName, so the header stays in step with
! the topologies actually tested.
write(u6,'(2X,*(A,:,", "))') (TopologyName(iTopology),iTopology=1,nTopology)
write(u6,*) ' INCLUDING DIAGONAL OPERATORS'
! The products E_pq E_rs are restricted to Gamma(p)Gamma(q) = Gamma(r)Gamma(s)
! (totally symmetric, as in H); each operator may change the symmetry.
write(u6,*) ' TOTALLY SYMMETRIC PRODUCTS (SINGLE OPERATORS OF ANY SYMMETRY)'
write(u6,*) '========================================='
write(u6,'(A,I0)') 'Number of graph levels: ', &
                    SGS(iState)%nLev
write(u6,'(A,I0)') 'Split-graph midlevel:   ', &
                    SGS(iState)%MidLev
write(u6,'(A,I0)') 'Number of CSFs:         ', &
                    nCSF
write(u6,*)

!
! Test all nine operator-product topologies:
!
!     UU: E_pq(upper) E_rs(upper)
!     LL: E_pq(lower) E_rs(lower)
!     UL: E_pq(upper) E_rs(lower)
!     LU: E_pq(lower) E_rs(upper)
!     UO: E_pq(upper) E_rs(open, i.e. crossing the midlevel)
!     OU: E_pq(open) E_rs(upper)
!     LO: E_pq(lower) E_rs(open)
!     OL: E_pq(open) E_rs(lower)
!     OO: E_pq(open) E_rs(open)
!
do iTopology=1,nTopology

   ! Empty ranges by default: every valid topology sets them below and the
   ! default case aborts, but some compilers cannot see that (Abend does
   ! not return) and warn that they may be used uninitialized.
   pMin = 1; pMax = 0
   qMin = 1; qMax = 0
   rMin = 1; rMax = 0
   sMin = 1; sMax = 0

   select case (iTopology)

   case (iUU)
      pMin=SGS(iState)%MidLev+1; pMax=SGS(iState)%nLev
      qMin=pMin; qMax=pMax
      rMin=SGS(iState)%MidLev+1; rMax=SGS(iState)%nLev
      sMin=rMin; sMax=rMax
   case (iLL)
      pMin=1; pMax=SGS(iState)%MidLev
      qMin=pMin; qMax=pMax
      rMin=1; rMax=SGS(iState)%MidLev
      sMin=rMin; sMax=rMax
   case (iUL)
      pMin=SGS(iState)%MidLev+1; pMax=SGS(iState)%nLev
      qMin=pMin; qMax=pMax
      rMin=1; rMax=SGS(iState)%MidLev
      sMin=rMin; sMax=rMax
   case (iLU)
      pMin=1; pMax=SGS(iState)%MidLev
      qMin=pMin; qMax=pMax
      rMin=SGS(iState)%MidLev+1; rMax=SGS(iState)%nLev
      sMin=rMin; sMax=rMax
   case (iUO)
      pMin=SGS(iState)%MidLev+1; pMax=SGS(iState)%nLev
      qMin=pMin; qMax=pMax
      rMin=1; rMax=SGS(iState)%nLev
      sMin=1; sMax=SGS(iState)%nLev
   case (iOU)
      ! First operator open (filtered below), second operator upper.
      pMin=1; pMax=SGS(iState)%nLev
      qMin=1; qMax=SGS(iState)%nLev
      rMin=SGS(iState)%MidLev+1; rMax=SGS(iState)%nLev
      sMin=rMin; sMax=rMax
   case (iLO)
      ! First operator lower, second operator open (filtered below).
      pMin=1; pMax=SGS(iState)%MidLev
      qMin=1; qMax=SGS(iState)%MidLev
      rMin=1; rMax=SGS(iState)%nLev
      sMin=1; sMax=SGS(iState)%nLev
   case (iOL)
      ! First operator open (filtered below), second operator lower.
      pMin=1; pMax=SGS(iState)%nLev
      qMin=1; qMax=SGS(iState)%nLev
      rMin=1; rMax=SGS(iState)%MidLev
      sMin=1; sMax=SGS(iState)%MidLev
   case (iOO)
      ! Both operators open (filtered below).
      pMin=1; pMax=SGS(iState)%nLev
      qMin=1; qMax=SGS(iState)%nLev
      rMin=1; rMax=SGS(iState)%nLev
      sMin=1; sMax=SGS(iState)%nLev
!     Restrict to only canonical case r > s
!     rMin=SGS(iState)%MidLev+1; rMax=SGS(iState)%nLev
!     sMin=1; sMax=SGS(iState)%MidLev

   case default

      write(u6,'(A,I0)') &
         'Illegal operator-product topology: ',iTopology
      call Abend()

   end select

   write(u6,'(A,A)') &
      'Testing operator-product topology: ', &
      TopologyName(iTopology)

   !
   ! Enumerate every ordered operator product E_pq E_rs
   ! belonging to the current topology.
   !
   do p=pMin,pMax
      do q=qMin,qMax

         if ((iTopology==iOU) .or. (iTopology==iOL) .or. (iTopology==iOO)) then
            ! Enforce that the first operator is open
            if (p==q) cycle
            if ((p>SGS(iState)%MidLev) .eqv. &
                (q>SGS(iState)%MidLev)) cycle
         end if

         FirstDiagonal  = p == q
         FirstCanonical = p > q

         !
         ! Diagonal operators are totally symmetric; the others may
         ! change the symmetry.
         !
         ! (Each operator may change the symmetry; only the product has
         ! to be totally symmetric, see the second operator.)

         do r=rMin,rMax
            do s=sMin,sMax

               if ((iTopology==iUO) .or. (iTopology==iLO) .or. (iTopology==iOO)) then
                  ! Enforce that the second operator is open
                  if (r==s) cycle
                  if ((r>SGS(iState)%MidLev) .eqv. &
                      (s>SGS(iState)%MidLev)) cycle
               end if

               SecondDiagonal  = r == s
               SecondCanonical = r > s

               ! Totally symmetric product: Gamma(p)Gamma(q) = Gamma(r)Gamma(s),
               ! as in H. The single operators may change the symmetry.
               if (Mul(SGS(iState)%ISM(r),SGS(iState)%ISM(s)) /= &
                   Mul(SGS(iState)%ISM(p),SGS(iState)%ISM(q))) cycle

               nOperatorPairs(iTopology) = &
                  nOperatorPairs(iTopology)+1

               !
               ! Classify according to diagonal status.
               !
               if (.not. FirstDiagonal .and. &
                   .not. SecondDiagonal) then

                  iDiagonalCase = iCC

                  if (FirstCanonical .and. &
                      SecondCanonical) then

                     iOrientationCase = iCanCan

                  else if (.not. FirstCanonical .and. &
                            SecondCanonical) then

                     iOrientationCase = iRevCan

                  else if (FirstCanonical .and. &
                           .not. SecondCanonical) then

                     iOrientationCase = iCanRev

                  else

                     iOrientationCase = iRevRev

                  end if

                  nPairOrientation( &
                     iTopology,iOrientationCase ) = &
                     nPairOrientation( &
                     iTopology,iOrientationCase )+1

               else if (FirstDiagonal .and. &
                       .not. SecondDiagonal) then

                  iDiagonalCase    = iDC
                  iOrientationCase = 0

               else if (.not. FirstDiagonal .and. &
                        SecondDiagonal) then

                  iDiagonalCase    = iCD
                  iOrientationCase = 0

               else

                  iDiagonalCase    = iDD
                  iOrientationCase = 0

               end if

               nPairDiagonal( &
                  iTopology,iDiagonalCase ) = &
                  nPairDiagonal( &
                  iTopology,iDiagonalCase )+1

               !
               ! Compare the CSF matrices: all (I,J) pairs, or a random
               ! sample of them if there are more than MaxPairs.
               !
               nPairsIJ = nCSF*nCSF
               Sampled  = nPairsIJ > MaxPairs
               ! Cap the draws relative to the block size, so that sampling a block
               ! just above MaxPairs never costs more than two full loops.
               nDrawMax = merge(min(MaxDraw,2*nPairsIJ),nPairsIJ,Sampled)
               nTakeZ   = 0
               nTakeNZ  = 0
               if (Sampled) nSampledOps(iTopology) = nSampledOps(iTopology)+1
               do iDraw=1,nDrawMax
                     if (Sampled) then
                        if (nTakeZ >= QuotaZ .and. nTakeNZ >= QuotaNZ) exit
                        I = Rand_Index(nCSF)
                        J = Rand_Index(nCSF)
                     else
                        I = (iDraw-1)/nCSF+1
                        J = mod(iDraw-1,nCSF)+1
                     end if

                     ! Skip the expensive reference when the occupations
                     ! already show that <I|EpqErs|J> is zero. The fast
                     ! routine is still called and must then return 0.
                     Screened = Quick_Zero(Occ(:,I),Occ(:,J),p,q,r,s)
                     if (Sampled) then
                        if (Screened) then
                           if (nTakeZ >= QuotaZ) cycle
                           nTakeZ = nTakeZ+1
                        else
                           if (nTakeNZ >= QuotaNZ) cycle
                           nTakeNZ = nTakeNZ+1
                        end if
                     end if
                     if (Screened) then
                        XRef = Zero
                        nScreened(iTopology) = nScreened(iTopology)+1
#ifdef _CHECK_QUICK_ZERO_
                        ! Verify the quick zero test itself (as expensive as
                        ! no screening at all; define it only to test it).
                        XChk = MatEl_EpqErs_Ref(LSM,I,J,p,q,r,s)
                        if (abs(XChk) > Threshold) then
                           write(u6,'(A,6I6,ES22.14)') &
                              'Quick_Zero failed: I,J,p,q,r,s,Ref = ', &
                              I,J,p,q,r,s,XChk
                           call Abend()
                        end if
#endif
                     else
                        XRef = MatEl_EpqErs_Ref( &
                                   LSM, &
                                   I,J, &
                                   p,q, &
                                   r,s )
                     end if

                     ! Set in every valid case below; the default aborts.
                     XTest = Zero
                     select case (iTopology)

                     case (iUU)

                        XTest = MatEl_EpqErs_UU( &
                                   LSM, &
                                   I,J, &
                                   p,q, &
                                   r,s )

                     case (iLL)

                        XTest = MatEl_EpqErs_LL( &
                                   LSM, &
                                   I,J, &
                                   p,q, &
                                   r,s )

                     case (iUL)

                        XTest = MatEl_EpqErs_UL( &
                                   LSM, &
                                   I,J, &
                                   p,q, &
                                   r,s )

                     case (iLU)

                        XTest = MatEl_EpqErs_LU( &
                                   LSM, &
                                   I,J, &
                                   p,q, &
                                   r,s )

                     case (iUO)

                        If (I==2 .and. J==10 .and. &
                            p==3 .and. q==4 .and. r==1 .and. s==3) Then
                           Write (u6,*) 'Test <I|Epq|K>'
                           Do K = 1, nCSF
                              If (K/=3) Cycle
                              XTest = MatEl_Epq(SGS(istate),CIS(istate),EXS(istate), &
                                                p,q,LSM,I,K,.True.)
                              Write (u6,*) 'I,p,q,K,<I|Epq|K>=',I,p,q,K,XTest
                           End Do
                           Write (u6,*) 'Test <K|Ers|J>'
                           Do K = 1, nCSF
                              If (K/=3) Cycle
                              XTest = MatEl_Epq(SGS(istate),CIS(istate),EXS(istate), &
                                                r,s,LSM,K,J,.True.)
                              Write (u6,*) 'K,r,s,J,<K|Ers|J>=',K,r,s,J,XTest
                           End Do
                        End If

                        XTest = MatEl_EpqErs_UO( LSM,I,J,p,q,r,s )

                     case (iOU)

                        XTest = MatEl_EpqErs_OU( LSM,I,J,p,q,r,s )

                     case (iLO)

                        XTest = MatEl_EpqErs_LO( LSM,I,J,p,q,r,s )

                     case (iOL)

                        XTest = MatEl_EpqErs_OL( LSM,I,J,p,q,r,s )

                     case (iOO)

                        XTest = MatEl_EpqErs_OO( LSM,I,J,p,q,r,s )

                     case default

                        write(u6,'(A,I0)') &
                           'Illegal topology in comparison: ', &
                           iTopology
                        call Abend()

                     end select

                     nMatrixElements(iTopology) = &
                        nMatrixElements(iTopology)+1

                     Difference = XRef-XTest

                     ! A screened element is exactly zero, so the fast
                     ! routine must return exactly zero, not just a small
                     ! number.
                     if (Screened) then
                        if (XTest == Zero) cycle
                     else
                        if (abs(Difference) <= Threshold) cycle
                     end if

                     nDiscrepancies(iTopology) = &
                        nDiscrepancies(iTopology)+1

                     nErrDiagonal( &
                        iTopology,iDiagonalCase ) = &
                        nErrDiagonal( &
                        iTopology,iDiagonalCase )+1

                     if (iDiagonalCase == iCC) then

                        nErrOrientation( &
                           iTopology,iOrientationCase ) = &
                           nErrOrientation( &
                           iTopology,iOrientationCase )+1

                     end if

                     !
                     ! Columns:
                     !
                     !   topology
                     !   diagonal-status class
                     !   orientation class
                     !   p q r s
                     !   I J
                     !   XRef XTest Difference
                     !
                     if (iDiagonalCase == iCC) then

                        write(u6, &
                           '(A,1X,A,1X,A,1X,4I4,' // &
                           '1X,2I7,3(1X,ES22.14))') &
                           TopologyName(iTopology), &
                           DiagonalName(iDiagonalCase), &
                           OrientationName( &
                              iOrientationCase), &
                           p,q,r,s, &
                           I,J, &
                           XRef,XTest,Difference

                     else

                        write(u6, &
                           '(A,1X,A,1X,A,1X,4I4,' // &
                           '1X,2I7,3(1X,ES22.14))') &
                           TopologyName(iTopology), &
                           DiagonalName(iDiagonalCase), &
                           '  ', &
                           p,q,r,s, &
                           I,J, &
                           XRef,XTest,Difference

                     end if

               end do

            end do
         end do

      end do
   end do

   write(u6,*)

end do

!
! Print separate summaries for UU, LL, UL, and LU.
!
write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' TWO-OPERATOR TEST SUMMARY'
write(u6,*) '========================================='

do iTopology=1,nTopology

   write(u6,*)
   write(u6,'(A,A)') &
      'Topology:               ', &
      TopologyName(iTopology)

   write(u6,'(A,I0)') &
      'Operator pairs:         ', &
      nOperatorPairs(iTopology)

   write(u6,'(A,I0)') &
      'Matrix elements:        ', &
      nMatrixElements(iTopology)

   write(u6,'(A,I0)') &
      'Screened as zero:       ', &
      nScreened(iTopology)

   write(u6,'(A,I0)') &
      'Operator pairs sampled: ', &
      nSampledOps(iTopology)

   write(u6,'(A,I0)') &
      'Discrepancies:          ', &
      nDiscrepancies(iTopology)

   write(u6,*)

   write(u6,'(A,I0,A,I0)') &
      'Closed x closed:        pairs = ', &
      nPairDiagonal(iTopology,iCC), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iCC)

   write(u6,'(A,I0,A,I0)') &
      'Diagonal x closed:      pairs = ', &
      nPairDiagonal(iTopology,iDC), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iDC)

   write(u6,'(A,I0,A,I0)') &
      'Closed x diagonal:      pairs = ', &
      nPairDiagonal(iTopology,iCD), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iCD)

   write(u6,'(A,I0,A,I0)') &
      'Diagonal x diagonal:    pairs = ', &
      nPairDiagonal(iTopology,iDD), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iDD)

   write(u6,*)

   write(u6,'(A,I0,A,I0)') &
      'Canonical x canonical:  pairs = ', &
      nPairOrientation(iTopology,iCanCan), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iCanCan)

   write(u6,'(A,I0,A,I0)') &
      'Reverse x canonical:    pairs = ', &
      nPairOrientation(iTopology,iRevCan), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iRevCan)

   write(u6,'(A,I0,A,I0)') &
      'Canonical x reverse:    pairs = ', &
      nPairOrientation(iTopology,iCanRev), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iCanRev)

   write(u6,'(A,I0,A,I0)') &
      'Reverse x reverse:      pairs = ', &
      nPairOrientation(iTopology,iRevRev), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iRevRev)

end do

write(u6,*)
if (SeedReported) then
   write(u6,'(A,I0,A)') 'Random sampling used, seed = ',iRandSeed, &
                        '  (reproduce with >export MOLCAS_RANDOM_SEED=<seed>)'
   write(u6,*)
end if
write(u6,*) '========================================='
write(u6,*) ' END OF TWO-OPERATOR TEST'
write(u6,*) '========================================='
write(u6,*)

! <I|H|J> from MatEl_Epq and the MatEl_EpqErs_* routines, against sg_h_psi.
if (sum(nDiscrepancies) == 0) call Check_MatEl_H()

call mma_deallocate(Occ)

if (sum(nDiscrepancies) /= 0) then
   write(u6,'(A,I0)') 'SG_Test_EpqErs_Ref: total discrepancies = ',sum(nDiscrepancies)
   call Abend()
end if

contains















#ifdef _NOT_USED_
subroutine Print_CSF_Decoding(CIS,ISYCI)

use Definitions, only: iwp,u6

implicit none

type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp) :: iCSF
integer(kind=iwp) :: MV
integer(kind=iwp) :: ISYUP
integer(kind=iwp) :: IUP
integer(kind=iwp) :: IDWN

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' CSF DECODING'
write(u6,*) '========================================='
write(u6,*)

write(u6,'(A)') &
'  ICSF    MV  ISYUP   IUP  IDWN'

do iCSF=1,CIS%NCSF(ISYCI)

   call Decode_CSF( &
        CIS,ISYCI, &
        iCSF, &
        MV,ISYUP,IUP,IDWN )

   write(u6,'(5I8)') &
        iCSF,MV,ISYUP,IUP,IDWN

end do

write(u6,*)

end subroutine Print_CSF_Decoding

subroutine Print_CI_Blocks(CIS,ISYCI)

use Definitions, only: iwp,u6
use Symmetry_Info, only: Mul

implicit none

type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp) :: MV
integer(kind=iwp) :: ISYUP
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: NUP
integer(kind=iwp) :: NDWN
integer(kind=iwp) :: NCI
integer(kind=iwp) :: IOFF

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' CI BLOCK STRUCTURE'
write(u6,*) '========================================='
write(u6,*)

write(u6,'(A)') &
' MV ISYUP ISYDWN   NUP  NDWN   NCI   IOCSF'

do MV=1,CIS%nMidV

   do ISYUP=1,size(CIS%NOCSF,1)

      NCI = CIS%NOCSF(ISYUP,MV,ISYCI)

      if (NCI == 0) cycle

      ISYDWN = Mul(ISYUP,ISYCI)

      NUP  = CIS%NOW(1,ISYUP ,MV)
      NDWN = CIS%NOW(2,ISYDWN,MV)

      IOFF = CIS%IOCSF(ISYUP,MV,ISYCI)

      write(u6,'(7I8)') &
           MV,ISYUP,ISYDWN, &
           NUP,NDWN,NCI,IOFF

   end do

end do

write(u6,*)

end subroutine Print_CI_Blocks
#endif

#ifdef _DEBUGPRINT_
subroutine Print_CSF_Walks(SGS,CIS,ISYCI)

use Symmetry_Info, only: Mul
use Definitions, only: iwp,u6

implicit none

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp) :: MV
integer(kind=iwp) :: ISYUP
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: IUP
integer(kind=iwp) :: IDWN

integer(kind=iwp) :: NUP
integer(kind=iwp) :: NDWN

integer(kind=iwp) :: IUW0
integer(kind=iwp) :: IDW0

integer(kind=iwp) :: ICUPOS
integer(kind=iwp) :: ICDPOS

integer(kind=iwp) :: ICUP
integer(kind=iwp) :: ICDWN

integer(kind=iwp) :: LEV
integer(kind=iwp) :: NNN
integer(kind=iwp) :: IC1, I

integer(kind=iwp) :: ICS(SGS%nLev)

character(len=1), parameter :: CODE(0:3)=['0','u','d','2']
character(len=200) :: LINE

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' CSF WALKS'
write(u6,*) '========================================='
write(u6,*)

I=0
do MV=1,CIS%nMidV

   do ISYUP=1,SGS%nSym

      ISYDWN = Mul(ISYUP,ISYCI)

      NUP  = CIS%NOW(1,ISYUP ,MV)
      NDWN = CIS%NOW(2,ISYDWN,MV)

      if (NUP*NDWN == 0) cycle

      write(u6,*)
      write(u6,'(A,4I6)') &
         'MV,ISYUP,ISYDWN,NUP=', &
         MV,ISYUP,ISYDWN,NUP

      write(u6,'(A,I6)') 'NDWN=',NDWN

      IUW0 = 1-CIS%nIpWlk + CIS%IOW(1,ISYUP ,MV)
      IDW0 = 1-CIS%nIpWlk + CIS%IOW(2,ISYDWN,MV)

      do IDWN=1,NDWN

         ICDPOS = IDW0 + IDWN*CIS%nIpWlk
         ICDWN  = CIS%ICase(ICDPOS)

         NNN = 0

         do LEV=1,SGS%MidLev

            NNN = NNN + 1

            if (NNN == nPack+1) then
               NNN = 1
               ICDPOS = ICDPOS + 1
               ICDWN  = CIS%ICase(ICDPOS)
            end if

            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1

         end do

         do IUP=1,NUP

            ICUPOS = IUW0 + IUP*CIS%nIpWlk
            ICUP   = CIS%ICase(ICUPOS)

            NNN = 0

            do LEV=SGS%MidLev+1,SGS%nLev

               NNN = NNN + 1

               if (NNN == nPack+1) then
                  NNN = 1
                  ICUPOS = ICUPOS + 1
                  ICUP   = CIS%ICase(ICUPOS)
               end if

               IC1 = ICUP/4
               ICS(LEV) = ICUP-4*IC1
               ICUP = IC1

            end do
            I = I + 1

            LINE=' '

            write(LINE(1:25),'(A,I3,A,I3,A,I3,A,I3,A,I3,A)') &
              'I=',I,'(',MV,':',ISYUP,':',IUP,'/',IDWN,')'

            do LEV=1,SGS%nLev
               LINE(30+LEV:30+LEV)=CODE(ICS(LEV))
            end do

            write(u6,'(A)') trim(LINE)

         end do

      end do

   end do
end do

write(u6,*)

end subroutine Print_CSF_Walks
#endif

#ifdef _NOT_USED_
subroutine Print_Block_CSF_Map(CIS,ISYCI,MV,ISYUP)

use Definitions, only: iwp,u6
use Symmetry_Info, only: Mul

type(CIStruct), intent(in) :: CIS

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: ISYUP

integer(kind=iwp) :: ISYDWN
integer(kind=iwp) :: NUP
integer(kind=iwp) :: NDWN

integer(kind=iwp) :: IUP
integer(kind=iwp) :: IDWN
integer(kind=iwp) :: iCSF

ISYDWN = Mul(ISYUP,ISYCI)

NUP  = CIS%NOW(1,ISYUP ,MV)
NDWN = CIS%NOW(2,ISYDWN,MV)

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' BLOCK -> CSF MAP'
write(u6,*) '========================================='

do IDWN=1,NDWN
   do IUP=1,NUP

      call Encode_CSF( &
           CIS,ISYCI, &
           MV,ISYUP,IUP,IDWN, &
           iCSF )

      write(u6,'(4I6,A,I6)') &
           MV,ISYUP,IUP,IDWN, &
           '  CSF=',iCSF

   end do
end do

end subroutine Print_Block_CSF_Map
#endif






#ifdef _NOT_USED_
subroutine Print_Epq_Matrix(iState,LSM,p,q)

use Constants, only: Zero
use Definitions, only: iwp, wp, u6

integer(kind=iwp), intent(in) :: iState
integer(kind=iwp), intent(in) :: LSM
integer(kind=iwp), intent(in) :: p,q

integer(kind=iwp) :: I,J
integer(kind=iwp) :: nCSF

real(kind=wp) :: X

nCSF = CIS(iState)%NCSF(LSM)

write(u6,*)
write(u6,*) '========================================='
write(u6,*) ' EPQ MATRIX'
write(u6,*) '========================================='
write(u6,'(A,2I6)') 'p,q=',p,q
write(u6,*)

write(u6,'(A)') &
'I        J                 <I|Epq|J>'

do I=1,nCSF
   do J=1,nCSF

      X = MatEl_Epq_Ref(iState,LSM,I,J,p,q)

      if (abs(X) < 1.0e-12_wp) cycle

      write(u6,'(2I8,F24.15)') I,J,X

   end do
end do

write(u6,*)

end subroutine Print_Epq_Matrix
#endif

function MatEl_Epq_Ref(iState,ISYI,I,J,p,q) result(Value)

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: iwp, wp
use stdalloc, only: mma_allocate, mma_deallocate

integer(kind=iwp), intent(in) :: iState
integer(kind=iwp), intent(in) :: ISYI
integer(kind=iwp), intent(in) :: I,J
integer(kind=iwp), intent(in) :: p,q

real(kind=wp) :: Value

integer(kind=iwp) :: iSYJ, nCSFI, nCSFJ

real(kind=wp), allocatable :: Ket(:)
real(kind=wp), allocatable :: Res(:)
logical(kind=iwp) :: Reuse_Save

Value = Zero

Reuse_Save=EXS(iState)%ReUse_SGTmp
EXS(iState)%ReUse_SGTmp=.False.

ISYJ  = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%iSM(q)),ISYI)

nCSFI = CIS(iState)%NCSF(iSYI)
nCSFJ = CIS(iState)%NCSF(iSYJ)

call mma_allocate(Ket,nCSFI,Label='Ket')
call mma_allocate(Res,nCSFJ,Label='Res')

Ket(:)=Zero
Ket(J)=One    ! C_J=1.00

Res(:)=Zero

call SG_Epq_Psi( &
     SGS(iState), &
     CIS(iState), &
     EXS(iState), &
     p,q,One,ISYI, &
     Ket,Res )

Value = Res(I)

call mma_deallocate(Ket)
call mma_deallocate(Res)

EXS(iState)%ReUse_SGTmp=Reuse_Save

end function MatEl_Epq_Ref








!
! <I|H|J> for I and J of symmetry LSM, in the form used by sg_h_psi:
!
!   H = Sum_{i>=j, G(ij)=1} k_ij Ehat_ij + Sum_{i>=j} Sum_{k>=l, G(kl)=G(ij)} 1/2 (kl|ij) Ehat_kl Ehat_ij
!
! with Ehat_ij = E_ij + E_ji (E_ii for i = j), k_ij = h_ij - 1/2 Sum_k (ik|kj),
! and triangularly packed TU (h) and TUVX ((kl|ij)), as in sg_h_psi.
!
function MatEl_H(I,J,TU,nTU,TUVX,nTUVX) result(H)

  use Index_functions, only: iTri
  use Constants, only: Half

  integer(kind=iwp), intent(in) :: I, J, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: H
  integer(kind=iwp) :: iO, jO, kO, lO, ij, kl, nL
  real(kind=wp) :: OneInt, TwoInt, G

  H = Zero
  nL = SGS(iState)%nLev
  do iO = 1, nL
    do jO = 1, iO
      ij = iTri(iO,jO)
      if (SGS(iState)%ISm(iO) == SGS(iState)%ISm(jO)) then
        TwoInt = Zero
        do kO = 1, nL
          TwoInt = TwoInt+TUVX(iTri(iTri(iO,kO),iTri(kO,jO)))
        end do
        OneInt = TU(ij)-Half*TwoInt
        H = H+OneInt*MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),iO,jO,LSM,I,J,.False.)
        if (iO /= jO) H = H+OneInt*MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),jO,iO,LSM,I,J,.False.)
      end if
      do kO = 1, nL
        do lO = 1, kO
          kl = iTri(kO,lO)
          if (Mul(SGS(iState)%ISm(kO),SGS(iState)%ISm(lO)) /= &
              Mul(SGS(iState)%ISm(iO),SGS(iState)%ISm(jO))) cycle
          G = Half*TUVX(iTri(kl,ij))
          if (G == Zero) cycle
          ! Ehat_kl Ehat_ij
          H = H+G*MatEl_EpqErs(I,J,kO,lO,iO,jO,Occ(:,I),Occ(:,J))
          if (iO /= jO) H = H+G*MatEl_EpqErs(I,J,kO,lO,jO,iO,Occ(:,I),Occ(:,J))
          if (kO /= lO) then
            H = H+G*MatEl_EpqErs(I,J,lO,kO,iO,jO,Occ(:,I),Occ(:,J))
            if (iO /= jO) H = H+G*MatEl_EpqErs(I,J,lO,kO,jO,iO,Occ(:,I),Occ(:,J))
          end if
        end do
      end do
    end do
  end do

end function MatEl_H





!
! Verify MatEl_H against sg_h_psi: with reproducible pseudo-random integrals
! (any values with the permutational symmetry of the triangular packing will
! do), compare <I|H|J> with (H e_J)_I for all I and for all columns J, or a
! random sample of MaxHCol columns. Aborts on a mismatch.
!
subroutine Check_MatEl_H()

  integer(kind=iwp), parameter :: MaxHCol = 20
  real(kind=wp), parameter :: TolH = 1.0e-10_wp
  integer(kind=iwp) :: nL, nTU, nTUVX, nC, nCol, iCol, I, J, k, nTestH, nFailH, nFailD
  integer(kind=int64) :: iInt
  real(kind=wp), allocatable :: TU(:), TUVX(:), Psi(:), Sgm(:)
  real(kind=wp) :: X, Dev, MaxDev, MaxDevD, T0, T1, W, TimH, TimD
  logical(kind=iwp) :: Reuse_Save, AllCols

  ! sg_h_psi uses this module, so it is called through a wrapper with
  ! intrinsic arguments only (an interface block with the derived types of
  ! the host module conflicts with the external procedure, Intel #8000)
  interface
    subroutine sg_h_psi_state(iState,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)
      import :: iwp, wp
      integer(kind=iwp), intent(in) :: iState, nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
      real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
      real(kind=wp), intent(out) :: Sigma(nCSFs)
    end subroutine sg_h_psi_state
  end interface

  nL = SGS(iState)%nLev
  nTU = nL*(nL+1)/2
  nTUVX = nTU*(nTU+1)/2
  nC = CIS(iState)%NCSF(LSM)

  ! Reproducible integrals in [-0.5,0.5), independent of the sampling seed.
  call mma_allocate(TU,nTU,Label='TU_H')
  call mma_allocate(TUVX,nTUVX,Label='TUVX_H')
  iInt = 271828183_int64
  do k = 1, nTU
    iInt = mod(48271_int64*iInt,2147483647_int64)
    TU(k) = real(iInt,kind=wp)/2147483647.0_wp-0.5_wp
  end do
  do k = 1, nTUVX
    iInt = mod(48271_int64*iInt,2147483647_int64)
    TUVX(k) = real(iInt,kind=wp)/2147483647.0_wp-0.5_wp
  end do
  ! new integrals for MatEl_H_Direct
  call Reset_MatEl_H()

  call mma_allocate(Psi,nC,Label='Psi_H')
  call mma_allocate(Sgm,nC,Label='Sgm_H')
  Reuse_Save = EXS(iState)%ReUse_SGTMP
  EXS(iState)%ReUse_SGTMP = .False.

  AllCols = nC <= MaxHCol
  nCol = merge(nC,MaxHCol,AllCols)
  nTestH = 0
  nFailH = 0
  nFailD = 0
  MaxDev = Zero
  MaxDevD = Zero
  TimH = Zero
  TimD = Zero
  do iCol = 1, nCol
    J = merge(iCol,Rand_Index(nC),AllCols)
    Psi(:) = Zero
    Psi(J) = One
    call sg_h_psi_state(iState,Psi,nC,LSM,Sgm,TUVX,nTUVX,TU,nTU)
    do I = 1, nC
      call CWTime(T0,W)
      X = MatEl_H(I,J,TU,nTU,TUVX,nTUVX)
      call CWTime(T1,W)
      TimH = TimH+T1-T0
      Dev = abs(X-Sgm(I))
      MaxDev = max(MaxDev,Dev)
      nTestH = nTestH+1
      if (Dev > TolH*max(One,abs(Sgm(I)))) then
        nFailH = nFailH+1
        if (nFailH <= 10) write(u6,'(A,2I8,3ES22.13)') 'MatEl_H mismatch: I,J,<I|H|J>,sg_h_psi,diff = ', &
                                                     I,J,X,Sgm(I),X-Sgm(I)
      end if
      call CWTime(T0,W)
      X = MatEl_H_Direct(I,J,Occ(:,I),Occ(:,J),TU,nTU,TUVX,nTUVX)
      call CWTime(T1,W)
      TimD = TimD+T1-T0
      Dev = abs(X-Sgm(I))
      MaxDevD = max(MaxDevD,Dev)
      if (Dev > TolH*max(One,abs(Sgm(I)))) then
        nFailD = nFailD+1
        if (nFailD <= 10) write(u6,'(A,2I8,3ES22.13)') 'MatEl_H_Direct mismatch: I,J,<I|H|J>,sg_h_psi,diff = ', &
                                                     I,J,X,Sgm(I),X-Sgm(I)
      end if
    end do
  end do

  EXS(iState)%ReUse_SGTMP = Reuse_Save
  call mma_deallocate(Sgm)
  call mma_deallocate(Psi)
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)

  write(u6,*)
  write(u6,'(A)') ' MatEl_H verification against sg_h_psi (pseudo-random integrals)'
  write(u6,'(A,I0,A,I0,A)') 'Columns tested         : ',nCol,' of ',nC, merge(' (all)   ',' (random)',AllCols)
  write(u6,'(A,I0)')   'Matrix elements tested : ',nTestH
  write(u6,'(A,2I8)')  'Failures (H, Direct)   : ',nFailH,nFailD
  write(u6,'(A,2ES12.4)') 'Maximum abs. deviation : ',MaxDev,MaxDevD
  write(u6,'(A,2F12.3)') 'CPU time (H, Direct)   : ',TimH,TimD
  if ((nFailH /= 0) .or. (nFailD /= 0)) then
    write(u6,'(A)') 'MatEl_H consistency check FAILED.'
    call Abend()
  end if
  write(u6,'(A)') 'MatEl_H consistency check PASSED.'

end subroutine Check_MatEl_H



!> Direct matrix element of a spin-free excitation operator between two CSFs.
!!
!! Computes EIJ = < I | E(IP,IQ) | J >, where I belongs to symmetry
!! Mul(ISm(IP),ISm(IQ),ISYCJ) and J belongs to ISYCJ.
!!
!! This is intentionally a reference implementation:
!!   * no SGTMP cache or saved state;
!!   * no sorting or mutation of ICOUP;
!!   * no packing, tiling, SIMD, GEMM, or scratch lists;
!!   * only the coupling entries needed for the requested CSFs are visited.
!!
!! I and J are 1-based CSF indices within their respective symmetry sectors.
!! CPQ is not included: the returned value is the bare <I|E_pq|J>.


!> Exhaustively compare MatEl_Epq and MatEl_Epq_Ref for one ket symmetry.
!!
!! The routine checks every orbital pair (p,q), every ket CSF J in
!! symmetry ISYCI, and every bra CSF I in the symmetry generated by E_pq.
!! It prints each inconsistency and a final summary.
subroutine Check_MatEl_Epq(istate,ISYCI,Tol,nTest,nFail,MaxDiff)
  use Symmetry_Info, only: Mul
  use Constants,     only: Zero
  use Definitions,   only: wp, iwp
#include "intent.fh"

  integer(kind=iwp), intent(in)  :: ISYCI, iState
  real(kind=wp),     intent(in)  :: Tol
  integer(kind=iwp), intent(out) :: nTest, nFail
  real(kind=wp),     intent(out) :: MaxDiff

  integer(kind=iwp) :: I, J, p, q
  integer(kind=iwp) :: nCSFI, nCSFJ, ISYCJ
  integer(kind=iwp) :: iMax, jMax, pMax, qMax
  real(kind=wp) :: Value, Value_Ref, Diff, Scale
  logical(kind=iwp) :: Reuse_Save, Verbose
  !
  ! Quick zero test: orbital occupations of the CSFs of every symmetry,
  ! OccAll(t,iOff(ISY)+I) = occupation of level t in CSF I of symmetry ISY.
  ! They are unpacked from the step vectors by CSF_Occupation, not taken
  ! from MatEl_Epq, which is what is being tested. (Taking them from the
  ! reference diagonal costs nLev*nCSF**2 and dominated large cases.)
  !
  integer(kind=iwp), allocatable :: OccAll(:,:)
  integer(kind=iwp) :: iOff(8), ISY, nScreen
#ifdef _CHECK_QUICK_ZERO_
  integer(kind=iwp) :: iLev
#endif
  logical(kind=iwp) :: Screened
  ! Sampling of (I,J) as in the two-operator test (parameters of the host).
  integer(kind=iwp) :: nPairsIJ, nDrawMax, iDraw, nTakeZ, nTakeNZ, nSampledPQ
  logical(kind=iwp) :: Sampled

  nTest   = 0_iwp
  nFail   = 0_iwp
  MaxDiff = Zero
  iMax    = 0_iwp
  jMax    = 0_iwp
  pMax    = 0_iwp
  qMax    = 0_iwp
  Verbose = .False.
! Verbose = .True.

  nCSFI = CIS(iState)%NCSF(ISYCI)

  Write (u6,*) '============================'
  Write (u6,*) '   MatEl_Epq Verification   '
  Write (u6,*) '============================'
  Write (u6,*) 'iState,iSYCI: ',iState,ISYCI
  Write (u6,*) 'SGS(iState)%iSM(:): ',SGS(iState)%iSM(:)

  ! Disable reuse of SGTMP
  Reuse_Save = EXS(iState)%ReUse_SGTMP
  EXS(iState)%ReUse_SGTMP = .False.

  iOff(:) = 0
  do ISY = 2, SGS(iState)%nSym
    iOff(ISY) = iOff(ISY-1)+CIS(iState)%NCSF(ISY-1)
  end do
  call mma_allocate(OccAll,SGS(iState)%nLev, &
                    iOff(SGS(iState)%nSym)+CIS(iState)%NCSF(SGS(iState)%nSym), &
                    Label='OccAll')
  do ISY = 1, SGS(iState)%nSym
    do I = 1, CIS(iState)%NCSF(ISY)
      call CSF_Occupation(ISY,I,OccAll(:,iOff(ISY)+I))
#ifdef _CHECK_QUICK_ZERO_
      ! Cross-check the unpacked occupations against the reference diagonal.
      do iLev = 1, SGS(iState)%nLev
        if (OccAll(iLev,iOff(ISY)+I) /= nint(MatEl_Epq_Ref(istate,ISY,I,I,iLev,iLev))) then
          write(u6,'(A,3I7)') 'CSF_Occupation failed: ISY,I,level =',ISY,I,iLev
          call Abend()
        end if
      end do
#endif
    end do
  end do
  nScreen = 0_iwp
  nSampledPQ = 0_iwp

  do q = 1, SGS(iState)%nLev
    do p = 1, SGS(iState)%nLev

      ISYCJ=Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYCI)
      nCSFJ = CIS(iState)%NCSF(ISYCJ)

      ! All (I,J) pairs, or a random sample of them if there are more than
      ! MaxPairs.
      nPairsIJ = nCSFI*nCSFJ
      Sampled  = nPairsIJ > MaxPairs
      ! Cap the draws relative to the block size, so that sampling a block
      ! just above MaxPairs never costs more than two full loops.
      nDrawMax = merge(min(MaxDraw,2*nPairsIJ),nPairsIJ,Sampled)
      nTakeZ   = 0
      nTakeNZ  = 0
      if (Sampled) nSampledPQ = nSampledPQ + 1_iwp
      do iDraw = 1, nDrawMax
          if (Sampled) then
            if (nTakeZ >= QuotaZ .and. nTakeNZ >= QuotaNZ) exit
            I = Rand_Index(nCSFI)
            J = Rand_Index(nCSFJ)
          else
            J = (iDraw-1)/nCSFI+1
            I = mod(iDraw-1,nCSFI)+1
          end if

          ! Skip the reference when <I|Epq|J> is zero by the occupations;
          ! MatEl_Epq must then return exactly zero.
          Screened = Epq_Is_Zero(OccAll(:,iOff(ISYCI)+I),OccAll(:,iOff(ISYCJ)+J),I,J,p,q)
          if (Sampled) then
            if (Screened) then
              if (nTakeZ >= QuotaZ) cycle
              nTakeZ = nTakeZ+1
            else
              if (nTakeNZ >= QuotaNZ) cycle
              nTakeNZ = nTakeNZ+1
            end if
          end if

          ! evaluate <I| Epq |J>

          Value = MatEl_Epq(SGS(iState), &
                            CIS(iState), &
                            EXS(iState), &
                            p,q,ISYCJ,I,J,Verbose)

          if (Screened) then
            Value_Ref = Zero
            nScreen = nScreen + 1_iwp
#ifdef _CHECK_QUICK_ZERO_
            if (abs(MatEl_Epq_Ref(istate,ISYCJ,I,J,p,q)) > Tol) then
              write(u6,'(A,6I7)') 'Epq_Is_Zero failed: ISYCI,I,ISYJ,J,p,q =', &
                                  ISYCI,I,ISYCJ,J,p,q
              call Abend()
            end if
#endif
          else
            Value_Ref = MatEl_Epq_Ref(istate,ISYCJ,I,J,p,q)
          end if

          nTest = nTest + 1_iwp
          Diff  = abs(Value-Value_Ref)
          Scale = max(1.0_wp,abs(Value),abs(Value_Ref))

          if (Diff > MaxDiff) then
            MaxDiff = Diff
            iMax = I
            jMax = J
            pMax = p
            qMax = q
          end if

          if ((Screened .and. Value /= Zero) .or. Diff > Tol*Scale) then
            nFail = nFail + 1_iwp
            write(u6,'(A,6I7,2ES24.14)') &
              'MATEL_EPQ mismatch: ISYCI,I,ISYJ,J,p,q,value,Value_Ref =', &
              ISYCI,I,ISYCJ,J,p,q,Value,Value_Ref
          end if

      end do
    end do
  end do

  EXS(iState)%ReUse_SGTMP = Reuse_Save
  call mma_deallocate(OccAll)

  write(u6,'(/,A,I0)') 'MatEl_Epq tests       : ',nTest
  write(u6,'(A,I0)')   'Screened as zero      : ',nScreen
  write(u6,'(A,I0)')   'Sampled (p,q) blocks  : ',nSampledPQ
  write(u6,'(A,I0)')   'MatEl_Epq failures    : ',nFail
  write(u6,'(A,ES12.4)') 'Relative tolerance    : ',Tol
  if (nFail > 0) then
    write(u6,'(A,ES24.14)') 'Maximum abs. diff.    : ',MaxDiff
    write(u6,'(A,4I8)') 'At I,J,p,q            : ', &
                       iMax,jMax,pMax,qMax
  end if

  if (nFail == 0) then
    write(u6,'(A)') 'MatEl_Epq consistency check PASSED.'
  else
    write(u6,'(A)') 'MatEl_Epq consistency check FAILED.'
  end if
  Write (u6,*) '============================'

end subroutine Check_MatEl_Epq

!
! Random index in 1..n from a Park-Miller generator (x <- 48271 x mod
! (2**31-1)) with the state iRandState of the host. Self-contained, so the
! sequence does not depend on other users of a random number generator.
!
function Rand_Index(n) result(iRand)
  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp) :: iRand
  if (.not. SeedReported) then
    if (SeedFromEnv) then
      write(u6,'(/,A,I0,A)') 'SG_Test_EpqErs_Ref: random sampling of (I,J) pairs, seed = ', &
                             iRandSeed,' (from MOLCAS_RANDOM_SEED)'
    else
      write(u6,'(/,A,I0,A)') 'SG_Test_EpqErs_Ref: random sampling of (I,J) pairs, seed = ', &
                             iRandSeed,' (from the clock)'
    end if
    write(u6,'(A,I0,/)')   'To reproduce this run: >export MOLCAS_RANDOM_SEED=',iRandSeed
    SeedReported = .true.
  end if
  iRandState = mod(48271_int64*iRandState,2147483647_int64)
  iRand = 1+int(mod(iRandState,int(n,kind=int64)),kind=iwp)
end function Rand_Index


!
! Seed of the generator: MOLCAS_RANDOM_SEED if set, otherwise the clock.
! (getSeed is not used, since inside verify it gives a deterministic seed.)
!
subroutine Init_Rand_Seed()
  character(len=72) :: Line
  integer(kind=int64) :: iCount
  call get_environment_variable('MOLCAS_RANDOM_SEED',Line)
  SeedFromEnv = Line /= ' '
  if (SeedFromEnv) then
    read(Line,*) iRandSeed
  else
    call system_clock(iCount)
    iRandSeed = 1_int64+mod(abs(iCount),2147483646_int64)
  end if
  ! The Park-Miller state must lie in 1..2**31-2.
  iRandState = 1_int64+mod(abs(iRandSeed),2147483646_int64)
  SeedReported = .false.
end subroutine Init_Rand_Seed

  !
  ! Quick zero test for <I|E_pq|J> from the occupations OccI, OccJ of I, J.
  ! Returns .true. only if the matrix element is exactly zero:
  !   p == q: E_pp is diagonal in the CSF basis, so zero unless I == J
  !           (the symmetries of I and J are then equal);
  !   p /= q: zero unless n_J(q) >= 1, n_J(p) <= 1 and
  !           n_I(t) = n_J(t) + d(t,p) - d(t,q) for every level t.
  !
  function Epq_Is_Zero(OccI,OccJ,I,J,p,q) result(IsZero)
    integer(kind=iwp), intent(in) :: OccI(:), OccJ(:), I, J, p, q
    logical(kind=iwp) :: IsZero
    integer(kind=iwp) :: t, nI
    IsZero = .true.
    if (p == q) then
      IsZero = I /= J
      return
    end if
    if (OccJ(q) == 0 .or. OccJ(p) == 2) return
    do t = 1, size(OccJ)
      nI = OccJ(t)
      if (t == p) nI = nI+1
      if (t == q) nI = nI-1
      if (nI /= OccI(t)) return
    end do
    IsZero = .false.
  end function Epq_Is_Zero

end subroutine SG_Test_EpqErs_Ref

!=======================================================================
! Matrix elements <I|E_pq|J>, <I|E_pq E_rs|J> and <I|H|J> between single
! CSFs (module level, used by SG_Test_EpqErs_Ref and by SG_MatEl_H). The
! state and the CSF symmetry are ME_iState and ME_LSM.
!=======================================================================

!
! Quick zero test for <I|E_pq E_rs|J>, from the orbital occupations OccI and
! OccJ of I and J (level by level). Returns .true. only if the matrix element is exactly zero:
!
!   1. |K> = E_rs|J> must exist: n_J(s) >= 1 and n_J(r) <= 1 for r /= s,
!      n_J(r) >= 1 for r == s (E_rr|J> = n_r|J>). Likewise E_pq|K> needs
!      n_K(p) >= 1 for p == q (for p /= q this follows from point 2).
!   2. n_I(t) = n_J(t) + d(t,p) - d(t,q) + d(t,r) - d(t,s) for every level t.
!
! A .false. result does not mean the element is nonzero (spin coupling can
! still make it vanish). Spatial symmetry is implied by point 2.
!
function Quick_Zero(OccI,OccJ,p,q,r,s) result(IsZero)

integer(kind=iwp), intent(in) :: OccI(:), OccJ(:), p, q, r, s
logical :: IsZero
integer(kind=iwp) :: t, nK, nI

IsZero = .true.

! Intermediate state K = E_rs|J>.
if (r == s) then
   if (OccJ(r) == 0) return
else
   if (OccJ(s) == 0) return
   if (OccJ(r) == 2) return
end if

! Occupations of I must equal those of E_pq E_rs |J>.
do t=1,SGS(ME_iState)%nLev
   nK = OccJ(t)
   if (t == r) nK = nK+1
   if (t == s) nK = nK-1
   if (p == q .and. t == p .and. nK == 0) return
   nI = nK
   if (t == p) nI = nI+1
   if (t == q) nI = nI-1
   if (nI /= OccI(t)) return
end do

IsZero = .false.

end function Quick_Zero

function MatEl_Upper_Half( &
           ISYUP,MV, &
           IUPI,IUPJ, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IUPI
integer(kind=iwp), intent(in) :: IUPJ
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF
integer(kind=iwp) :: ICP
integer(kind=iwp) :: IFrom
integer(kind=iwp) :: ITo
integer(kind=iwp) :: iValue

logical :: Reverse

Value = Zero

if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return
if (p >  SGS(ME_iState)%nLev) return
if (q >  SGS(ME_iState)%nLev) return

!
! Diagonal upper-half operator.
!
if (p == q) then
   if (IUPI /= IUPJ) return
   Value = Upper_Diagonal_Value( &
              ISYUP,MV,IUPI,p )
   return
end if

Reverse = p < q

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYUP,MV)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYUP,MV)

do ICP=1,NCP
   if (.not. Reverse) then
      IFrom = EXS(ME_iState)%ICoup(1,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(2,IOFF+ICP)
   else
      IFrom = EXS(ME_iState)%ICoup(2,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(1,IOFF+ICP)
   end if

   if (IFrom /= IUPJ) cycle
   if (ITo   /= IUPI) cycle

   iValue = EXS(ME_iState)%ICoup(3,IOFF+ICP)
   Value = Value + EXS(ME_iState)%VTab(iValue)
end do

end function MatEl_Upper_Half

function MatEl_Lower_Half( &
           ISYCI,ISYUP,MV, &
           IDWNI,IDWNJ, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IDWNI
integer(kind=iwp), intent(in) :: IDWNJ
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: ISYDWN
integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF
integer(kind=iwp) :: ICP
integer(kind=iwp) :: IFrom
integer(kind=iwp) :: ITo
integer(kind=iwp) :: iValue

logical :: Reverse


Value = Zero

if (p < 1) return
if (q < 1) return
if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (p == q) then
   if (IDWNI /= IDWNJ) return

   Value = Lower_Diagonal_Value( &
              ISYUP,MV,IDWNI,p )

   return
end if

!
! The lower coupling tables are indexed by the symmetry of
! the lower partial walk.
!
ISYDWN = MUL(ISYCI,ISYUP)

Reverse = p < q
call Decode_Lower_Operator( SGS(ME_iState),p,q,INDEO )

!
! Off-diagonal lower coupling paths are indexed by the
! symmetry carried by the lower partial walk.
!
NCP = EXS(ME_iState)%NOCP(INDEO,ISYDWN,MV)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYDWN,MV)

do ICP=1,NCP


   if (Reverse) then
      IFrom = EXS(ME_iState)%ICoup(2,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(1,IOFF+ICP)
   else
      IFrom = EXS(ME_iState)%ICoup(1,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(2,IOFF+ICP)
   end if

   if (IFrom /= IDWNJ) cycle
   if (ITo   /= IDWNI) cycle

   iValue = EXS(ME_iState)%ICoup(3,IOFF+ICP)
   Value = Value + EXS(ME_iState)%VTab(iValue)

end do

end function MatEl_Lower_Half

function Lower_Diagonal_Value( &
           ISYUP,MV,IDWN,p ) result(Value)

use Constants,     only: Zero
use Definitions,   only: wp, iwp
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IDWN
integer(kind=iwp), intent(in) :: p

real(kind=wp) :: Value

integer(kind=iwp) :: ISYDOWN
integer(kind=iwp) :: IOLW
integer(kind=iwp) :: IPSHFT
integer(kind=iwp) :: IPPOW
integer(kind=iwp) :: LLW
integer(kind=iwp) :: JC
integer(kind=iwp) :: ICS

Value = Zero

!
! Evaluate
!
!     <I_L|E_pp|I_L>
!
! for a lower-half walk.
!
if (p < 1) return
if (p > SGS(ME_iState)%MidLev) return

!
! The symmetry of the lower walk follows from the upper-walk
! symmetry and the symmetry of the complete CI state.
!
ISYDOWN = Mul(ISYUP,ME_LSM)

!
! Starting position of the lower walks in this split-graph block.
!
IOLW = CIS(ME_iState)%IOW(2,ISYDOWN,MV)

!
! Each step vector occupies two packed bits.
!
IPSHFT = 2*(p-1)

!
! Locate the packed integer containing the step vector at p.
!
LLW = 1 + IOLW - CIS(ME_iState)%nIpWlk + &
      IPSHFT/(2*nPack)

IPSHFT = mod(IPSHFT,2*nPack)
IPPOW  = 2**IPSHFT

!
! Extract the step vector belonging to lower walk IDWN.
!
JC = CIS(ME_iState)%ICase( &
        LLW + IDWN*CIS(ME_iState)%nIpWlk )

ICS = mod(JC/IPPOW,4)

!
! Diagonal segment values:
!
!     step 0 -> 0
!     step 1 -> 1
!     step 2 -> 1
!     step 3 -> 2
!
if (ICS /= 0) then
   Value = real((1+ICS)/2,kind=wp)
end if

end function Lower_Diagonal_Value

function MatEl_Lower_One( &
           ISYCI, &
           I,J, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF


logical :: Reverse

Value = Zero

!
! This helper handles one non-diagonal lower operator.
!
if (p == q) return

if (p < 1) return
if (q < 1) return

if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! A lower operator leaves the upper walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! Lower coupling tables are indexed by the symmetry carried
! by the lower partial walk.
!
ISYDWN = Mul(ISYCI,ISYUPI)

Reverse = p < q

call Decode_Lower_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYDWN,MVI)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYDWN,MVI)

! the transition IDWNJ -> IDWNI (walk index)
Value = Coup_Value(EXS(ME_iState),IOFF,NCP,merge(2,1,Reverse),IDWNJ,IDWNI)

end function MatEl_Lower_One


function Upper_Diagonal_Value( &
           ISYUP,MV,IUP,p ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IUP
integer(kind=iwp), intent(in) :: p

real(kind=wp) :: Value

integer(kind=iwp) :: IOUW
integer(kind=iwp) :: IPSHFT
integer(kind=iwp) :: IPPOW
integer(kind=iwp) :: LUW
integer(kind=iwp) :: IC
integer(kind=iwp) :: ICS

Value = Zero

!
! This routine evaluates the upper-walk diagonal matrix element
!
!     <I_U|E_pp|I_U>.
!
if (p <= SGS(ME_iState)%MidLev) return
if (p >  SGS(ME_iState)%nLev)   return

!
! Starting position of the upper walks belonging to this
! split-graph block.
!
IOUW = CIS(ME_iState)%IOW(1,ISYUP,MV)

!
! Each step vector occupies two packed bits.
!
IPSHFT = 2*(p-1-SGS(ME_iState)%MidLev)

!
! Select the packed integer containing the step vector at level p.
!
LUW = 1 + IOUW - CIS(ME_iState)%nIpWlk + &
      IPSHFT/(2*nPack)

IPSHFT = mod(IPSHFT,2*nPack)
IPPOW  = 2**IPSHFT

!
! Extract the step vector of upper walk IUP at level p.
!
IC = CIS(ME_iState)%ICase( &
        LUW + IUP*CIS(ME_iState)%nIpWlk )

ICS = mod(IC/IPPOW,4)

!
! The diagonal segment values are:
!
!     step 0 -> 0
!     step 1 -> 1
!     step 2 -> 1
!     step 3 -> 2
!
if (ICS /= 0) then
   Value = real((1+ICS)/2,kind=wp)
end if

end function Upper_Diagonal_Value

function MatEl_Upper_One( &
           ISYCI, &
           I,J, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF


logical :: Reverse

Value = Zero

!
! This helper handles only a non-diagonal upper operator.
!
if (p == q) return

if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! An upper operator leaves the lower walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

Reverse = p < q

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYUPI,MVI)

if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYUPI,MVI)

! the transition IUPJ -> IUPI (walk index)
Value = Coup_Value(EXS(ME_iState),IOFF,NCP,merge(2,1,Reverse),IUPJ,IUPI)

end function MatEl_Upper_One

function MatEl_EpqErs_Ref(ISYCI,I,J,p,q,r,s) result(Value)

   integer(kind=iwp), intent(in) :: ISYCI,I,J
   integer(kind=iwp), intent(in) :: p,q,r,s

   real(kind=wp) :: Value

   integer(kind=iwp) :: nCSFI,nCSFK
#ifdef _DEBUGPRINT_
   integer(kind=iwp) :: K
   integer(kind=iwp) :: MVK,ISYUPK,IUPK,IDWNK
#endif
   integer(kind=iwp) :: ISYK

   real(kind=wp), allocatable :: Ket(:)
   real(kind=wp), allocatable :: Tmp(:)
   real(kind=wp), allocatable :: Res(:)
   logical(kind=iwp) :: Reuse_Save

   Reuse_Save=EXS(ME_iState)%ReUse_SGTMP
   EXS(ME_iState)%ReUse_SGTMP=.False.
   !----------------------------------------------------------
   ! Symmetry after application of E_rs.
   !
   ! |K> = E_rs |J>
   !
   ! Gamma(K)=Gamma(r)Gamma(s)Gamma(J)
   !----------------------------------------------------------

   ISYK = Mul(SGS(ME_iState)%ISm(r), &
              Mul(SGS(ME_iState)%ISm(s),ISYCI))

   nCSFI = CIS(ME_iState)%NCSF(ISYCI)
   nCSFK = CIS(ME_iState)%NCSF(ISYK)

   call mma_allocate(Ket,nCSFI,Label='Ket')
   call mma_allocate(Tmp,nCSFK,Label='Tmp')
   call mma_allocate(Res,nCSFI,Label='Res')

   Ket(:) = Zero
   Ket(J) = One

   Tmp(:) = Zero

   call SG_Epq_Psi( SGS(ME_iState), &
                    CIS(ME_iState), &
                    EXS(ME_iState), &
                    r,s,One,ISYCI, &
                    Ket,Tmp )

#ifdef _DEBUGPRINT_
   write(6,*)
   write(6,*) 'INTERMEDIATE STATES AFTER E_rs'
   write(6,*) 'I,J,p,q,r,s=',I,J,p,q,r,s
   write(6,*) 'ISYK=',ISYK

   do K=1,nCSFK

      if (abs(Tmp(K)) > 1.0d-12) then

         call Decode_CSF( CIS(ME_iState),ISYK,K, &
                          MVK,ISYUPK,IUPK,IDWNK )

         write(6,'(I6,F16.8,4I6)') &
              K,Tmp(K),MVK,ISYUPK,IUPK,IDWNK

      end if

   end do
#endif

   Res(:) = Zero

   call SG_Epq_Psi( SGS(ME_iState), &
                    CIS(ME_iState), &
                    EXS(ME_iState), &
                    p,q,One,ISYK, &
                    Tmp,Res )

   Value = Res(I)

   call mma_deallocate(Ket)
   call mma_deallocate(Tmp)
   call mma_deallocate(Res)

   EXS(ME_iState)%ReUse_SGTMP=Reuse_Save

end function MatEl_EpqErs_Ref


function MatEl_LL( &
           I,J, &
           INDEO_Left,Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           ISYRS ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: INDEO_Left
integer(kind=iwp), intent(in) :: INDEO_Right

logical, intent(in) :: Reverse_Left
logical, intent(in) :: Reverse_Right

! Symmetry of the right operator, Gamma(r)Gamma(s)
integer(kind=iwp), intent(in) :: ISYRS

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ




Value = Zero

call Decode_ME( &
     ME_LSM, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ME_LSM, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Both operators act entirely on the lower half.
! The upper half-walk is a spectator.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! For
!
!     <I|E_Left E_Right|J>,
!
! the rightmost operator acts first:
!
!     J_L --E_Right--> K_L --E_Left--> I_L.
!
!
! Lower coupling blocks are indexed by the symmetry carried
! by the lower partial walk, not by the upper-walk symmetry.
!
ISYDWN = Mul(ME_LSM,ISYUPI)

Value = Walk_Path_Value( &
           INDEO_Right,Reverse_Right,merge(Mul(ISYRS,ISYDWN),ISYDWN,Reverse_Right), &
           INDEO_Left, Reverse_Left, merge(ISYDWN,Mul(ISYRS,ISYDWN),Reverse_Left), &
           MVI,IDWNJ,IDWNI )

end function MatEl_LL

function MatEl_UU( &
           ISYCI, &
           I,J, &
           INDEO_Left,Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           ISYRS ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: INDEO_Left
integer(kind=iwp), intent(in) :: INDEO_Right

logical, intent(in) :: Reverse_Left
logical, intent(in) :: Reverse_Right

! Symmetry of the right operator, Gamma(r)Gamma(s)
integer(kind=iwp), intent(in) :: ISYRS

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ




Value = Zero

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Both operators act entirely on the upper half.
! The lower half-walk must therefore be unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

!
! For
!
!     <I|E_Left E_Right|J>,
!
! the rightmost operator acts first:
!
!     J_U --E_Right--> K_U --E_Left--> I_U.
!
Value = Walk_Path_Value( &
           INDEO_Right,Reverse_Right,merge(Mul(ISYRS,ISYUPI),ISYUPI,Reverse_Right), &
           INDEO_Left, Reverse_Left, merge(ISYUPI,Mul(ISYRS,ISYUPI),Reverse_Left), &
           MVI,IUPJ,IUPI )

end function MatEl_UU

!
! Sum over the intermediate walk K of the coupling-coefficient products
! along J --E_Right--> K --E_Left--> I, for the walks IFromJ and IToI of a
! half graph (upper or lower) at midvertex MV. The blocks of both operators
! are given by INDEO, orientation (Reverse: transposed stored transitions)
! and the symmetry under which the block is stored (that of its ICoup(1)
! walk). Only the transitions starting at IFromJ are followed, located
! by the walk index (Path_Value).
!
function Walk_Path_Value(INDEO_Right,Reverse_Right,ISYRight,INDEO_Left,Reverse_Left,ISYLeft,MV,IFromJ,IToI) result(Value)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: INDEO_Right, ISYRight, INDEO_Left, ISYLeft, MV, IFromJ, IToI
logical, intent(in) :: Reverse_Right, Reverse_Left
real(kind=wp) :: Value
integer(kind=iwp) :: NRight, NLeft

NRight = EXS(ME_iState)%NOCP(INDEO_Right,ISYRight,MV)
NLeft  = EXS(ME_iState)%NOCP(INDEO_Left, ISYLeft, MV)
Value = Path_Value(EXS(ME_iState), &
                   EXS(ME_iState)%IOCP(INDEO_Right,ISYRight,MV),NRight,merge(2,1,Reverse_Right),IFromJ, &
                   EXS(ME_iState)%IOCP(INDEO_Left, ISYLeft, MV),NLeft, merge(2,1,Reverse_Left), IToI)

end function Walk_Path_Value

subroutine Decode_Upper_Operator(SGS,p,q,INDEO)

type(SGStruct), intent(in) :: SGS
integer(kind=iwp), intent(in)  :: p,q
integer(kind=iwp), intent(out) :: INDEO

logical(kind=iwp) :: Swap

Swap = q > p

INDEO = merge( 2*SGS%nLev + (q*(q-1))/2 + p, &
               2*SGS%nLev + (p*(p-1))/2 + q, &
               Swap )

end subroutine Decode_Upper_Operator

subroutine Decode_Lower_Operator(SGS,p,q,INDEO)

implicit none

type(SGStruct), intent(in) :: SGS

integer(kind=iwp), intent(in)  :: p
integer(kind=iwp), intent(in)  :: q
integer(kind=iwp), intent(out) :: INDEO

integer(kind=iwp) :: iHigh
integer(kind=iwp) :: iLow

!
! The closed-loop table is indexed by the unordered orbital pair.
! Store the larger orbital index in iHigh and the smaller one
! in iLow.
!
iHigh = max(p,q)
iLow  = min(p,q)

INDEO = 2*SGS%nLev + &
        (iHigh*(iHigh-1))/2 + iLow

end subroutine Decode_Lower_Operator

function MatEl_EpqErs_UU( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

integer(kind=iwp) :: INDEO_Left
integer(kind=iwp) :: INDEO_Right

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

logical :: Diagonal_Left
logical :: Diagonal_Right

logical :: Reverse_Left
logical :: Reverse_Right

real(kind=wp) :: XClosed
real(kind=wp) :: WLeft
real(kind=wp) :: WRight

real(kind=wp) :: Value

Value = Zero

!
! Both operators must lie completely in the upper graph.
!
if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return
if (r <= SGS(ME_iState)%MidLev) return
if (s <= SGS(ME_iState)%MidLev) return

Diagonal_Left  = p == q
Diagonal_Right = r == s

!
! Decode the bra and ket CSFs. The decoded upper walks are
! required for diagonal weights.
!
call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Every UU operator product leaves the lower walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

!
! ============================================================
! Case 1: diagonal x diagonal
!
!     <I|E_pp E_rr|J>
!
! Both operators are diagonal, so the complete upper walk must
! also be unchanged.
! ============================================================
!
if (Diagonal_Left .and. Diagonal_Right) then

   if (IUPI /= IUPJ) return

   WLeft = Upper_Diagonal_Value( &
              ISYUPI,MVI,IUPI,p )

   if (WLeft == Zero) return

   WRight = Upper_Diagonal_Value( &
               ISYUPI,MVI,IUPI,r )

   Value = WLeft*WRight

   return

end if

!
! ============================================================
! Case 2: diagonal x closed
!
!     <I|E_pp E_rs|J>
!
! Since E_pp is the left operator, its weight is evaluated on
! the final bra-side upper walk I_U:
!
!     w_p(I_U) <I|E_rs|J>.
! ============================================================
!
if (Diagonal_Left) then

   WLeft = Upper_Diagonal_Value( &
              ISYUPI,MVI,IUPI,p )

   if (WLeft == Zero) return

   XClosed = MatEl_Upper_One( &
                ISYCI, &
                I,J, &
                r,s )

   Value = WLeft*XClosed

   return

end if

!
! ============================================================
! Case 3: closed x diagonal
!
!     <I|E_pq E_rr|J>
!
! Since E_rr is the rightmost operator, its weight is evaluated
! on the initial ket-side upper walk J_U:
!
!     <I|E_pq|J> w_r(J_U).
! ============================================================
!
if (Diagonal_Right) then

   XClosed = MatEl_Upper_One( &
                ISYCI, &
                I,J, &
                p,q )

   if (XClosed == Zero) return

   WRight = Upper_Diagonal_Value( &
               ISYUPJ,MVJ,IUPJ,r )

   Value = XClosed*WRight

   return

end if

!
! ============================================================
! Case 4: closed x closed
!
! Use the already validated general ordered contraction.
! ============================================================
!
Reverse_Left  = p < q
Reverse_Right = r < s

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO_Left )

call Decode_Upper_Operator( &
     SGS(ME_iState),r,s,INDEO_Right )

Value = MatEl_UU( &
           ISYCI, &
           I,J, &
           INDEO_Left, Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)) )

end function MatEl_EpqErs_UU

function MatEl_EpqErs_LL( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

integer(kind=iwp) :: INDEO_Left
integer(kind=iwp) :: INDEO_Right

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

logical :: Diagonal_Left
logical :: Diagonal_Right

logical :: Reverse_Left
logical :: Reverse_Right

real(kind=wp) :: XClosed
real(kind=wp) :: WLeft
real(kind=wp) :: WRight

real(kind=wp) :: Value

Value = Zero

!
! Both operators must lie completely inside the lower graph.
!
if (p < 1 .or. q < 1) return
if (r < 1 .or. s < 1) return

if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (r > SGS(ME_iState)%MidLev) return
if (s > SGS(ME_iState)%MidLev) return

Diagonal_Left  = p == q
Diagonal_Right = r == s

!
! Decode the bra and ket CSFs.
!
call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Every LL product leaves the upper walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! ============================================================
! Case 1: diagonal x diagonal
!
!     <I|E_pp E_rr|J>
!
! Both operators are diagonal. Consequently, the lower walk
! must also remain unchanged.
! ============================================================
!
if (Diagonal_Left .and. Diagonal_Right) then

   if (IDWNI /= IDWNJ) return

   WLeft = Lower_Diagonal_Value( &
              ISYUPI,MVI,IDWNI,p )

   if (WLeft == Zero) return

   WRight = Lower_Diagonal_Value( &
               ISYUPI,MVI,IDWNI,r )

   Value = WLeft*WRight

   return

end if

!
! ============================================================
! Case 2: diagonal x closed
!
!     <I|E_pp E_rs|J>
!
! E_pp is the left operator. Its weight is therefore evaluated
! on the final bra-side lower walk:
!
!     w_p(I_L) <I|E_rs|J>.
! ============================================================
!
if (Diagonal_Left) then

   WLeft = Lower_Diagonal_Value( &
              ISYUPI,MVI,IDWNI,p )

   if (WLeft == Zero) return

   XClosed = MatEl_Lower_One( &
                ISYCI, &
                I,J, &
                r,s )

   Value = WLeft*XClosed

   return

end if

!
! ============================================================
! Case 3: closed x diagonal
!
!     <I|E_pq E_rr|J>
!
! E_rr is the rightmost operator and acts first. Its weight is
! evaluated on the initial ket-side lower walk:
!
!     <I|E_pq|J> w_r(J_L).
! ============================================================
!
if (Diagonal_Right) then

   XClosed = MatEl_Lower_One( &
                ISYCI, &
                I,J, &
                p,q )

   if (XClosed == Zero) return

   WRight = Lower_Diagonal_Value( &
               ISYUPJ,MVJ,IDWNJ,r )

   Value = XClosed*WRight

   return

end if

!
! ============================================================
! Case 4: closed x closed
!
! Use the already validated oriented LL contraction.
! ============================================================
!
Reverse_Left  = p < q
Reverse_Right = r < s

call Decode_Lower_Operator( &
     SGS(ME_iState),p,q,INDEO_Left )

call Decode_Lower_Operator( &
     SGS(ME_iState),r,s,INDEO_Right )

Value = MatEl_LL( &
           I,J, &
           INDEO_Left, Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)) )

end function MatEl_EpqErs_LL

!
! Compute <I|EpqErs|J> as Sum_K <I|Epq|K><K|Ers|J>
!
! Where <I|Epq|K> is an upper closed loop: Upq(Iu,Ku) x Delta(Il,Kl)
! and   <K|Ers|J> is an open loop: A_r(Ku,Ju) x A_s(Kl,Jl)
!
! Which reduce to
!
!  A_s(Il,Jl) x Sum_Ku Upq(Iu,Ku) A_r(Ku,Ju)
!
function MatEl_EpqErs_UO(ISYCJ,I,J,p,q,r_in,s_in) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r_in,s_in

real(kind=wp) :: Value
real(kind=wp) :: XUpper,XLower,XRef

integer(kind=iwp) :: r,s
integer(kind=iwp) :: MVI,MVJ,MVK,MVX,MV
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK
integer(kind=iwp) :: IUI,IUJ
integer(kind=iwp) :: IDI,IDJ
integer(kind=iwp) :: ISYDI,ISYDJ
integer(kind=iwp) :: iBand
integer(kind=iwp) :: INDEO_R,INDEO_S,INDEO_PQ
integer(kind=iwp) :: NOU,NOL,NU
integer(kind=iwp) :: OffOU,OffOL,OffU
integer(kind=iwp) :: ISYCI, ISYP, ISYQ, ISYPQ, ISYR, ISYS, ISYRS, ISYPQRS

logical(kind=iwp) :: Reverse_U
logical(kind=iwp) :: Swap_rs, Verbose

Verbose=.False.
Value=Zero

Swap_rs = s_in > r_in
r = merge(s_in,r_in,Swap_rs)
s = merge(r_in,s_in,Swap_rs)

ISYP =SGS(ME_iState)%ISm(p)
ISYQ =SGS(ME_iState)%ISm(q)
ISYPQ=Mul(ISYP,ISYQ)
ISYR =SGS(ME_iState)%ISm(r)
ISYS =SGS(ME_iState)%ISm(s)
ISYRS=Mul(ISYR,ISYS)
ISYPQRS=Mul(ISYPQ,ISYRS)
ISYCI = Mul(ISYPQRS,ISYCJ)

if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

if (p<=SGS(ME_iState)%MidLev .or. q<=SGS(ME_iState)%MidLev) return
if (r<=SGS(ME_iState)%MidLev .or. s<1 .or. s>SGS(ME_iState)%MidLev) return
if (max(p,q,r)>SGS(ME_iState)%nLev) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYUK=Mul(ISYPQ,ISYUI)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

If (Verbose) Then
   Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
   Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
   Write (u6,*) '        IUP, JUP  =', IUI, IUJ
   Write (u6,*) '        IDN, JDN  =', IDI, IDJ
   Write (u6,*) '        ISYCI, ISYCJ=', ISYCI, ISYCJ
End If


if (p==q) then
   ! <I|EppErs|J> = Sum_K <I|Epp|K><K|Ers|J>
   !
   ! Epp|K> = np|K> gives
   !
   ! <I|EppErs> = <I|Epp|I> <I|Ers|J>
   Value=MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),p,p,ISYCI,I,I,.False.)* &
         MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),r_in,s_in,ISYCJ,I,J,.False.)
   if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
   XRef=MatEl_EpqErs_Ref(ISYCJ,I,J,p,q,r_in,s_in)
   if (abs(Value-XRef) > 1.0e-12_wp) then
      write(u6,*) '<I|EppErs|J> discrepancy:',I,J,p,r_in,s_in,Value,XRef
      call Abend()
   end if
   end if
   return
end if

! For the open-loop, A_r(Ku(m),Ju(m')) and m=MVK find the connecting m'=MVJ
! Since we have delta(Il,Kl) from the closed loop this implies that MVK=MVI.

MVK=MVI
iBand=0
do MV = 1, 2
   MVX=MERGE(EXS(ME_iState)%MVR(MVK,MV),EXS(ME_iState)%MVL(MVK,MV),Swap_rs)
   If (MVX == 0 .or. MVX /= MVJ) cycle
   iBand = MV
end do

If (iBand==0) Return
MV=iBand

INDEO_R=(MV-1)*SGS(ME_iState)%nLev+r  ! Open-Loop Upper

INDEO_S=(MV-1)*SGS(ME_iState)%nLev+s  ! Open-Loop Lower
If (Verbose) Write (u6,*) 'INDEO_R,INDEO_S=',INDEO_R,INDEO_S

! Symmetry relations for the two canonical open-loop factors.

! The open loop E_rs connects J to the intermediate K (not to I), whose
! upper walk has the symmetry ISYUK = Gamma(pq) Gamma(Iu).
if (ISYUJ /= Mul(SGS(ME_iState)%ISm(r),ISYUK)) return

ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

if (ISYDJ /= Mul(SGS(ME_iState)%ISm(s),ISYDI)) return

! ---------------------------------------------------------------
! Upper contraction
!
!   XUpper(Iu(m),Ju(m') =
!
!     sum_Ku U_pq(Iu(m),Ku(m)) x A_r(Ku(m),Ju(m'))
! ---------------------------------------------------------------

! Process: A_r(Ku(m),Ju(m'))
NOU   = MERGE(EXS(ME_iState)%NOCP(INDEO_R,ISYUK,MVK), &
              EXS(ME_iState)%NOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
OffOU = MERGE(EXS(ME_iState)%IOCP(INDEO_R,ISYUK,MVK), &
              EXS(ME_iState)%IOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOU=',NOU
if (NOU==0) return

! Retrieve the upper closed-loop operator E_pq: Upq(Iu,Ku)
Reverse_U=p<q

call Decode_Upper_Operator(SGS(ME_iState),p,q,INDEO_PQ)

! A closed-loop block is indexed by the symmetry of its ICoup(1) (left)
! walk: K in canonical orientation, I if reversed (they differ for
! Gamma(pq) /= 1).
NU  =EXS(ME_iState)%NOCP(INDEO_PQ,merge(ISYUI,ISYUK,Reverse_U),MVI)
OffU=EXS(ME_iState)%IOCP(INDEO_PQ,merge(ISYUI,ISYUK,Reverse_U),MVI)
if (NU==0) return
If (Verbose) Write (u6,*) 'NU=',NU

! Sum_Ku U_pq(Iu,Ku) A_r(Ku,Ju): only the transitions of A_r starting at
! Ju and of U_pq starting at Ku (walk index, Path_Value).
XUpper=Path_Value(EXS(ME_iState),OffOU,NOU,merge(2,1,Swap_rs),IUJ,OffU,NU,merge(2,1,Reverse_U),IUI)
If (Verbose) Write (u6,*) 'XUpper=',XUpper
if (XUpper==Zero) return

! ---------------------------------------------------------------
! Lower canonical open-loop factor
!
!   XLower = A_s(Il(m),Jl(m')).
!
! ---------------------------------------------------------------

NOL  = MERGE(EXS(ME_iState)%NOCP(INDEO_S,ISYDI,MVI), &
             EXS(ME_iState)%NOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
OffOL= MERGE(EXS(ME_iState)%IOCP(INDEO_S,ISYDI,MVI), &
             EXS(ME_iState)%IOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
If (NOL==0) return
If (Verbose) Write (u6,*) 'NOL=',NOL

XLower=Coup_Value(EXS(ME_iState),OffOL,NOL,merge(2,1,Swap_rs),IDJ,IDI)
If (Verbose) Write (u6,*) 'XLower=',XLower

Value=XUpper*XLower
If (Verbose) Write (u6,*) 'Value=',Value

end function MatEl_EpqErs_UO

!
! <I|E_pq E_rs|J> with E_pq open (crossing the midlevel) and E_rs in the
! upper half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|(E_pq E_rs)^+|I> = <J|E_sr E_qp|I>,
! and E_sr E_qp is an upper x open product, so this is MatEl_EpqErs_UO with
! bra and ket exchanged and the indices (p,q,r,s) -> (s,r,q,p). (Unlike the
! UL/LU case the two operators may share an upper level, so they do not
! commute in general.)
!
! MatEl_EpqErs_UO takes the symmetry of its ket, here I:
! Gamma(I) = Gamma(p) Gamma(q) Gamma(r) Gamma(s) Gamma(J).
!
function MatEl_EpqErs_OU(ISYCJ,I,J,p,q,r,s) result(Value)

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s
real(kind=wp) :: Value
integer(kind=iwp) :: ISYCI

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
Value = MatEl_EpqErs_UO(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_OU

!
! Compute <I|EpqErs|J> as Sum_K <I|Epq|K><K|Ers|J>
!
! Where <I|Epq|K> is a lower closed loop: Lpq(Il,Kl) x Delta(Iu,Ku)
! and   <K|Ers|J> is an open loop: A_r(Ku,Ju) x A_s(Kl,Jl)
!
! Which reduce to
!
!  A_r(Iu,Ju) x Sum_Kl Lpq(Il,Kl) A_s(Kl,Jl)
!
! This is the mirror image of MatEl_EpqErs_UO with the roles of the upper
! and lower halves exchanged. ISYCJ is the symmetry of the ket J.
!
function MatEl_EpqErs_LO(ISYCJ,I,J,p,q,r_in,s_in) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r_in,s_in

real(kind=wp) :: Value
real(kind=wp) :: XUpper,XLower,XRef

integer(kind=iwp) :: r,s
integer(kind=iwp) :: MVI,MVJ,MVK,MVX,MV
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK
integer(kind=iwp) :: IUI,IUJ
integer(kind=iwp) :: IDI,IDJ
integer(kind=iwp) :: ISYDI,ISYDJ,ISYDK
integer(kind=iwp) :: iBand
integer(kind=iwp) :: INDEO_R,INDEO_S,INDEO_PQ
integer(kind=iwp) :: NOU,NOL,NL
integer(kind=iwp) :: OffOU,OffOL,OffL
integer(kind=iwp) :: ISYCI, ISYP, ISYQ, ISYPQ, ISYR, ISYS, ISYRS, ISYPQRS

logical(kind=iwp) :: Reverse_L
logical(kind=iwp) :: Swap_rs, Verbose

Verbose=.False.
Value=Zero

! Canonical open loop: r in the upper, s in the lower half.
Swap_rs = s_in > r_in
r = merge(s_in,r_in,Swap_rs)
s = merge(r_in,s_in,Swap_rs)

ISYP =SGS(ME_iState)%ISm(p)
ISYQ =SGS(ME_iState)%ISm(q)
ISYPQ=Mul(ISYP,ISYQ)
ISYR =SGS(ME_iState)%ISm(r)
ISYS =SGS(ME_iState)%ISm(s)
ISYRS=Mul(ISYR,ISYS)
ISYPQRS=Mul(ISYPQ,ISYRS)
ISYCI = Mul(ISYPQRS,ISYCJ)

if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

if (p<1 .or. q<1 .or. p>SGS(ME_iState)%MidLev .or. q>SGS(ME_iState)%MidLev) return
if (r<=SGS(ME_iState)%MidLev .or. r>SGS(ME_iState)%nLev .or. s<1 .or. s>SGS(ME_iState)%MidLev) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

! The intermediate K has the upper walk and the midvertex of I, and the
! lower walk symmetry Gamma(Kl) = Gamma(pq) Gamma(Il).
MVK  =MVI
ISYUK=ISYUI
ISYDK=Mul(ISYPQ,ISYDI)

If (Verbose) Then
   Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
   Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
   Write (u6,*) '        IUP, JUP  =', IUI, IUJ
   Write (u6,*) '        IDN, JDN  =', IDI, IDJ
   Write (u6,*) '        ISYCI, ISYCJ=', ISYCI, ISYCJ
End If

if (p==q) then
   ! E_pp is diagonal in the CSF basis: <I|EppErs|J> = <I|Epp|I> <I|Ers|J>
   Value=MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),p,p,ISYCI,I,I,.False.)* &
         MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),r_in,s_in,ISYCJ,I,J,.False.)
   if (CI_Verify('ANY')) then  ! sanity check, only when verifying CI developments
   XRef=MatEl_EpqErs_Ref(ISYCJ,I,J,p,q,r_in,s_in)
   if (abs(Value-XRef) > 1.0e-12_wp) then
      write(u6,*) '<I|EppErs|J> discrepancy (LO):',I,J,p,r_in,s_in,Value,XRef
      call Abend()
   end if
   end if
   return
end if

! For the open loop A_r(Ku(m),Ju(m')) with m=MVK find the connecting m'=MVJ.
iBand=0
do MV = 1, 2
   MVX=MERGE(EXS(ME_iState)%MVR(MVK,MV),EXS(ME_iState)%MVL(MVK,MV),Swap_rs)
   If (MVX == 0 .or. MVX /= MVJ) cycle
   iBand = MV
end do
If (iBand==0) Return
MV=iBand

INDEO_R=(MV-1)*SGS(ME_iState)%nLev+r  ! Open-Loop Upper
INDEO_S=(MV-1)*SGS(ME_iState)%nLev+s  ! Open-Loop Lower
If (Verbose) Write (u6,*) 'INDEO_R,INDEO_S=',INDEO_R,INDEO_S

! Symmetry relations for the two open-loop factors.
if (ISYUJ /= Mul(ISYR,ISYUK)) return
if (ISYDJ /= Mul(ISYS,ISYDK)) return

! ---------------------------------------------------------------
! Upper open-loop factor
!
!   XUpper = A_r(Iu(m),Ju(m'))      (Ku = Iu)
! ---------------------------------------------------------------

NOU  = MERGE(EXS(ME_iState)%NOCP(INDEO_R,ISYUK,MVK), &
             EXS(ME_iState)%NOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
OffOU= MERGE(EXS(ME_iState)%IOCP(INDEO_R,ISYUK,MVK), &
             EXS(ME_iState)%IOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOU=',NOU
If (NOU==0) return

XUpper=Coup_Value(EXS(ME_iState),OffOU,NOU,merge(2,1,Swap_rs),IUJ,IUI)
If (Verbose) Write (u6,*) 'XUpper=',XUpper
if (XUpper==Zero) return

! ---------------------------------------------------------------
! Lower contraction
!
!   XLower = sum_Kl L_pq(Il(m),Kl(m)) x A_s(Kl(m),Jl(m'))
! ---------------------------------------------------------------

! Process: A_s(Kl(m),Jl(m'))
NOL  = MERGE(EXS(ME_iState)%NOCP(INDEO_S,ISYDK,MVK), &
             EXS(ME_iState)%NOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
OffOL= MERGE(EXS(ME_iState)%IOCP(INDEO_S,ISYDK,MVK), &
             EXS(ME_iState)%IOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOL=',NOL
if (NOL==0) return

! The lower closed-loop operator E_pq: Lpq(Il,Kl). As in MatEl_Epq, the
! table is indexed by the bra (Il) symmetry for p < q and by the ket (Kl)
! symmetry otherwise.
Reverse_L=p<q
call Decode_Lower_Operator(SGS(ME_iState),p,q,INDEO_PQ)
NL  =merge(EXS(ME_iState)%NOCP(INDEO_PQ,ISYDI,MVI), &
           EXS(ME_iState)%NOCP(INDEO_PQ,ISYDK,MVI),Reverse_L)
OffL=merge(EXS(ME_iState)%IOCP(INDEO_PQ,ISYDI,MVI), &
           EXS(ME_iState)%IOCP(INDEO_PQ,ISYDK,MVI),Reverse_L)
If (Verbose) Write (u6,*) 'NL=',NL
if (NL==0) return

! Sum_Kl L_pq(Il,Kl) A_s(Kl,Jl): walk index, Path_Value.
XLower=Path_Value(EXS(ME_iState),OffOL,NOL,merge(2,1,Swap_rs),IDJ,OffL,NL,merge(2,1,Reverse_L),IDI)
If (Verbose) Write (u6,*) 'XLower=',XLower

Value=XUpper*XLower
If (Verbose) Write (u6,*) 'Value=',Value

end function MatEl_EpqErs_LO

!
! <I|E_pq E_rs|J> with E_pq open (crossing the midlevel) and E_rs in the
! lower half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|E_sr E_qp|I>, and E_sr E_qp is a
! lower x open product, so this is MatEl_EpqErs_LO with bra and ket
! exchanged and the indices (p,q,r,s) -> (s,r,q,p). The operators may share
! a lower level, so they do not commute in general.
!
! MatEl_EpqErs_LO takes the symmetry of its ket, here I:
! Gamma(I) = Gamma(p) Gamma(q) Gamma(r) Gamma(s) Gamma(J).
!
function MatEl_EpqErs_OL(ISYCJ,I,J,p,q,r,s) result(Value)

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s
real(kind=wp) :: Value
integer(kind=iwp) :: ISYCI

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
Value = MatEl_EpqErs_LO(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_OL

!
! Compute <I|EpqErs|J> with both E_pq and E_rs open (crossing the midlevel)
! as Sum_K <I|Epq|K><K|Ers|J>, where each open loop factorises into an
! upper and a lower coupling coefficient:
!
!   <I|Epq|K> = A(Iu,Ku) x A(Il,Kl),   <K|Ers|J> = A(Ku,Ju) x A(Kl,Jl)
!
! so that the sum over K separates for each intermediate midvertex m_K:
!
!   Sum_mK [ Sum_Ku A(Iu,Ku) A(Ku,Ju) ] x [ Sum_Kl A(Il,Kl) A(Kl,Jl) ]
!
! m_K runs over the (at most two) midvertices reached from m_I by the band
! of E_pq, from which m_J is reached by the band of E_rs. Both operators may
! have either orientation (upper index first or second); the coupling
! tables are accessed as in the open-loop case of MatEl_Epq. The
! symmetries of the bra, the ket and K are kept separate, so this also holds
! for symmetry-changing operators. ISYCJ is the symmetry of the ket J.
!
function MatEl_EpqErs_OO(ISYCJ,I,J,p,q,r,s) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s

real(kind=wp) :: Value
real(kind=wp) :: XU, XL

integer(kind=iwp) :: MidLev, nLev
integer(kind=iwp) :: ISYCI
integer(kind=iwp) :: MVI,MVJ,MVK,MV1,MV2
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK,ISYDI,ISYDJ,ISYDK
integer(kind=iwp) :: IUI,IUJ,IDI,IDJ
integer(kind=iwp) :: U1,L1,U2,L2
integer(kind=iwp) :: INDEO
integer(kind=iwp) :: N1,Off1,N2,Off2
integer(kind=iwp) :: iK
logical(kind=iwp) :: Swap1, Swap2

Value=Zero
MidLev=SGS(ME_iState)%MidLev
nLev  =SGS(ME_iState)%nLev

! Both operators must be open.
if (min(p,q,r,s) < 1 .or. max(p,q,r,s) > nLev) return
if ((p > MidLev) .eqv. (q > MidLev)) return
if ((r > MidLev) .eqv. (s > MidLev)) return

! Upper (U) and lower (L) index of each operator, and its orientation.
Swap1 = q > p
U1 = max(p,q); L1 = min(p,q)
Swap2 = s > r
U2 = max(r,s); L2 = min(r,s)

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

! Symmetries of the halves of K: E_pq takes K to I.
ISYUK=Mul(SGS(ME_iState)%ISm(U1),ISYUI)
ISYDK=Mul(SGS(ME_iState)%ISm(L1),ISYDI)
! and E_rs takes J to K.
if (ISYUJ /= Mul(SGS(ME_iState)%ISm(U2),ISYUK)) return
if (ISYDJ /= Mul(SGS(ME_iState)%ISm(L2),ISYDK)) return

do MV1 = 1, 2
   ! Intermediate midvertex m_K reached from m_I by E_pq.
   MVK = merge(EXS(ME_iState)%MVR(MVI,MV1),EXS(ME_iState)%MVL(MVI,MV1),Swap1)
   if (MVK == 0) cycle
   ! Band of E_rs from m_K that reaches m_J.
   MV2 = 0
   do iK = 1, 2
      if (merge(EXS(ME_iState)%MVR(MVK,iK),EXS(ME_iState)%MVL(MVK,iK),Swap2) == MVJ) MV2 = iK
   end do
   if (MV2 == 0) cycle

   ! ------------------------------------------------------------
   ! Upper half: XU = Sum_Ku A_pq(Iu,Ku) A_rs(Ku,Ju)
   ! ------------------------------------------------------------
   INDEO = (MV2-1)*nLev+U2
   N2   = merge(EXS(ME_iState)%NOCP(INDEO,ISYUK,MVK),EXS(ME_iState)%NOCP(INDEO,ISYUJ,MVJ),Swap2)
   Off2 = merge(EXS(ME_iState)%IOCP(INDEO,ISYUK,MVK),EXS(ME_iState)%IOCP(INDEO,ISYUJ,MVJ),Swap2)
   if (N2 == 0) cycle
   INDEO = (MV1-1)*nLev+U1
   N1   = merge(EXS(ME_iState)%NOCP(INDEO,ISYUI,MVI),EXS(ME_iState)%NOCP(INDEO,ISYUK,MVK),Swap1)
   Off1 = merge(EXS(ME_iState)%IOCP(INDEO,ISYUI,MVI),EXS(ME_iState)%IOCP(INDEO,ISYUK,MVK),Swap1)
   if (N1 == 0) cycle
   ! Walk index (Path_Value): only the transitions of E_rs starting at the
   ! ket walk IUJ and of E_pq starting at K.
   XU = Path_Value(EXS(ME_iState),Off2,N2,merge(2,1,Swap2),IUJ,Off1,N1,merge(2,1,Swap1),IUI)
   if (XU == Zero) cycle

   ! ------------------------------------------------------------
   ! Lower half: XL = Sum_Kl A_pq(Il,Kl) A_rs(Kl,Jl)
   ! ------------------------------------------------------------
   INDEO = (MV2-1)*nLev+L2
   N2   = merge(EXS(ME_iState)%NOCP(INDEO,ISYDK,MVK),EXS(ME_iState)%NOCP(INDEO,ISYDJ,MVJ),Swap2)
   Off2 = merge(EXS(ME_iState)%IOCP(INDEO,ISYDK,MVK),EXS(ME_iState)%IOCP(INDEO,ISYDJ,MVJ),Swap2)
   if (N2 == 0) cycle
   INDEO = (MV1-1)*nLev+L1
   N1   = merge(EXS(ME_iState)%NOCP(INDEO,ISYDI,MVI),EXS(ME_iState)%NOCP(INDEO,ISYDK,MVK),Swap1)
   Off1 = merge(EXS(ME_iState)%IOCP(INDEO,ISYDI,MVI),EXS(ME_iState)%IOCP(INDEO,ISYDK,MVK),Swap1)
   if (N1 == 0) cycle
   ! Walk index (Path_Value): only the transitions of E_rs starting at the
   ! ket walk IDJ and of E_pq starting at K.
   XL = Path_Value(EXS(ME_iState),Off2,N2,merge(2,1,Swap2),IDJ,Off1,N1,merge(2,1,Swap1),IDI)

   Value = Value+XU*XL
end do

end function MatEl_EpqErs_OO

!
! <I|E_pq E_rs|J> for I and J of symmetry LSM, by topology: each operator is
! U (both levels upper), L (both lower) or O (open). Screened with the exact
! quick zero test first, from the occupations OccI and OccJ of I and J.
!
function MatEl_EpqErs(I,J,p,q,r,s,OccI,OccJ) result(X)

  integer(kind=iwp), intent(in) :: I, J, p, q, r, s, OccI(:), OccJ(:)
  real(kind=wp) :: X

  X = Zero
  if (Quick_Zero(OccI,OccJ,p,q,r,s)) return
  X = MatEl_EpqErs_Top(I,J,p,q,r,s)

end function MatEl_EpqErs

!> <I|E_pq E_rs|J> by the topology of the two operators, without the
!> occupation test of MatEl_EpqErs: for callers whose operators are already
!> allowed by the occupations of I and J (MatEl_H_Direct, H_Double).
function MatEl_EpqErs_Top(I,J,p,q,r,s) result(X)

  integer(kind=iwp), intent(in) :: I, J, p, q, r, s
  real(kind=wp) :: X
  character(len=2) :: Top

  X = Zero
  Top = Op_Class(p,q)//Op_Class(r,s)
  select case (Top)
    case ('UU')
      X = MatEl_EpqErs_UU(ME_LSM,I,J,p,q,r,s)
    case ('LL')
      X = MatEl_EpqErs_LL(ME_LSM,I,J,p,q,r,s)
    case ('UL')
      X = MatEl_EpqErs_UL(ME_LSM,I,J,p,q,r,s)
    case ('LU')
      X = MatEl_EpqErs_LU(ME_LSM,I,J,p,q,r,s)
    case ('UO')
      X = MatEl_EpqErs_UO(ME_LSM,I,J,p,q,r,s)
    case ('OU')
      X = MatEl_EpqErs_OU(ME_LSM,I,J,p,q,r,s)
    case ('LO')
      X = MatEl_EpqErs_LO(ME_LSM,I,J,p,q,r,s)
    case ('OL')
      X = MatEl_EpqErs_OL(ME_LSM,I,J,p,q,r,s)
    case ('OO')
      X = MatEl_EpqErs_OO(ME_LSM,I,J,p,q,r,s)
  end select

end function MatEl_EpqErs_Top

!> 'U', 'L' or 'O' for an operator E_pq with both levels in the upper half,
!> both in the lower half, or crossing the midlevel.
function Op_Class(p,q) result(C)
  integer(kind=iwp), intent(in) :: p, q
  character(len=1) :: C
  if ((p > SGS(ME_iState)%MidLev) .and. (q > SGS(ME_iState)%MidLev)) then
    C = 'U'
  else if ((p <= SGS(ME_iState)%MidLev) .and. (q <= SGS(ME_iState)%MidLev)) then
    C = 'L'
  else
    C = 'O'
  end if
end function Op_Class

!
! <I|H|J> as MatEl_H, but only from the operator strings that can connect J
! to I, as given by the occupation difference: with the creations P and the
! annihilations M (multisets, I <- J) the term E_pq E_rs contributes only if
! {p,r} = P + X and {q,s} = M + X, X of size 2-|P|:
!
!   |P| > 2 : zero
!   |P| = 2 : the orderings of P and M, at most 2 loop walks (H_Double)
!   |P| = 1 : a <- b: <I|E_ab|J> times k_ab and the analytic part of the
!             strings with X = {x} (with [E_ax,E_xb] = E_ab):
!               x = a, b : E_ab E_xx, E_xx E_ab = n_x(J), n_x(I) E_ab
!               x other  : (ab|xx) n_x E_ab, and E_ax E_xb + E_xb E_ax =
!                          E_ab (x empty), -E_ab (x doubly occupied),
!                          2 E_ax E_xb - E_ab (x open: one loop walk)
!   |P| = 0 : the diagonal Sum_t k_tt n_t + 1/2 Sum_xy (xx|yy) n_x n_y for
!             I = J, and the exchange terms E_xy E_yx for x /= y: by loop
!             walks only for pairs of open shells, otherwise (I = J only)
!             analytically
!
! with H = Sum_pq k_pq E_pq + 1/2 Sum_pqrs (pq|rs) E_pq E_rs (the form of
! MatEl_H). The cost per element is O(o) MatEl_EpqErs calls for a single
! and O(o**2) (o open shells) for equal occupations, instead of O(n**4).
!
function MatEl_H_Direct(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX) result(H)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: I, J, OccI(:), OccJ(:), nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: H
  integer(kind=iwp) :: nL, t, d, k, nP, nM, x, y, a, b
  real(kind=wp) :: E1, W, G
  integer(kind=iwp) :: Pl(2), Mi(2)
  logical(kind=iwp) :: Hit

  H = Zero
  nL = SGS(ME_iState)%nLev

  ! creations (Pl) and annihilations (Mi) taking J to I
  nP = 0
  nM = 0
  do t = 1, nL
    d = OccI(t)-OccJ(t)
    do k = 1, d
      nP = nP+1
      if (nP > 2) return
      Pl(nP) = t
    end do
    do k = 1, -d
      nM = nM+1
      if (nM > 2) return
      Mi(nM) = t
    end do
  end do

  select case (nP)
    case (2)
      H = H_Double(I,J,Pl(1),Pl(2),Mi(1),Mi(2),TUVX,nTUVX)
    case (1)
      a = Pl(1)
      b = Mi(1)
      E1 = MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),a,b,ME_LSM,I,J,.False.)
      ! coefficient W of <I|E_ab|J>: it depends only on the occupations, so it
      ! is kept for the CSF pairs of the same configuration pair (ME_HW).
      if (ME_HW_Valid .and. (ME_HW_ab(1) == a) .and. (ME_HW_ab(2) == b)) then
        Hit = all(ME_HW_Occ(1:nL) == OccJ(1:nL))
      else
        Hit = .false.
      end if
      if (Hit) then
        W = ME_HW
      else
        W = H_Int1(a,b,TU,nTU,TUVX,nTUVX)+Half*H_Int2(a,b,a,a,TUVX,nTUVX)*(OccI(a)+OccJ(a))+ &
            Half*H_Int2(a,b,b,b,TUVX,nTUVX)*(OccI(b)+OccJ(b))
        do x = 1, nL
          if ((x == a) .or. (x == b)) cycle
          W = W+H_Int2(a,b,x,x,TUVX,nTUVX)*OccJ(x)+merge(Half,-Half,OccJ(x) == 0)*H_Int2(a,x,x,b,TUVX,nTUVX)
        end do
        call Keep_Occ(ME_HW_Occ,OccJ,nL)
        ME_HW_ab(:) = [a,b]
        ME_HW = W
        ME_HW_Valid = .true.
      end if
      ! open shells x: the spin coupling needs the loop walk
      do x = 1, nL
        if ((OccJ(x) /= 1) .or. (x == a) .or. (x == b)) cycle
        G = H_Int2(a,x,x,b,TUVX,nTUVX)
        H = H+G*MatEl_EpqErs_Top(I,J,a,x,x,b)
      end do
      H = H+W*E1
    case (0)
      ! the part of <I|H|I> that depends only on the occupations (ME_HD):
      ! one-electron and Coulomb terms, and the exchange of the pairs that
      ! are not open-open (see below)
      if (I == J) then
        if (ME_HD_Valid) then
          Hit = all(ME_HD_Occ(1:nL) == OccJ(1:nL))
        else
          Hit = .false.
        end if
        if (Hit) then
          H = ME_HD
        else
          do x = 1, nL
            if (OccJ(x) == 0) cycle
            H = H+H_Int1(x,x,TU,nTU,TUVX,nTUVX)*OccJ(x)
            do y = 1, nL
              H = H+Half*H_Int2(x,x,y,y,TUVX,nTUVX)*OccJ(x)*OccJ(y)
            end do
          end do
          do x = 1, nL
            do y = x+1, nL
              if ((OccJ(x) == 1) .and. (OccJ(y) == 1)) cycle
              H = H+Half*H_Int2(x,y,y,x,TUVX,nTUVX)*(OccJ(x)+OccJ(y)-OccJ(x)*OccJ(y))
            end do
          end do
          call Keep_Occ(ME_HD_Occ,OccJ,nL)
          ME_HD = H
          ME_HD_Valid = .true.
        end if
      end if
      ! Exchange terms. With E_xy E_yx = n_x - 1/2 n_x n_y - 2 S_x.S_y (x /= y)
      ! and S_t = 0 for a closed or empty level, a pair that is not open-open
      ! contributes only to I = J, analytically:
      ! E_xy E_yx + E_yx E_xy = n_x + n_y - n_x n_y.
      ! For an open-open pair E_yx E_xy = E_xy E_yx - (E_xx - E_yy) has the
      ! same matrix element (n_x = n_y), so one loop walk suffices.
      do x = 1, nL
        if (OccJ(x) /= 1) cycle
        do y = x+1, nL
          if (OccJ(y) /= 1) cycle
          H = H+H_Int2(x,y,y,x,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,x,y,y,x)
        end do
      end do
  end select

end function MatEl_H_Direct

! Double excitation a1,a2 <- b1,b2 (a1 <= a2, b1 <= b2, no a equal to a b):
! the sum of 1/2 (pq|rs) <I|E_pq E_rs|J> over the distinct orderings {p,r}
! of {a1,a2} and {q,s} of {b1,b2}. The two operators of each string have no
! common index, so they commute ([E_pq,E_rs] = delta_qr E_ps - delta_ps E_rq
! = 0) and the strings with the order reversed give the same value: at most
! two loop walks.
function H_Double(I,J,a1,a2,b1,b2,TUVX,nTUVX) result(V)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: I, J, a1, a2, b1, b2, nTUVX
  real(kind=wp), intent(in) :: TUVX(nTUVX)
  real(kind=wp) :: V

  if ((a1 == a2) .and. (b1 == b2)) then
    ! the single string E_ab E_ab
    V = Half*H_Int2(a1,b1,a1,b1,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,a1,b1,a1,b1)
  else if ((a1 == a2) .or. (b1 == b2)) then
    ! E_a1b1 E_a2b2 and its reverse
    V = H_Int2(a1,b1,a2,b2,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,a1,b1,a2,b2)
  else
    ! E_a1b1 E_a2b2, E_a1b2 E_a2b1 and their reverses
    V = H_Int2(a1,b1,a2,b2,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,a1,b1,a2,b2)+ &
        H_Int2(a1,b2,a2,b1,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,a1,b2,a2,b1)
  end if

end function H_Double

! iTri of Index_functions, here in the module so that the compiler can
! inline it in the integral lookups of MatEl_H_Direct.
pure function Tri_SG(i,j) result(ij)

  integer(kind=iwp), intent(in) :: i, j
  integer(kind=iwp) :: ij

  ij = merge(j*(j-1)/2+i,i*(i-1)/2+j,j > i)

end function Tri_SG

! (pq|rs) from the triangular packing of TUVX.
function H_Int2(p,q,r,s,TUVX,nTUVX) result(V)

  integer(kind=iwp), intent(in) :: p, q, r, s, nTUVX
  real(kind=wp), intent(in) :: TUVX(nTUVX)
  real(kind=wp) :: V

  V = TUVX(Tri_SG(Tri_SG(p,q),Tri_SG(r,s)))

end function H_Int2

! k_pq = h_pq - 1/2 Sum_k (pk|kq).
function H_Int1(p,q,TU,nTU,TUVX,nTUVX) result(V)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: p, q, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: V
  integer(kind=iwp) :: kO

  V = Zero
  do kO = 1, SGS(ME_iState)%nLev
    V = V+H_Int2(p,kO,kO,q,TUVX,nTUVX)
  end do
  V = TU(Tri_SG(p,q))-Half*V

end function H_Int1

function MatEl_EpqErs_LU( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J
integer(kind=iwp), intent(in) :: p,q,r,s

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ
integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ
integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

real(kind=wp) :: XL
real(kind=wp) :: XU
real(kind=wp) :: Value

Value = Zero

!
! Left operator is lower, right operator is upper.
!
if (p < 1 .or. q < 1) return
if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (r <= SGS(ME_iState)%MidLev) return
if (s <= SGS(ME_iState)%MidLev) return
if (r >  SGS(ME_iState)%nLev) return
if (s >  SGS(ME_iState)%nLev) return

call Decode_ME( &
     ISYCI,I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI,J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! The upper operator changes the upper-walk symmetry by Gamma(r)Gamma(s)
! (and the lower one the lower-walk symmetry by Gamma(p)Gamma(q), the
! same for a totally symmetric product). A coupling block is indexed by
! the symmetry of its ICoup(1) walk: the source (J) in canonical
! orientation, the destination (I) if reversed.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYUPJ)) return

XL = MatEl_Lower_Half( &
        ISYCI,merge(ISYUPI,ISYUPJ,p < q),MVI, &
        IDWNI,IDWNJ, &
        p,q )

if (abs(XL) <= 1.0e-15_wp) return

XU = MatEl_Upper_Half( &
        merge(ISYUPI,ISYUPJ,r < s),MVI, &
        IUPI,IUPJ, &
        r,s )

Value = XL*XU

end function MatEl_EpqErs_LU

function MatEl_EpqErs_UL( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

!
! <I|E_pq E_rs|J> with E_pq in the upper and E_rs in the lower half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|(E_pq E_rs)^+|I> = <J|E_sr E_qp|I>,
! and E_sr E_qp is a lower x upper product, so this is MatEl_EpqErs_LU with
! bra and ket exchanged and the indices (p,q,r,s) -> (s,r,q,p). (Since the
! two operators share no level they also commute, E_pq E_rs = E_rs E_pq.)
!
! ISYCI is the symmetry of both bra and ket, as in MatEl_EpqErs_LU: the
! present implementation only treats operators that individually preserve
! spatial symmetry. Once that is generalised, the symmetry passed here must
! be that of I, the ket of the exchanged matrix element.
!
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

real(kind=wp) :: Value

Value = MatEl_EpqErs_LU(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_UL

function MatEl_Epq(SGS,CIS,EXS,IP,IQ,ISYJ,I_CSF,J_CSF,Verbose) result(EIJ)
  use Symmetry_Info, only: Mul
  use Constants,     only: Zero
  use Definitions,   only: wp, iwp
#include "intent.fh"
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(in) :: EXS
  integer(kind=iwp), intent(in) :: IP, IQ, ISYJ, I_CSF, J_CSF
  logical(kind=iwp) :: Verbose
  real(kind=wp) :: EIJ

  integer(kind=iwp) :: ISYP, ISYQ, ISYPQ, ISYI
  integer(kind=iwp) :: MVI, MVJ, ISYUI, ISYUJ
  integer(kind=iwp) :: IUP, IDN, JUP, JDN
  logical(kind=iwp) :: Swap, Lower, Upper
  integer(kind=iwp) :: INDEO, NCP, LICP
  integer(kind=iwp) :: INDEO_P, INDEO_Q
  integer(kind=iwp) :: MV, MVX, NCP1, NCP2, LICP1, LICP2, iBand
  integer(kind=iwp) :: ISYDI, ISYUC, ISYDC, JP, JQ, ISYDJ
  integer(kind=iwp) :: IOUW, IOLW, IPSHFT, IPPOW, IC, ICS
  real(kind=wp) :: Upper_Factor, Lower_Factor

  EIJ  = Zero
  ISYP = SGS%ISm(IP)      ! Gamma(p)
  ISYQ = SGS%ISm(IQ)      ! Gamma(q)
  ISYPQ= Mul(ISYP,ISYQ)  ! Gamma(pq)
  ISYI = Mul(ISYPQ,ISYJ)  ! Gamma(I)=Gamma(pq) x Gamma(J)
  If (Verbose) Then
     Write (u6,*)
     Write (u6,*) 'Enter: MatEl_epq'
     Write (u6,*) 'Testing: I,p,q,J=',i_CSF,ip,iq,j_CSF
     Write (u6,*) '        ISYI, ISYJ  =', ISYI, ISYJ
     Write (u6,*) '        I_CSF,j_CSF=', I_CSF,j_CSF
     Write (u6,*) '        nCSFi,nCSFj=',CIS%NCSF(ISYI),CIS%NCSF(ISYJ)
  End If

  if (I_CSF < 1 .or. I_CSF > CIS%NCSF(ISYI)) return
  if (J_CSF < 1 .or. J_CSF > CIS%NCSF(ISYJ)) return

  call Decode_CSF(CIS,ISYI,I_CSF,MVI,ISYUI,IUP,IDN)
  call Decode_CSF(CIS,ISYJ,J_CSF,MVJ,ISYUJ,JUP,JDN)
  ISYDI = Mul(ISYI,ISYUI)   ! Gamma(Il) = Gamma(I) x Gamma(Iu)
  ISYDJ = Mul(ISYJ,ISYUJ)   ! Gamma(Jl) = Gamma(J) x Gamma(Ju)

  If (Verbose) Then
     Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
     Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
     Write (u6,*) '        IUP, JUP  =', IUP, JUP
     Write (u6,*) '        IDN, JDN  =', IDN, JDN
     Write (u6,*) '        ISYI, ISYJ=', ISYI, ISYJ
  End If

  if (IP == IQ) then
    If (Verbose) Write (u6,*) 'Diagonal case'
    ! E_pp is diagonal in this representation.

    if (ISYI /= ISYJ .or. I_CSF /= J_CSF) return

    if (IP > SGS%MidLev) then
      IOUW = CIS%IOW(1,ISYUI,MVI)
      IPSHFT = 2*(IP-1-SGS%MidLev)
      IC = CIS%ICase(1+IOUW-CIS%nIpWlk+IPSHFT/(2*nPack) &
                     +IUP*CIS%nIpWlk)
    else
      IOLW = CIS%IOW(2,ISYDI,MVI)
      IPSHFT = 2*(IP-1)
      IC = CIS%ICase(1+IOLW-CIS%nIpWlk+IPSHFT/(2*nPack) &
                     +IDN*CIS%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*nPack)
    IPPOW = 2**IPSHFT
    ICS = mod(IC/IPPOW,4)
    if (ICS /= 0) EIJ = real((1+ICS)/2,kind=wp)
    If (Verbose) Then
    Upper_Factor = 1.0_wp
    Lower_Factor = 1.0_wp
    if (IP > SGS%MidLev) then
      Upper_Factor = EIJ
    else
      Lower_Factor = EIJ
    end if
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p=q)'
    write(u6,'(A,5I8)') '  I,J,p,q,ICS       =', &
                        I_CSF,J_CSF,IP,IQ,ICS
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  contribution      =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If
    return
  end if

  Swap = IQ > IP
  If (Verbose) Write (u6,*) 'SWAP=',SWAP
  Lower = merge(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
  Upper = merge(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

  if (Lower) then
    If (Verbose) Write (u6,*) 'Non-diagonal lower closed-loop case'
    ! Both orbital levels are in the lower half.  The upper walk is unchanged.
    if (MVI /= MVJ .or. ISYUI /= ISYUJ .or. IUP /= JUP) return
    INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP, &
                  2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
    If (Verbose) Write (u6,*) ' INDEO_PQ =',INDEO

    ! Gamma(pq) x Gamma(Il)
    ISYDC  = Mul(ISYPQ,ISYDI)
    If (ISYDC /= ISYDJ) Call Abend()

    NCP = merge(EXS%NOCP(INDEO,ISYDI,MVI), &
                EXS%NOCP(INDEO,ISYDJ,MVI),Swap)
    LICP = merge(EXS%IOCP(INDEO,ISYDI,MVI), &
                 EXS%IOCP(INDEO,ISYDJ,MVI),Swap)
    If (Verbose) Write (u6,*) ' NCP, LICP=',NCP,LICP
    Upper_Factor = 1.0_wp
    ! the transition JDN -> IDN (walk index)
    Lower_Factor = Coup_Value(EXS,LICP,NCP,merge(2,1,Swap),JDN,IDN)
    EIJ = Lower_Factor*Upper_Factor
    If (Verbose) Then
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p/=q)'
    write(u6,'(A)') '  case              = lower-half excitation'
    write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
    Write(u6,'(A,L1)') 'Swap=',Swap
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If

  else if (Upper) then
    If (Verbose) Write (u6,*) 'Non-diagonal upper closed-loop case'
    ! Both orbital levels are in the upper half.  The lower walk is unchanged.
    ! Walk indices are local to their symmetry block, so the symmetry of the
    ! lower walk must be compared as well as its index.
    if (MVI /= MVJ .or. ISYDI /= ISYDJ .or. IDN /= JDN) return
    INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP, &
                  2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
    If (Verbose) Write (u6,*) ' INDEO_PQ =',INDEO
    NCP = merge(EXS%NOCP(INDEO,ISYUI,MVI), &
                EXS%NOCP(INDEO,ISYUJ,MVI),Swap)
    LICP = merge(EXS%IOCP(INDEO,ISYUI,MVI), &
                 EXS%IOCP(INDEO,ISYUJ,MVI),Swap)
    If (Verbose) Write (u6,*) ' NCP, LICP=',NCP,LICP
    Lower_Factor = 1.0_wp
    ! the transition JUP -> IUP (walk index)
    Upper_Factor = Coup_Value(EXS,LICP,NCP,merge(2,1,Swap),JUP,IUP)
    EIJ = Lower_Factor*Upper_Factor
    If (Verbose) Then
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p/=q)'
    write(u6,'(A)') '  case              = upper-half excitation'
    write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
    Write(u6,'(A,L1)') 'Swap=',Swap
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If

  else
    If (Verbose) Write (u6,*) 'Non-diagonal open-loop case'
    ! The excitation crosses MIDLEV.  Contract the two coupling tables
    ! directly, retaining only paths from J to I.

    ! Evaluate <I(IUP,IDN,MVI)|Epq|J(JUP,JDN,MVJ)> =
    !          A_p(IUP,MVI,JUP,MVJ) x A_q(IDN,MVI,JDN,MVJ)

    JP=merge(IQ,IP,Swap)
    JQ=merge(IP,IQ,Swap)

    ! Find if the midvertex MVI connects with MVJ for delta b = -1 or +1
    iBand=0
    do MV=1,2
      ! Find the Midvertex to the left of MVI
      MVX=merge(EXS%MVR(MVI,MV),EXS%MVL(MVI,MV),Swap)
      If (Verbose) Then
         Write (u6,*) 'MVX=',MVX
         Write (u6,*) 'MVX/=MVJ',MVX /= MVJ
         Write (u6,*) 'Skip= ',MVX == 0 .or. MVX /= MVJ
      End If
      ! branch out if there are no such vertex or if the index don't match
      if (MVX == 0 .or. MVX /= MVJ) cycle
      iBand=MV
    end do

    If (iBand==0) Return
    MV=iBand
    If (Verbose) Write (u6,*) 'MV=',MV

    !Gamma(Ju) = Gamma(p) x Gamma(Iu)
    ISYUC=merge(Mul(ISYQ,ISYUI),Mul(ISYP,ISYUI),Swap)
    if (ISYUC /= ISYUJ) Return

    ! Process A^p(Iu(MVI),Ju(MVJ))
    ! Accessed as (INDEO_P,ISYUJ,MVJ)
    INDEO_P=(MV-1)*SGS%nLev+JP
    If (Verbose) Write (u6,*) 'INDEO_P=',INDEO_P
    NCP1= merge(EXS%NOCP(INDEO_P,ISYUI,MVI), &
                EXS%NOCP(INDEO_P,ISYUJ,MVJ),Swap)
    LICP1=merge(EXS%IOCP(INDEO_P,ISYUI,MVI), &
                EXS%IOCP(INDEO_P,ISYUJ,MVJ),Swap)
    ! the transition JUP -> IUP (walk index)
    Upper_Factor=Coup_Value(EXS,LICP1,NCP1,merge(2,1,Swap),JUP,IUP)

    If (Upper_Factor==Zero) Return

    !Gamma(Jl) = Gamma(q) x Gamma(Il)
    ISYDC=merge(Mul(ISYP,ISYDI),Mul(ISYQ,ISYDI),Swap)
    if (iSYDC /= ISYDJ) Return

    ! Process A^q(Il(MVI),Jl(MVJ))
    ! Accessed as (INDEO_Q,ISYDJ,MVJ)
    INDEO_Q=(MV-1)*SGS%nLev+JQ
    If (Verbose) Write (u6,*) 'INDEO_Q=',INDEO_Q
    NCP2= merge(EXS%NOCP(INDEO_Q,ISYDI,MVI), &
                EXS%NOCP(INDEO_Q,ISYDJ,MVJ),Swap)
    LICP2=merge(EXS%IOCP(INDEO_Q,ISYDI,MVI), &
                EXS%IOCP(INDEO_Q,ISYDJ,MVJ),Swap)
    ! the transition JDN -> IDN (walk index)
    Lower_Factor=Coup_Value(EXS,LICP2,NCP2,merge(2,1,Swap),JDN,IDN)

    If (Lower_Factor==Zero) Return

    EIJ= Upper_Factor*Lower_Factor

    If (Verbose) Then
       write(u6,'(/,A)') 'MatEl_Epq: OPEN LOOP (p/=q)'
       write(u6,'(A)') '  case              = crossing excitation'
       write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
       Write(u6,'(A,L1)') 'Swap=',Swap
       write(u6,'(A,2I8)') '  MVI, MVX          =',MVI,MVX
       write(u6,'(A,I8)') '  crossing branch   =',MV
       write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
       write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
       write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
       Write (u6,*) 'Exit: MatEl_epq'
       Write (u6,*)
    EndIf

  end if


end function MatEl_Epq

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

!> Set up the single-CSF matrix elements for the CSFs of symmetry LSM of
!> state iState (SG_CSF_Occ, SG_MatEl_H).
subroutine SG_MatEl_Init(iState,LSM)
  integer(kind=iwp), intent(in) :: iState, LSM
  ME_iState = iState
  ME_LSM = LSM
  call Reset_Decode_ME()
  if (.not. allocated(CIS(iState)%ICase)) call MkCOT(SGS(iState),CIS(iState))
  if (.not. allocated(EXS(iState)%USGN)) call MkSgNum(LSM,SGS(iState),CIS(iState),EXS(iState))
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

!> <I|H|J> between the GUGA CSFs I and J with the occupations OccI and OccJ
!> (SG_CSF_Occ), H = Sum_pq h_pq E_pq + 1/2 Sum_pqrs (pq|rs) (E_pq E_rs -
!> delta_qr E_ps) with the triangularly packed TU and TUVX, without core
!> energy (MatEl_H_Direct, after SG_MatEl_Init); I and J of the same
!> symmetry.
function SG_MatEl_H(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX) result(H)
  integer(kind=iwp), intent(in) :: I, J, OccI(:), OccJ(:), nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: H
  H = MatEl_H_Direct(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX)
end function SG_MatEl_H

end module SGUGA
