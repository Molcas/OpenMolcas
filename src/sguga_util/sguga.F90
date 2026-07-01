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

module SGUGA

use Molcas, only: MxLev, MxSym, MxGas
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
private

integer(kind=iwp) :: iq

type SGStruct

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
  integer(kind=iwp) :: NSym = 0, nActEl = 0, IFRAS = 0
  integer(kind=iwp) :: IA0, IB0, IC0, iSpin, nLev, nVert, nVert0, MidLev, MVSta, MVEnd, MXUP, MXDWN
  integer(kind=iwp), allocatable :: ISm(:), DRT(:,:), DRT0(:,:), Down(:,:), Down0(:,:), Up(:,:), Ver(:), MAW(:,:), LTV(:), &
                                    DAW(:,:), RAW(:,:)
  integer(kind=iwp), pointer :: DRTP(:,:), DOWNP(:,:)
  integer(kind=iwp) :: L2ACT(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: LEVEL(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: NRAS(MxSym,MxGAS), NRASEL(MxGAS), nRsPrt=0
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
! MVL(:,:) : Right partner midvertex for delta(b)=-1/+1
! USGN(:,:): map from upper reverese-direct-arc-weight sums to upper walk numbers
! LSGN(:,:): map from lower direct-arc-weight sums to lower walk numbers
! VTAB(:): pool of distinct numerical coupling values.
type EXStruct
  integer(kind=iwp) :: MxEO, nICoup
  integer(kind=iwp), allocatable :: NOCP(:,:,:), IOCP(:,:,:), ICoup(:,:), MVL(:,:), MVR(:,:), USGN(:,:), LSGN(:,:)
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
integer(kind=iwp), parameter :: nSegBase   = 26
integer(kind=iwp), parameter :: nSegWeight = 8
integer(kind=iwp), parameter :: nSegTot    = nSegBase + nSegWeight
integer(kind=iwp), parameter :: nSeg  = nSegBase

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

public :: SGStruct, CIStruct, EXStruct, MkCOT, MkSgNum, SG_Free, SG_Init, SG_Init_Simple, SG_Epq_Psi

! Set nPack to the number of cases (2 bit per case) that can be packed in one integer.
#ifdef SIZE_INITIALIZATION
integer(kind=iwp), parameter:: nPack=Storage_size(1_iwp)/2-1
#elif defined (_I8_)
integer(kind=iwp), parameter:: nPack=32-1
#else
integer(kind=iwp), parameter:: nPack=16-1
#endif

public :: nPack

integer(kind=iwp) :: NCP_Max

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
    call RmVert(SGS)

    ! REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)

    call mkDRT(SGS)

    ! IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED DRT TABLE

  end if

  nullify(SGS%DOWNP,SGS%DRTP)

  ! COMPUTE DOWNCHAIN TABLE AND ARC WEIGHT

  call MKDAW(SGS)

  ! COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS

  call MKRAW(SGS)

  ! COMPUTE LTV TABLES.

  call MKLTV(SGS)

  ! COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.

  call MKMID(SGS)

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

  allocate(vmap(0:maxKey-1))
  vmap(:) = 0

  ! ------------------------------------------
  ! helper: inline key encoding
  ! ------------------------------------------
#define KEY(A,B,C) ((A)*(Bmax+1)*(Cmax+1) + (B)*(Cmax+1) + (C))

  ! TOP vertex
  SGS%DRTP(1,LTAB) = SGS%nLev
  SGS%DRTP(1,NTAB) = SGS%NACTEL
  SGS%DRTP(1,ATAB) = SGS%IA0
  SGS%DRTP(1,BTAB) = SGS%IB0
  SGS%DRTP(1,CTAB) = SGS%IC0

  key = KEY(SGS%IA0,SGS%IB0,SGS%IC0)
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

        key = KEY(ADWN,BDWN,CDWN)

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

  deallocate(vmap)

#undef KEY

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

    SGS%MidLev=Min(SGS%nLev,1)

    MINW = 1000000 ! Initiate to a large number
    do IL=1,SGS%nLev-1 ! Loop from level 1 (one above the root vertex) to the
                       ! level below the head vertex.

      ! SGS%RAW(IV,4) is the number of possible paths from the root vertex up to vertex IV
      ! SGS%DAW(IV,4) is the number of possible paths from the head vertex down to vertex IV
      NW = 0
      do IV=SGS%LTV(IL),SGS%LTV(IL-1)-1  ! Loop over all vertex at this level
        NW = NW+ABS(SGS%RAW(IV,4)-SGS%DAW(IV,4))
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
      SGS%MxUp =Max(SGS%MxUp ,SGS%RAW(MV,4))
      SGS%MxDwn=Max(SGS%MxDwn,SGS%DAW(MV,4))
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
  implicit none
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

        if (ID > 0 .and. SGS%Ver(ID) == 0) then
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
        if (ID > 0 .and. SGS%Ver(ID) == 1) then
          CONN(ID) = 1
        end if
      end do
    end do

    ! remove vertices not connected from above
    do IV = 1, SGS%nVert
      if (SGS%Ver(IV) == 1 .and. CONN(IV) == 0) then
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

subroutine SG_Init(nSym,nActEl,iSpin,SGS,CIS,                      &
                   nRas,nRasEl,nRsPrt,                             &
                   EXS,xLevel,xL2Act,xNLEV,xNSM)

  integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin
  type(SGStruct), intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt),nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xnLev, xNSM(MxLev)
  type(EXStruct), optional, intent(inout) :: EXS

  type(TRStruct) :: TRS

  call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,        &
                      nRas,nRasEl,nRsPrt,               &
                      EXS,xLevel,xL2Act,xnLev,xNSM)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

  call MKMAW(SGS)

  if (present(EXS)) then
!     FORM VARIOUS OFFSET TABLES:

!     CONSTRUCT THE CASE LIST

    call MKCOT(SGS,CIS)

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

    call MKSEG(SGS,CIS,EXS)

    !Create the transition infrastructure.
    Call MkTrans(SGS,CIS,TRS)

    ! Count coupling coefficients in compressed blocks indexed by excitation/operator type, symmetry and midvertex.
    call MkNRCOUP(SGS,CIS,EXS,TRS)

    ! Explicitly generates the compressed coupling tuples '(left walk, right walk, value index)' and compacts repeated
    ! numerical values into 'VTab'.
    call MKCOUP(SGS,CIS,EXS,TRS)

    Call Trans_Free(TRS)
  end if

end subroutine SG_Init

subroutine SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,              &
                          nRas,nRasEl,nRsPrt,                     &
                          EXS,xLevel,xL2Act,xNLEV,xNSM,Do_MkSGUGA)

  integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin
  type(SGStruct), intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  type(EXStruct), optional, intent(inout) :: EXS
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xNLEV, xNSM(MxLev)
  logical(kind=iwp), optional, intent(in) :: Do_MkSGUGA
  integer(kind=iwp) :: iSym

  SGS%IFRAS=0
  If (nRsPrt==3) SGS%IFRAS=1
  Do iSym = 1, nSym
     If (Sum(nRas(iSym,1:nRsPrt))/=0) SGS%IFRAS=SGS%IFRAS+1
  End Do
  SGS%nRasEL(1:nRsPrt)=nRasEl(1:nRsPrt)
  SGS%nRas(1:MxSym,1:nRsPrt)=nRas(1:MxSym,1:nRsPrt)
  SGS%nRsPrt=nRsPrt

  ! Make sure that we start from a clean slate.
  if (present(EXS)) then
   ! Here if the extended parameter list was used.
    call SG_Free(SGS,CIS,EXS)
  else
   ! Here if the terse parameter list was used.
    call SG_Free(SGS,CIS)
  end if

  if (nSym < 1 .or. nSym > 8) then
    write(u6,*) ' SG_Init_Simple: illegal nSym value:',nSym
    call Abend()
  end if
  if (iSpin < 1) then
    write(u6,*) ' SG_Init_Simple: illegal iSpin value:',iSpin
    call Abend()
  end if
  if (nActEl < 0) then
    write(u6,*) ' SG_Init_Simple: illegal nActEl value:',nActEl
    call Abend()
  end if

  SGS%nSym=nSym
  SGS%iSpin=iSpin
  SGS%nActEl=nActEl

  if (present(xLevel)) SGS%Level(:) = xLevel(:)
  if (present(xL2Act)) SGS%L2Act(:) = xL2Act(:)
! Initiate if not already set externally.
  if (SGS%LEVEL(1) == 0) SGS%LEVEL(1:SGS%nLev) = [(iq,iq=1,SGS%nLev)]
  if (SGS%L2Act(1) == 0) SGS%L2Act(1:SGS%nLev) = [(iq,iq=1,SGS%nLev)]

! CREATE THE SYMMETRY INDEX VECTOR

SGS%NLEV = xnLEV
! Allocate Level to Symmetry table ISm:
call mma_allocate(SGS%ISM,SGS%nLev,Label='SGS%ISM')
SGS%ISM(1:SGS%nLev) = xNSM(1:SGS%nLev)

  if (present(Do_MkSGUGA)) then
    if (Do_MkSGUGA) call MkSGUGA(SGS,CIS)
  else
    call MkSGUGA(SGS,CIS)
  end if

end subroutine SG_Init_Simple

subroutine SG_Free(SGS,CIS,EXS)
! PURPOSE: FREE THE SGUGA TABLES

type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), optional, intent(inout) :: EXS

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
nullify(SGS%DRTP,SGS%DOWNP)

call mma_deallocate(CIS%NOW,safe='*')
call mma_deallocate(CIS%IOW,safe='*')
call mma_deallocate(CIS%NCSF,safe='*')
call mma_deallocate(CIS%NOCSF,safe='*')
call mma_deallocate(CIS%IOCSF,safe='*')
call mma_deallocate(CIS%ICase,safe='*')
call mma_deallocate(CIS%VSGM,safe='*')
call mma_deallocate(CIS%IVR,safe='*')
call mma_deallocate(CIS%ISGM,safe='*')

if (present(EXS)) then
 call mma_deallocate(EXS%NOCP,safe='*')
 call mma_deallocate(EXS%IOCP,safe='*')
 call mma_deallocate(EXS%ICoup,safe='*')
 call mma_deallocate(EXS%VTab,safe='*')
 call mma_deallocate(EXS%SGTMP,safe='*')
 call mma_deallocate(EXS%MVL,safe='*')
 call mma_deallocate(EXS%MVR,safe='*')
 call mma_deallocate(EXS%USGN,safe='*')
 call mma_deallocate(EXS%LSGN,safe='*')
 call mma_deallocate(EXS%I1list,safe='*')
 call mma_deallocate(EXS%I2list,safe='*')
 call mma_deallocate(EXS%Xlist,safe='*')
end if

end subroutine SG_Free

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
integer(kind=iwp) :: IHALF, ILND, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, LEV, LEV1, LEV2, MV
integer(kind=iwp) :: INIT, IC, IPOS, L, LL
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3
integer(kind=iwp), allocatable :: SCR(:,:)
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: IS, NUW
# endif

Do INIT=0,1
If (INIT==0) Then
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
Else
!call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr',safe='*')  ! First index referenced by IVERT, SYM, and ISTEP
call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')
! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY
CIS%NOW(:,:,:) = 0
EndIf

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
      IF (ISTP==1 .or. ISTP==2) ISML = SGS%ISm(LEV)

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
      If (INIT==1) Then
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
      End If

      ! BACK UP ONE LEVEL AND EXPLORE NEW WALKS:
      LEV = LEV+1
    end do
  end do
end do

If (INIT==0) THEN
! Generate an off-set array for the CIS%NOW array
Call Mk_IOW(CIS,SGS)
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
End If

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
integer(kind=iwp) :: IA, IAL, IB, IBL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB, IVLT, IVRB, IVRT, LEV, MV, MVLL
integer(kind=iwp) :: INL, IN
real(kind=wp) :: V

call mma_allocate(CIS%IVR,SGS%nVert,2,Label='CIS%IVR')
call mma_allocate(CIS%ISGM,SGS%nVert,nSeg,Label='CIS%ISGM')
call mma_allocate(CIS%VSGM,SGS%nVert,nSeg,Label='CIS%VSGM')
call mma_allocate(EXS%MVL,CIS%nMidV,2,Label='EXS%MVL')
call mma_allocate(EXS%MVR,CIS%nMidV,2,Label='EXS%MVR')

! CONSTRUCT THE IVR TABLE.
! Make list of upper right vertices which connects to a left vertex on the same level
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
      IN = 2*IA+IB
!
!     The number of particles in the node decrease left to right. If the difference between two vertices
!     is higher than 1 then there is no reason to explore further vertices to the right.

      If (INL-IN>1) cycle

!     Test if right vertex is consistent with a tail or an intermediate segment
!     Delta(b) = B(left vertex) - B(right vertex)
!     The intermediate segments are divided into five with delta(b)=-1, and five with delat(b)=+1
!     The tail segments are divided into two with delta(b)=-1, and two with delta(b)=+1
!
!     2a+b=N, where N is the number of electrons in the CSF described by the vertex
!     b=2S, where S is the total electonic spinf of the CSF described by the vertex
!
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
write(u6,1234) (MV,EXS%MVL(MV,1),EXS%MVL(MV,2),MV=1,CIS%nMidV)
1234 format('  MV=',I2,'    UPPER WALKS:',8I6)
write(u6,*) ' MVR TABLE:'
write(u6,1234) (MV,EXS%MVR(MV,1),EXS%MVR(MV,2),MV=1,CIS%nMidV)
write(u6,*)
write(u6,*) ' VERTEX PAIR TABLE IVR IN MKSEG:'
write(u6,1234) (IVL,CIS%IVR(IVL,1),CIS%IVR(IVL,2),IVL=1,SGS%nVert)
#endif

! INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
CIS%ISGM(:,:) = 0
CIS%VSGM(:,:) = Zero

! FOR EACH VERTEX LOOP OVER ALL POSSIBLE SEGMENTS IT CAN BE A PART OF
do IVLT=1,SGS%nVert     ! Upper left vertex
  do ISGT=1,nSeg
    ITT = ITVPT(ISGT)   ! 0-3
    IVRT = IVLT               ! Upper right vertex is the same as the upper left vertex.
    If (ITT==1 .or. ITT==2)IVRT = CIS%IVR(IVLT,ITT)  ! Pick up the associated upper right vertex, depends on delta(b)
    if (IVRT == 0) cycle          ! Branch out if there is no vertex that will contruct the top of the segment.
    ! Pick up the vertex index for the left lower vertex as a function of the case
    IVLB = SGS%Down(IVLT,IC1(ISGT))
    if (IVLB == 0) cycle          ! Branch out if there is no such vertex
    ! Pick up the vertex index for the right lower vertex as a function of the case
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
                     MV, MV1, MV2, MV3, MV4, MV5, MXDWN, MXUP, NDWNS1, NSGMX, NSGTMP, NT1TMP, NT2TMP, NT3TMP, NT4TMP, NT5TMP, &
                     NUPS1, INDEO0
integer(kind=iwp), allocatable :: NRL(:,:,:)
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
Do band=1,nOpenBands
   If (ActiveBand(band)) EXS_OpenBlock=EXS_OpenBlock+SGS%nLev
End Do

EXS%MxEO = NRL_OpenBlock + (SGS%nLev*(SGS%nLev+1))/2
Memory= EXS_OpenBlock + (SGS%nLev*(SGS%nLev+1))/2

call mma_allocate(EXS%NOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
call mma_allocate(EXS%IOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
!call mma_allocate(EXS%NOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
!call mma_allocate(EXS%IOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
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

do IVLT = 1, SGS%MVSta-1    ! This formally denotes the upper left node of an segment
  LEV = SGS%DRT(IVLT,LTAB)

  ! loop over all possible transition classes: upper walk, delta(b)=-1/+1, or lower walk
  ! A segment is classified by its top and bottom transition class
  do ICLASS = 0, TRS%nClass-1

    IT0 = TRS%ITR0(IVLT,ICLASS)   ! Offset to list of valid segments
    NT  = TRS%NTR(IVLT,ICLASS)    ! # of valid segments for this vertex

    if (NT == 0) cycle

    do K = 1, NT                  ! loop over all segment
      ITR  = IT0 + K              ! counter
      ISGT = TRS%ISGT(ITR)        ! segment index
      IVLB = TRS%IVLB(ITR)        ! vertex index of lower left node of the segments
      ITOP = TRS%ITOP(ITR)        ! index of top connection type (delta(b)=+1 or -1), ITOP=ICLASS
      IBOT = TRS%IBOT(ITR)        ! dito lower connection type
      ISYM = TRS%ISYM(ITR)        ! the symmetry index of the segment

      do ITSYM = 1, SGS%nSym      ! loop over the symmetry index of the source vertex
        IBSYM = Mul(ITSYM,ISYM)   ! get the symmetry index of target vertex given the symmetry index of the segment

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
  do LFTSYM=1,SGS%nSym             ! Loop over symmetries
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
        EXS%NOCP(INDEO,LFTSYM,MV) = Max(EXS%NOCP(INDEO,LFTSYM,MV),NRL(LFTSYM,IVLT,INDEO))
      end do

    end do

    ! ---- closed loops ----
    do IP = 1, SGS%nLev
    do IQ = 1, IP
      IPQ = IP*(IP-1)/2 + IQ
      INDEO_EXS = IPQ + EXS_OpenBlock
      INDEO_NRL = IPQ + NRL_OpenBlock
      EXS%NOCP(INDEO_EXS,LFTSYM,MV) = MAX(EXS%NOCP(INDEO_EXS,LFTSYM,MV),NRL(LFTSYM,IVLT,INDEO_NRL))
    end do
    end do

  end do
end do

EXS%nICOup = 0
do INDEO=1,SIZE(EXS%IOCP,1)
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      EXS%IOCP(INDEO,LFTSYM,MV) = EXS%nICOup
      EXS%nICOup = EXS%nICOup +  EXS%NOCP(INDEO,LFTSYM,MV)
    end do
  end do
end do

Call Mk_IOW(CIS,SGS)

call CSFCOUNT(CIS,SGS)

! INSERT FOR USE IN SIGMA ROUTINE

NSGMX  = 1
NDWNS1 = 0  ! Dummy initialization
NSGTMP = MXUP*MXDWN   ! Max size of intermediate sigma block

do MV3=1,CIS%nMidV
  MV1 = EXS%MVL(MV3,2)
  MV2 = EXS%MVL(MV3,1)

  MV4 = EXS%MVR(MV3,1)
  MV5 = EXS%MVR(MV3,2)

  do ISYUS1=1,SGS%nSym
    NUPS1 =    CIS%NOW(1,ISYUS1,MV3)
    do ISYDS1=1,SGS%nSym
      NDWNS1 = CIS%NOW(2,ISYDS1,MV3)
      NSGMX = max(NSGMX,CIS%NOCSF(ISYUS1,MV3,ISYDS1))

      if (MV1 /= 0) then
        NT4TMP = NUPS1*CIS%NOW(2,ISYDS1,MV1)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV1)
        NSGTMP = max(NSGTMP,NT4TMP,NT5TMP)
      end if

      if (MV2 /= 0) then
        NT3TMP = NUPS1*CIS%NOW(2,ISYDS1,MV2)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV2)
        NSGTMP = max(NSGTMP,NT3TMP,NT5TMP)
      end if

      if (MV4 /= 0) then
        NT1TMP = NUPS1*CIS%NOW(2,ISYDS1,MV4)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV4)
        NSGTMP = max(NSGTMP,NT1TMP,NT5TMP)
      end if

      if (MV5 /= 0) then
        NT2TMP = NUPS1*CIS%NOW(2,ISYDS1,MV5)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV5)
        NSGTMP = max(NSGTMP,NT2TMP,NT5TMP)
      end if

    end do
  end do
end do

! Set up the code to wheater or not intermediate sigma vectors are to be reused.
Call mma_maxDBLE(NT1TMP)
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
1234 format('  MV=',I2,'    UPPER WALKS:',8I6)
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

  implicit none
  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(inout) :: CIS
  type(EXStruct), intent(inout) :: EXS
  type(TRStruct), intent(in)    :: TRS

  integer(kind=iwp), allocatable :: VHashKey(:), VHashVal(:)
  integer(kind=iwp) :: hsize
  integer(kind=iwp) :: i, i1, i2
  integer(kind=iwp) :: IAWS, IC, ICL, ICOP, ICR
  integer(kind=iwp) :: IHALF, ILND, INDEO, IP, IPOS, IQ
  integer(kind=iwp) :: IS, ISG, ISGT, ISYM, IT, ITYP, ITYPT
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

  call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
  call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')

  MxEO_Block = nOpenBands * SGS%nLev
  ! Special case
  if (SGS%nLev == 1) then
    NVTAB_FINAL = 0
    call mma_allocate(EXS%VTab,NVTAB_FINAL,Label='EXS%VTab')
    return
  end if

  ! NOW is reused as a counter array and restored later by higher-level logic
  do IHALF = 1, 2
    do MV = 1, CIS%nMidV
      do IS = 1, SGS%nSym
        CIS%NOW(IHALF,IS,MV) = 0
      end do
    end do
  end do

  ! Same idea for NOCP
  do INDEO = 1, SIZE(EXS%NOCP,1)
    do MV = 1, CIS%nMidV
      do IS = 1, SGS%nSym
        EXS%NOCP(INDEO,IS,MV) = 0
      end do
    end do
  end do

  call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
  call mma_allocate(ISGPTH,[1,8],[0,SGS%nLev],Label='ISGPTH')
  call mma_allocate(val,[0,SGS%nLev],Label='val')
  call mma_allocate(VTab,nVTab,Label='VTab')
  hsize = 2*nVTab + 1
  call mma_allocate(VHashKey,[0,hsize-1],Label='VHashKey')
  call mma_allocate(VHashVal,[0,hsize-1],Label='VHashVal')
  VHashKey(:) = -huge(1_iwp)
  VHashVal(:) = 0

  ! Coupling coefficient value table
  NVTAB_FINAL = 2
  VTab(1) =  One
  VTab(2) = -One

  NCHECK = 0

  do IHALF = 1, 2

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

    do IVTOP = IVTSTA, IVTEND

      do ITYP = 0, ITYPMX

        IVRTOP = IVTOP
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
          if (TRS%IPRT(ITR) /= 0) then
            IVRT = CIS%IVR(IVLT,TRS%IPRT(ITR))
          end if

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

            INDEO = SGS%nLev*(IT-1) + IP
            if (IT == 3) INDEO = MxEO_Block + (IP*(IP-1))/2 + IQ

            ICOP = 1 + EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
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

! ============================================================
! M5 BEGIN
! Keep all of the following logic identical to the original
! MKCOUP until the transition-based path traversal is fully
! validated:
!   - VTab lookup / insertion
!   - ICoup write
!   - ILNDW-based Lund renumbering
!   - final EXS%VTab materialization
! ============================================================

! ------------------------------------------------------------
! M5: keep VTab logic identical to the original MKCOUP
! Do not optimize or refactor this yet.
! ------------------------------------------------------------
            C = val(LEV2)

            do i = 1, NVTAB_FINAL
              IVTAB = i
              if (abs(C - VTab(i)) < 1.0e-10_wp) exit
            end do

            if (i > NVTAB_FINAL) then
              NVTAB_FINAL = NVTAB_FINAL + 1
              if (NVTAB_FINAL > nVTab) then
                write(u6,*) 'MKCOUP: NVTAB_FINAL exceeded nVTab'
                call Abend()
              end if
              VTab(NVTAB_FINAL) = C
              IVTAB = NVTAB_FINAL
            end if

! ------------------------------------------------------------
! M5: keep ICoup write logic identical to the original MKCOUP
! ------------------------------------------------------------

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

  NCP_Max=Max(1,MaxVal(EXS%NOCP(:,:,:)))
! Write (u6,*) 'NCP_MAX=',NCP_MAX
! Write (u6,*) 'SIZE(EXS%NOCP)=',SIZE(EXS%NOCP)
! Write (u6,*) 'EXS%NOCP=',EXS%NOCP
  Call mma_allocate(EXS%I1List,NCP_Max,Label='I1List')
  Call mma_allocate(EXS%I2List,NCP_Max,Label='I2List')
  Call mma_allocate(EXS%XList,NCP_Max,Label='XList')
  call mma_deallocate(VHashKey)
  call mma_deallocate(VHashVal)
end subroutine MKCOUP

subroutine MKSGNUM(STSYM,SGS,CIS,EXS)
! PURPOSE: FOR ALL UPPER AND LOWER WALKS
!          COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
!          REVERSE ARC WEIGHT SUM, RESPECTIVELY.
!          STORE THE DATA IN THE TABLES USGN AND LSGN

implicit none
integer(kind=iwp), intent(in) :: STSYM
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IC, ICODE, ICONF, IDAWSUM, ILOFF, ILW, IPOS, IRAWSUM, ISTEP, ISYM, IUOFF, IUW, JPOS, JSYM, LEV, LV, MIDV, &
                     NLW, NUW
integer(kind=iwp), allocatable :: ISTEPVEC(:)

call mma_allocate(EXS%USGN,SGS%MxUp,CIS%nMidV,Label='EXS%USGN')
call mma_allocate(EXS%LSGN,SGS%MxDwn,CIS%nMidV,Label='EXS%LSGN')
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

type(CIStruct), intent(inout) :: CIS
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: ISYDWN, ISYTOT, ISYUP, MV, N



! CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
! SEPARATED BY MIDVERTICES AND SYMMETRY.
! FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.

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
type(CIStruct), intent(inout) :: CIS
type(SGStruct), intent(inout) :: SGS

integer(kind=iwp) :: ISYM, MV

! CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
! SEPARATED FOR EACH MIDVERTEX AND SYMMETRY
!
! CIS%IOW(1/2,ISYM,MV) is an off-set vector associated with CIS%NOW(1/2,ISYM,MV)
! CIS%NOW(1/2,ISYM,MV) is the number of upper/lower walks of symmetry ISYM that ends/starts in mid vertex MV (relative indexation).
! 1/2 indicate if the off-set is for upper or lower walks.
!
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

pure logical(kind=iwp) function NeedsRightPartner(ICLASS)
  integer(kind=iwp), intent(in) :: ICLASS

  NeedsRightPartner = .false.
  If (ICLASS==1 .or. ICLASS==2) NeedsRightPartner = .true.
end function NeedsRightPartner

pure integer(kind=iwp) function PartnerSlot(ICLASS)
  integer(kind=iwp), intent(in) :: ICLASS

  PartnerSlot = 0
  select case (ICLASS)
  case (1)
    PartnerSlot = 1
  case (2)
    PartnerSlot = 2
  end select
end function PartnerSlot

pure integer(kind=iwp) function OpenBand(ICLASS)
  integer(kind=iwp), intent(in) :: ICLASS

  OpenBand = 0
  select case (ICLASS)
  case (1)
    OpenBand = 1
  case (2)
    OpenBand = 2
  case (3)
    OpenBand = 3
  case (4)
    OpenBand = 4
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

  integer(kind=iwp) :: IVLT, ISGT, ICLASS
  integer(kind=iwp) :: LEV, IVLB
  integer(kind=iwp) :: ICL, ICR, ITOP, IBOT
  integer(kind=iwp) :: IPRT, ISYM, IOBAND, IFLAG
  integer(kind=iwp) :: ITR, N
  integer(kind=iwp), allocatable :: IPOS(:,:)
  real(kind=wp) :: VSEG

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

  do IVLT = 1, SGS%nVert
    LEV = SGS%DRT(IVLT,LTAB)

    do ISGT = 1, nSeg

      IVLB = CIS%ISGM(IVLT,ISGT)
      if (IVLB == 0) cycle

      ITOP = ITVPT(ISGT)   ! top connection class
      IBOT = IBVPT(ISGT)   ! bottom connection class

      ICL = IC1(ISGT)      ! left case
      ICR = IC2(ISGT)      ! right cas

      ! symmetry index of the segment
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

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout), target :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
real(kind=wp), intent(in) :: CPQ, CI(*)
real(kind=wp), intent(_OUT_) :: SGM(*)
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX, nCSFs
real(kind=wp) :: X

! declarations to facilitate the reuse of sigma vectors if possible
integer(kind=iwp), save:: i_save_p=0, i_save_q=0
integer(kind=iwp), save:: i_save_p_sym=-1, i_save_q_sym=-1
logical(kind=iwp) :: Reuse_Sigma
real(kind=wp) :: CI_ID=Zero
integer(kind=iwp) ::  iOff, jOff

!***********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!***********************************************************************

nCSFs=CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
ISYP = SGS%ISm(IP)
ISYQ = SGS%ISm(IQ)
ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)

if (IQ < IP) then
  ! THEN THIS IS AN EXCITING OPERATOR.
  if (IP <= SGS%MidLev) then

    ! EXCITING CASE, IQ<IP<=MIDLEV.
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym

        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
        ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)       ; if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
        NCP = EXS%NOCP(INDEO,ISYDC,MVSGM)    ; if (NCP == 0) cycle
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP = EXS%IOCP(INDEO,ISYDC,MVSGM)      ! get the off-set to the block of compressed coupling coefficients.
        ! CASE IS: LOWER HALF, EXCITE:
        call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,.false.)
        call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.false.)

      end do
    end do

  else if (SGS%MidLev < IQ) then

    ! EXCITING CASE, MIDLEV<IQ<IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)  ; if (NS1 == 0) cycle
        ISYUC = Mul(ISYPQ,ISYUSG)
        NUPC = CIS%NOW(1,ISYUC,MVSGM)         ; if (NUPC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYUC,MVSGM)     ; if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUC,MVSGM)
        ! CASE IS: UPPER HALF, EXCITE:
        call apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.false.)
      end do
    end do

  else

    ! EXCITING CASE, IQ<=MIDLEV<IP

    iOff=1
    Reuse_Sigma=.False.
    If (EXS%Reuse_SGTMP .and. i_save_p==IP .and. i_save_q_sym==ISYQ) Then
       Reuse_Sigma = CI_ID==Sum(CI(1:nCSFs))+DBLE(nCSFs)
    End If

    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
        MVX = EXS%MVL(MVSGM,MV) ; if (MVX == 0) cycle
        do ISYUSG=1,SGS%nSym
          NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
          ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
          ISYDSG = Mul(ISYUSG,ISYSGM)
          ISYUC = Mul(ISYP,ISYUSG)
          ISYDC = Mul(ISYQ,ISYDSG)

          NUPC = CIS%NOW(1,ISYUC,MVX)         ; if (NUPC == 0) cycle
          NDWNC = CIS%NOW(2,ISYDC,MVX)        ; if (NDWNC == 0) cycle

          INDEO = merge(IP, SGS%nLev+IP, MV==1)
          NCP1 = EXS%NOCP(INDEO,ISYUC,MVX)  ; if (NCP1 == 0) cycle

          ! CASE IS: UPPER HALF, EXCITE:
          LICP = EXS%IOCP(INDEO,ISYUC,MVX)
          IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

          ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
          If (EXS%Reuse_SGTMP .and. .NOT.Reuse_Sigma) Then
             NTMP = NUPSG*NDWNC
             EXS%SGTMP(iOff:iOff+NTMP-1) = Zero
             call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP+1),swap=.false.)
          End If

          jOff=iOff
          If (EXS%Reuse_SGTMP) iOff = iOff + NUPSG*NDWNC

          INDEO = merge(IQ, SGS%nLev+IQ, MV==1)
          NCP2 = EXS%NOCP(INDEO,ISYDC,MVX) ; if (NCP2 == 0) cycle

          If (.NOT.EXS%Reuse_SGTMP) Then
             NTMP = NUPSG*NDWNC
             EXS%SGTMP(jOff:jOff+NTMP-1) = Zero
             call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP+1),swap=.false.)
          End If

          ! CASE IS: LOWER HALF, EXCITE:
          NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
          LICP = EXS%IOCP(INDEO,ISYDC,MVX)
          call Apply_col(CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1),swap=.false.)

        end do
      end do
    end do
    i_save_p=IP
    i_save_q_sym=ISYQ
    i_save_q=0
    i_save_p_sym=-1
    If (EXS%Reuse_SGTMP .and. .Not.Reuse_Sigma) CI_ID=Sum(CI(1:nCSFs))+DBLE(nCSFs)

  end if
else if (IP < IQ) then
  ! THEN THIS IS A DEEXCITING OPERATOR.
  if (IQ <= SGS%MidLev) then

    ! DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
    INDEO = 2*SGS%nLev+(IQ*(IQ-1))/2+IP
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
        ISYDSG = Mul(ISYUSG,ISYSGM)
        ISYDC = Mul(ISYPQ,ISYDSG)
        NDWNC = CIS%NOW(2,ISYDC,MVSGM)       ; if (NDWNC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYDSG,MVSGM)   ; if (NCP == 0) cycle
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
        ! CASE IS: LOWER HALF, DEEXCITE:
        call sort_icoup_block(EXS%ICOUP(1,LICP+1),NCP,.true.)
        call Apply_col(CPQ,NUPSG,NDWNC,CI(IOC+1),NDWNSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.True.)
      end do
    end do

  else if (SGS%MidLev < IP) then

    ! DEEXCITING OPERATOR, MIDLEV<IP<IQ
    INDEO = 2*SGS%nLev+(IQ*(IQ-1))/2+IP
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)  ; if (NS1 == 0) cycle
        ISYUC = Mul(ISYPQ,ISYUSG)
        NUPC = CIS%NOW(1,ISYUC,MVSGM)         ; if (NUPC == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
        NCP = EXS%NOCP(INDEO,ISYUSG,MVSGM)    ; if (NCP == 0) cycle
        LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
        ! CASE IS: UPPER HALF, DEEXCITE:
        call Apply_row(CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA+1),NCP,EXS%ICOUP(1,LICP+1),swap=.true.)
      end do
    end do

  else

    ! DEEXCITING CASE, IP<=MIDLEV<IQ.
    iOff=1
    Reuse_Sigma=.False.
    If (EXS%Reuse_SGTMP .and. i_save_q==IQ .and. i_save_p_sym==ISYP) Then
       Reuse_Sigma = CI_ID==Sum(CI(1:nCSFs))+DBLE(nCSFs)
    End If

    do MVSGM=1,CIS%nMidV
      do MV = 1, 2
         MVX = EXS%MVR(MVSGM,MV) ; if (MVX == 0) cycle

         do ISYUSG=1,SGS%nSym
           NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
           ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
           NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
           ISYDSG = Mul(ISYUSG,ISYSGM)
           ISYUC = Mul(ISYQ,ISYUSG)
           ISYDC = Mul(ISYP,ISYDSG)

           NUPC = CIS%NOW(1,ISYUC,MVX)          ; if (NUPC == 0) cycle
           NDWNC = CIS%NOW(2,ISYDC,MVX)         ; if (NDWNC == 0) cycle

           INDEO = merge(IQ, SGS%nLev+IQ, MV==1)
           NCP1 = EXS%NOCP(INDEO,ISYUSG,MVSGM)  ; if (NCP1 == 0) cycle

           ! CASE IS: UPPER HALF, DEEXCITE:
           LICP = EXS%IOCP(INDEO,ISYUSG,MVSGM)
           IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

           ! IN CASE OF REUSE COMPUTE THE TEMPORARY SIGMA VECTOR REGARDLESS OF THE NCP2 VALUE.
           If (EXS%Reuse_SGTMP .and. .NOT.Reuse_Sigma) Then
              NTMP = NUPSG*NDWNC
              EXS%SGTMP(iOff:iOff+NTMP-1) = Zero
              call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP+1),swap=.true.)
           End If

           jOff=iOff
           If (EXS%Reuse_SGTMP) iOff = iOff + NUPSG*NDWNC

           INDEO = merge(IP, SGS%nLev+IP, MV==1)
           NCP2 = EXS%NOCP(INDEO,ISYDSG,MVSGM) ; if (NCP2 == 0) cycle

           ! CASE IS: UPPER HALF, DEEXCITE:
           If (.NOT.EXS%Reuse_SGTMP) Then
              NTMP = NUPSG*NDWNC
              EXS%SGTMP(jOff:jOff+NTMP-1) = Zero
              call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP+1),swap=.true.)
           End If

           ! CASE IS: LOWER HALF, DEEXCITE:
           NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
           LICP = EXS%IOCP(INDEO,ISYDSG,MVSGM)
           call Apply_col(CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP+1),swap=.true.)

        end do

      end do
    end do
    i_save_q=IQ
    i_save_p_sym=ISYP
    i_save_p=0
    i_save_q_sym=-1
    If (EXS%Reuse_SGTMP .and. .Not.Reuse_Sigma) CI_ID=Sum(CI(1:nCSFs))+DBLE(nCSFs)

  end if
else
  ! THEN THIS IS A SPECIAL CASE.
  if (IP > SGS%MidLev) then

    ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
    ! IP=IQ>MIDLEV
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
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
          ICS = mod(IC/IPPOW,4) ; if (ICS == 0) cycle
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
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM) ; if (NS1 == 0) cycle
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
          ICS = mod(JC/IPPOW,4) ; if (ICS == 0) cycle
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



subroutine sort_icoup_block(ICOUP, NCP, swap)
  use Definitions, only: iwp
  implicit none
  integer(kind=iwp), intent(in) :: NCP
  integer(kind=iwp), intent(inout) :: ICOUP(3,NCP)
  logical(kind=iwp), intent(in) :: swap
  integer(kind=iwp) :: i, j
  integer(kind=iwp) :: temp(3)

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
     end do
     ICOUP(:,j+1) = temp
  end do
end subroutine sort_icoup_block

subroutine apply_col(CPQ, NUP, NDWNC, CI, NDWNSG, SIGMA, NCP, ICOUP, swap)
  use Definitions, only: wp, iwp
  implicit none
  integer(kind=iwp), intent(in) :: NUP, NDWNC, NDWNSG, NCP
  integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUP,NDWNC)
  real(kind=wp), intent(inout) :: SIGMA(NUP,NDWNSG)
  logical(kind=iwp), intent(in) :: swap

  integer(kind=iwp) :: ICP, start, finish, start2, finish2
  integer(kind=iwp) :: nk, nk2, i, I2, I2b, offset, block

  integer(kind=iwp), parameter :: KBLOCK = 16

  integer(kind=iwp) :: I1a(KBLOCK), I1b(KBLOCK)
  real(kind=wp) :: Xa(KBLOCK), Xb(KBLOCK)
  real(kind=wp) :: ABLOCK(NUP,KBLOCK)
  real(kind=wp) :: W(KBLOCK,2)
  real(kind=wp) :: TEMP(NUP,2)

  real(kind=wp), pointer :: VTAB(:)
  VTAB => EXS%VTab

  ICP = 1
  start2=0
  i2b=0

  if (swap) then

  do while (ICP <= NCP)

    I2 = ICOUP(1,ICP)
    start = ICP
    finish = ICP
    do while (finish <= NCP .and. ICOUP(1,finish)==I2)
      finish = finish + 1
    end do

    nk = finish-start

    if (finish <= NCP) then
      I2b = ICOUP(1,finish)
      start2 = finish
      finish2 = start2
      do while (finish2 <= NCP .and. ICOUP(1,finish2)==I2b)
        finish2 = finish2 + 1
      end do
      nk2 = finish2-start2
    else
      nk2 = -1
    end if

    ! structural GEMM only if perfectly safe AND fits buffers
    if (nk == nk2 .and. nk > 0 .and. nk <= KBLOCK) then
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

        call DGEMM_('N','N', NUP, 2, nk, 1.0_wp, ABLOCK, NUP, W, KBLOCK, 0.0_wp, TEMP, NUP)

        SIGMA(:,I2)  = SIGMA(:,I2)  + TEMP(:,1)
        SIGMA(:,I2b) = SIGMA(:,I2b) + TEMP(:,2)

        ICP = finish2
        cycle
      end if
    end if

    ! fallback blocked GEMV
    do offset = 1, nk, KBLOCK
      block = min(KBLOCK, nk-offset+1)

      do i=1,block
        I1a(i)=ICOUP(2,start+offset+i-2)
        Xa(i)=CPQ*VTAB(ICOUP(3,start+offset+i-2))
        ABLOCK(:,i)=CI(:,I1a(i))
      end do

      call DGEMM_('N','N', NUP, 1, block, 1.0_wp, ABLOCK, NUP, Xa, KBLOCK, 0.0_wp, TEMP, NUP)

      SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
    end do

    ICP = finish

  end do

  else

  do while (ICP <= NCP)

    I2 = ICOUP(2,ICP)
    start = ICP
    finish = ICP
    do while (finish <= NCP .and. ICOUP(2,finish)==I2)
      finish = finish + 1
    end do

    nk = finish-start

    if (finish <= NCP) then
      I2b = ICOUP(2,finish)
      start2 = finish
      finish2 = start2
      do while (finish2 <= NCP .and. ICOUP(2,finish2)==I2b)
        finish2 = finish2 + 1
      end do
      nk2 = finish2-start2
    else
      nk2 = -1
    end if

    if (nk == nk2 .and. nk > 0 .and. nk <= KBLOCK) then
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

        call DGEMM_('N','N', NUP, 2, nk, 1.0_wp, ABLOCK, NUP, W, KBLOCK, 0.0_wp, TEMP, NUP)

        SIGMA(:,I2)  = SIGMA(:,I2)  + TEMP(:,1)
        SIGMA(:,I2b) = SIGMA(:,I2b) + TEMP(:,2)

        ICP = finish2
        cycle
      end if
    end if

    do offset = 1, nk, KBLOCK
      block = min(KBLOCK, nk-offset+1)

      do i=1,block
        I1a(i)=ICOUP(1,start+offset+i-2)
        Xa(i)=CPQ*VTAB(ICOUP(3,start+offset+i-2))
        ABLOCK(:,i)=CI(:,I1a(i))
      end do

      call DGEMM_('N','N', NUP, 1, block, 1.0_wp, ABLOCK, NUP, Xa, KBLOCK, 0.0_wp, TEMP, NUP)

      SIGMA(:,I2) = SIGMA(:,I2) + TEMP(:,1)
    end do

    ICP = finish

  end do

  end if

  VTAB => null()
end subroutine apply_col


subroutine apply_row(CPQ, NDWN, NUPC, CI, NUPSG, SIGMA, NCP, ICOUP, swap)
  integer(kind=iwp), intent(in) :: NDWN, NUPC, NUPSG, NCP
  integer(kind=iwp), intent(in) :: ICOUP(3,NCP)
  real(kind=wp), intent(in) :: CPQ, CI(NUPC,NDWN)
  real(kind=wp), intent(inout) :: SIGMA(NUPSG,NDWN)
  logical(kind=iwp), intent(in) :: swap

  integer(kind=iwp), parameter :: KBLOCK=16
  logical, parameter :: USE_OPT=.true.

  integer(kind=iwp) :: ICP, IDWN, i, j, blk, nblk

  real(kind=wp), pointer :: VTAB(:), XLIST(:)
  integer(kind=iwp), pointer :: I1LIST(:), I2LIST(:)

  real(kind=wp) :: X

  ! small tiles
  real(kind=wp) :: CI_blk(KBLOCK,KBLOCK)

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
!!$OMP SIMD
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
!!$OMP SIMD
          do IDWN=1,nblk
            SIGMA(I2list(i+ICP-1), j+IDWN-1) =               SIGMA(I2list(i+ICP-1), j+IDWN-1) + X*CI_blk(IDWN,ICP)
          end do
        end do

      end do
    end do

  end if

  VTAB => null()

end subroutine apply_row


end subroutine SG_Epq_Psi
end module SGUGA
