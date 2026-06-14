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

module SGUGA

use Molcas, only: MxLev, MxSym, MxGas
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
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

! This lists nSeg different types of segments, i=1,...,nSeg
!  1- 4: segments of the head walk from the loop head to the graph head
!  5- 8: head segments
!  9-13: intermediate segments for the case delta(b)=-1
! 14-18: intermediate segments for the case delta(b)=+1
! 19-22: tails segments
! 23-26: segments of the tail walk from the loop tail to the graph tail

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
integer(kind=iwp), parameter :: nSeg=26
integer(kind=iwp), parameter ::                                                                                                    &
                                ITVPT(nSeg) = [ 0, 0, 0, 0,  0, 0, 0, 0,  1, 1, 1, 1, 1,  2, 2, 2, 2, 2,  1, 1, 2, 2,  3, 3, 3, 3],&
                                IBVPT(nSeg) = [ 0, 0, 0, 0,  1, 1, 2, 2,  1, 1, 2, 1, 1,  2, 2, 1, 2, 2,  3, 3, 3, 3,  3, 3, 3, 3],&
                                IC1(nSeg)   = [ 0, 1, 2, 3,  0, 2, 0, 1,  0, 1, 1, 2, 3,  0, 1, 2, 2, 3,  1, 3, 2, 3,  0, 1, 2, 3],&
                                IC2(nSeg)   = [ 0, 1, 2, 3,  1, 3, 2, 3,  0, 1, 2, 2, 3,  0, 1, 1, 2, 3,  0, 2, 0, 1,  0, 1, 2, 3],&
                                ISVC(nSeg)  = [ 1, 1, 1, 1,  1, 7, 8, 4,  1, 2, 9,10, 2,  1, 2,11,12, 2,  1, 5, 6, 3,  1, 1, 1, 1]


integer(kind=iwp), parameter :: TR_WALK  = 1
integer(kind=iwp), parameter :: TR_OPEN  = 2
integer(kind=iwp), parameter :: TR_MID   = 4
integer(kind=iwp), parameter :: TR_CLOSE = 8
integer(kind=iwp), parameter :: TR_TAIL  = 16

public :: SGStruct, CIStruct, EXStruct, MkCOT, MkSgNum, SG_Free, SG_Init, SG_Init_Simple, TRStruct
!public :: SGStruct, CIStruct, EXStruct, MkCOT, MkSgNum, SG_Free, SG_Init, SG_Init_Simple, TRStruct, MkTrans, Trans_Free

! Set nPack to the number of cases (2 bit per case) that can be packed in one integer.
#ifdef SIZE_INITIALIZATION
integer(kind=iwp), parameter:: nPack=Storage_size(1_iwp)/2-1
#elif defined (_I8_)
integer(kind=iwp), parameter:: nPack=32-1
#else
integer(kind=iwp), parameter:: nPack=16-1
#endif

public :: nPack

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
  ! PURPOSE: CONSTRUCT THE UNRESTRICTED GUGA TABLE


    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: ADWN, AUP, BDWN, BUP, CDWN, CUP, LEV, ISTEP, VDWN, VEND, VSTA, VUP, VNEW
    integer(kind=iwp), parameter:: Steps(4,0:3)=Reshape([0,0,1,0, 0,1,0,1, 1,-1,1,1, 1,0,0,2],[4,4])
#   ifdef _DEBUGPRINT_
    integer(kind=iwp) :: VERT
#   endif

    SGS%DRTP(:,:)=0
    SGS%DownP(:,:) = 0  ! Initiate step as void
    SGS%NACTEL = 2*SGS%IA0+SGS%IB0
    SGS%nLev = SGS%IA0+SGS%IB0+SGS%IC0  ! Number of levels

    ! EACH ROW IN THE DESTINCT ROW TABLE (DRT) CORRESPONDS TO A VERTEX IN THE SHAVITT GRAPH
    ! Each row corresponds to five pretabulated values, which depend of the final spin of
    ! the CSFSsi (S), the total number of electrons(N_n), and the total number of levels -- active orbitals -- (n_n)

    ! a_n + b_n + c_n = n_n
    ! N_n=2a_n + b_n, or a_n = N_n/2 - S
    ! b_n = 2S

    ! SET UP TOP ROW

    SGS%DRTP(1,LTAB) = SGS%nLev   ! Level index, n_n
    SGS%DRTP(1,NTAB) = SGS%NACTEL ! Number of electrons, N_n
    SGS%DRTP(1,ATAB) = SGS%IA0    ! a_n
    SGS%DRTP(1,BTAB) = SGS%IB0    ! b_n
    SGS%DRTP(1,CTAB) = SGS%IC0    ! c_n

    VSTA = 1   ! On the top level we have only one node.
    VEND = 1

#   ifdef _DEBUGPRINT_
    write(u6,*) 'A0,B0,C0,NVERT=',SGS%IA0,SGS%IB0,SGS%IC0,SGS%nVert
#   endif

    ! LOOP OVER ALL LEVELS, TOP-DOWN

    do LEV=SGS%nLev,1,-1               ! Loop over all levels

      ! LOOP OVER VERTICES AT LEVEL LEV

      VNEW=0
      do VUP=VSTA,VEND
        ! Pick up the a, b, and c values of this node
        AUP = SGS%DRTP(VUP,ATAB)
        BUP = SGS%DRTP(VUP,BTAB)
        CUP = SGS%DRTP(VUP,CTAB)

        ! LOOP OVER CASES (STEP VECTORS)
        ! AND STORE ONLY VALID CASE NUMBERS WITH ADDRESSES

        ! 0 <= a_k <= a_n
        ! 0 <= c_k <= c_n
        ! 0 <= b_k <= b_n + Min(a_n-a_k,c_n-c_k)

STEP:   do ISTEP=0,3    ! loop over the step vectors

          ! Check that this node can be reached by this step
          ADWN = AUP-STEPS(1,ISTEP)
          if (ADWN < 0) cycle
          BDWN = BUP-STEPS(2,ISTEP)
          if (BDWN < 0) cycle
          CDWN = CUP-STEPS(3,ISTEP)
          if (CDWN < 0) cycle

          ! Check if this is a new vertex, compare with the new vertices accumulated so far for the level below ...
          DO VDWN = VEND+1, VEND+VNEW
             IF (ADWN==SGS%DRTP(VDWN,ATAB) .AND. BDWN==SGS%DRTP(VDWN,BTAB) .AND. CDWN==SGS%DRTP(VDWN,CTAB)) THEN
                SGS%DownP(VUP,ISTEP)=VDWN
                CYCLE STEP
             END IF
          END DO

          VNEW=VNEW+1
          SGS%DownP(VUP,ISTEP)=VEND+VNEW
          SGS%DRTP(VEND+VNEW,ATAB)=ADWN
          SGS%DRTP(VEND+VNEW,BTAB)=BDWN
          SGS%DRTP(VEND+VNEW,CTAB)=CDWN
        end do STEP

      end do
      VSTA=VEND+1
      VEND=VEND+VNEW

    end do
    ! End of loop over levels.

    ! COMPLETE DRT TABLE BY ADDING NO. OF ORBITALS AND ELECTRONS
    ! INTO THE FIRST AND SECOND COLUMN

    SGS%DRTP(1:VEND,LTAB) = SGS%DRTP(1:VEND,ATAB)+SGS%DRTP(1:VEND,BTAB)+SGS%DRTP(1:VEND,CTAB)
    SGS%DRTP(1:VEND,NTAB) = 2*SGS%DRTP(1:VEND,ATAB)+SGS%DRTP(1:VEND,BTAB)
#   ifdef _DEBUGPRINT_
    do VERT=1,VEND
      write(u6,*) 'DRT0(:,LTAB)=',SGS%DRTP(VERT,LTAB)
      write(u6,*) 'DRT0(:,NTAB)=',SGS%DRTP(VERT,NTAB)
    end do
#   endif

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
! Purpose: Remove vertices from a DRT table.

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, iRO, iSy, IV, L, Lev, N, NCHANGES, NLD, NV
logical(kind=iwp) :: Test
integer(kind=iwp), allocatable :: CONN(:), Lim(:)

! Construct a restricted graph.
call mma_allocate(Lim,SGS%nLev,Label='Lim')
Lim(:) = 0
! Fill in the occupation limit table:
! This tells us what is the minimum of particles that has to be found at the levels when we move from on
! partioning, RASx, to another. Note that we assume that the levels of the RASx spaces are consecutive
! and that these levels are ordered in the sequence RAS1, RAS2, and RAS3 levels.

Lev = 0
do iRO=1,SGS%nRsPrt   ! Loop over the partitionings, RAS1, RAS2, RAS3
  do iSy=1,SGS%nSym   ! Loop over the irreps and accumulate the number of levels associated with the RASx space.
    Lev = Lev+SGS%nRas(iSy,iRO)
  end do
  ! If the number of levels are non-zero set the limit at the level which marks the end of a RASx space to be
  ! equivalent to the number of electrons in that subspace
  if (Lev > 0) Lim(Lev) = SGS%nRasEl(iRO)
end do

call mma_allocate(SGS%Ver,SGS%nVert0,Label='SGS%Ver')
SGS%Ver(:)=1 ! Default to that all vertices are included.
call mma_allocate(CONN,SGS%nVert,Label='CONN')

! KILL VERTICES THAT DO NOT OBEY RESTRICTIONS.
do IV=1,SGS%nVert-1  ! loop over all vertices but the tail vertex
  L = SGS%DRT0(IV,LTAB)       ! Pick up the level index
  N = SGS%DRT0(IV,NTAB)       ! pick up the number of particles expressed by the vertex
  ! Test if the number of partciles is too fe
  if (N < Lim(L)) SGS%Ver(IV) = 0
end do

NCHANGES = 1 ! Initiate first loop
do while (NCHANGES > 0)
  ! REMOVE ARCS HAVING A DEAD UPPER OR LOWER VERTEX.
  ! COUNT THE NUMBER OF ARCS REMOVED OR VERTICES KILLED.
  NCHANGES = 0
  ! Loop over all vertices but the tail/root vertex.
  do IV=1,SGS%nVert-1
    if (SGS%Ver(IV) == 0) then
      ! Vertex IV is eliminated
      ! Loop over the step vectors
      do IC=0,3
        ID = SGS%Down0(IV,IC)  ! Pick up the index of the lower level vertex which it connects to with this step
        ! If it is a viable node/vertex disable the arc
        ! If ID > 0 that means that vertices (IV,ID) are connected by an arc described by the step vector IC
         if (ID > 0) then
          SGS%Down0(IV,IC) = 0
          NCHANGES = NCHANGES+1
        end if
      end do
    else
      !Vertex IV is NOT eliminated
      NLD = 0
      ! Loop over all step vectors
      do IC=0,3
        ID = SGS%Down0(IV,IC)  ! Pick up the node it is connected to
        ! If vertex IV is connected to vertex ID through the step IC then do
        if (ID > 0) then
          ! Check if vertex ID has been removed
          if (SGS%Ver(ID) == 0) then
            ! If ID has been removed change the arc to be void
            SGS%Down0(IV,IC) = 0
            NCHANGES = NCHANGES+1
          else
            NLD = NLD+1
          end if
        end if
      end do
      ! If NLD==0 that means that vertex IV is not connected to any vertex below.
      if (NLD == 0) then
        SGS%Ver(IV) = 0   ! Remove vertex IV
        NCHANGES = NCHANGES+1
      end if
    end if
  end do
  ! ALSO CHECK ON CONNECTIONS FROM ABOVE:
  CONN(:) = 0
  CONN(1) = SGS%Ver(1)
  ! Loop over all vertex but the last one -- the tail/root vertex
  do IV=1,SGS%nVert-1
    ! Process if the vertex is viable
    if (SGS%Ver(IV) == 1) then
      ! Loop over all step vectors
      do IC=0,3
        ID = SGS%Down0(IV,IC) ! Pick up the lower level vertex this step vector connects to
        Test = ID > 0     ! Test=.TRUE. if the arc is viable
        ! If the arc is viable check that the vertex is not removed
        if (Test) Test = SGS%Ver(ID) == 1
        ! if the arc is viable and the vertex is not removed mark vertex ID up as a vertex that is
        ! connected to an vertex at the level above.
        if (Test) CONN(ID) = 1
      end do
    end if
  end do
  ! Loop over all vertice
  do IV=1,SGS%nVert
    ! Remove any vertex that is kept but doesn't connect with any vertex above.
    if ((SGS%Ver(IV) == 1) .and. (CONN(IV) == 0)) then
      SGS%Ver(IV) = 0
      NCHANGES = NCHANGES+1
    end if
  end do

end do

! IF NO CHANGES, THE REMAINING GRAPH IS VALID.
! EVERY VERTEX OBEYS THE RESTRICTIONS. EVERY VERTEX IS
! CONNECTED ABOVE AND BELOW (EXCEPTING THE TOP AND BOTTOM)
! TO OTHER CONFORMING VERTICES.
! THE PROCEDURE IS GUARANTEED TO FIND A STABLE SOLUTIONS,
! SINCE EACH ITERATION REMOVES ARCS AND/OR VERTICES FROM THE
! FINITE NUMBER WE STARTED WITH.

if (SGS%Ver(1) == 0) then
  write(u6,*) 'RASSI/RMVERT: Too severe restrictions.'
  write(u6,*) 'Not one single configuration is left.'
  call ABEND()
end if

! Now get the new indexation of the remaining vertices.
NV = 0
do IV=1,SGS%nVert
  if (SGS%Ver(IV) == 1) then
    NV = NV+1
    SGS%Ver(IV) = NV
  end if
end do
SGS%nVert = NV  ! Update nVert

call mma_deallocate(CONN)
call mma_deallocate(Lim)

end subroutine RMVERT

end subroutine MKSGUGA

subroutine SG_Init(nSym,nActEl,iSpin,SGS,CIS,TRS,                  &
                   nRas,nRasEl,nRsPrt,                             &
                   EXS,xLevel,xL2Act,xNLEV,xNSM)

  integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin
  type(SGStruct), intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  type(TRStruct), intent(inout) :: TRS
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt),nRasEl(nRsPrt)
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xnLev, xNSM(MxLev)
  type(EXStruct), optional, intent(inout) :: EXS

  call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,TRS,    &
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
!   Call MkTrans(SGS,CIS,EXS,TRS)
    Call MkTrans(SGS,CIS,TRS)

    ! Count coupling coefficients in compressed blocks indexed by excitation/operator type, symmetry and midvertex.
    call MkNRCOUP(SGS,CIS,EXS)

    ! Explicitly generates the compressed coupling tuples '(left walk, right walk, value index)' and compacts repeated
    ! numerical values into 'VTab'.
    call MKCOUP(SGS,CIS,EXS)
  end if

end subroutine SG_Init

subroutine SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,TRS,          &
                          nRas,nRasEl,nRsPrt,                     &
                          EXS,xLevel,xL2Act,xNLEV,xNSM,Do_MkSGUGA)

  integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin
  type(SGStruct), intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  type(TRStruct), intent(inout) :: TRS
  integer(kind=iwp), intent(in) :: nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt)
  type(EXStruct), optional, intent(inout) :: EXS
  integer(kind=iwp), optional, intent(in) :: xLevel(MxLev), xL2Act(MxLev), &
                                             xNLEV, xNSM(MxLev)
  logical(kind=iwp), optional, intent(in) :: Do_MkSGUGA
  integer(kind=iwp) :: iSym

Select case(nRsPrt)
Case(1)
   SGS%IFRAS=0
Case(3)
   SGS%IFRAS=1
End Select
Do iSym = 1, nSym
   If (Sum(nRas(iSym,1:nRsPrt))/=0) SGS%IFRAS=SGS%IFRAS+1
End Do
SGS%nRasEL(1:nRsPrt)=nRasEl(1:nRsPrt)
SGS%nRas(1:MxSym,1:nRsPrt)=nRas(1:MxSym,1:nRsPrt)
SGS%nRsPrt=nRsPrt

! Make sure that we start from a clean slate.
  if (present(EXS)) then
   ! Here if the extended parameter list was used.
    call SG_Free(SGS,CIS,TRS,EXS)
  else
   ! Here if the terse parameter list was used.
    call SG_Free(SGS,CIS,TRS)
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

subroutine SG_Free(SGS,CIS,TRS,EXS)
! PURPOSE: FREE THE SGUGA TABLES

type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
type(TRStruct), intent(inout) :: TRS
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

call Trans_Free(TRS)

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
call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr',safe='*')  ! First index referenced by IVERT, SYM, and ISTEP
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
      SELECT CASE(ISTP)
!       CASE(0,3)
!         ISML = 1
        CASE(1,2)
          ISML = SGS%ISm(LEV)
        CASE DEFAULT
          ISML = 1
      END SELECT

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

end do ! IHALF

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
    SELECT CASE(ITT)
!     Case(0,3)
!       IVRT = IVLT               ! Upper right vertex is the same as the upper left vertex.
      Case(1,2)
        IVRT = CIS%IVR(IVLT,ITT)  ! Pick up the associated upper right vertex, depends on delta(b)
      Case Default
        IVRT = IVLT               ! Upper right vertex is the same as the upper left vertex.
    End Select
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
      case default
        V = Zero ! Dummy assignment
        call Abend()
    end select
    CIS%VSGM(IVLT,ISGT) = V   ! Store away the segment value.
  end do
end do

end subroutine MKSEG

subroutine MkNRCOUP(SGS,CIS,EXS)

type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IBSYM, ICL, INDEO, INDEOB, INDEOT, IP, IPQ, IQ, ISGT, ISYDS1, ISYM, ISYUS1, ITSYM, IVLB, IVLT, LEV, LFTSYM, &
                     MV, MV1, MV2, MV3, MV4, MV5, MXDWN, MXUP, N, NDWNS1, NSGMX, NSGTMP, NT1TMP, NT2TMP, NT3TMP, NT4TMP, NT5TMP, &
                     NUPS1, INDEO0
integer(kind=iwp), allocatable :: NRL(:,:,:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IS, IST, NCP, NUW
#endif
call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW',safe='*')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW',safe='*')

call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF',safe='*')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF',safe='*')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF',safe='*')

EXS%MxEO = SGS%nLev + SGS%nLev + (SGS%nLev*(SGS%nLev+1))/2

call mma_allocate(EXS%NOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
call mma_allocate(EXS%IOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')

! NRL(sym,vertex,indeo): number of valid segment paths of a given symmetry and operators class arriving at a given vertex
call mma_allocate(NRL,[1,SGS%nSym],[1,SGS%nVert],[0,EXS%MxEO],Label='NRL')
! indeo=0 denotes a ordinary half-walk with no open- or closed-loop attached
INDEO0=0

! For upper walks
NRL(:,1:SGS%MVEnd,:) = 0
NRL(1,1,0) = 1     ! Initiate that for symmetry=1, the root vertex and indeo=0 that the number of valid setment paths are 1

do IVLT=1,SGS%MVSta-1                     ! Loop from head vertex to last vertex before the midvertices, IVLT
  LEV = SGS%DRT(IVLT,LTAB)                ! Pick up the level of the vertex, LEV
  do ISGT=1,nSeg                          ! loop over all segments, ISGT
    IVLB = CIS%ISGM(IVLT,ISGT)            ! Pick up the lower-left destination vertex, IVLB
    if (IVLB == 0) cycle                  ! Cycle if not a valid vertex
    ICL = IC1(ISGT)                       ! Pick up the bra case

    Select case (ICL)                     ! Pick up the symmetry label of the segment
      CASE(1,2)
        ISYM = SGS%ISm(LEV)
      CASE DEFAULT
        ISYM = 1
    End Select

    do ITSYM=1,SGS%nSym                   ! Loop over all symmetries
      IBSYM = Mul(ITSYM,ISYM)             ! If the walk was of symmetry ITSYM at the level above this segment will render it now
                                          ! to be of symmetry IBSYM

      ! Accumulate contributions to vertices below
      select case (ISGT)
        case (1:4)
          ! THIS IS AN UPPER WALK.
          NRL(IBSYM,IVLB,INDEO0) = NRL(IBSYM,IVLB,INDEO0)+NRL(ITSYM,IVLT,INDEO0)
        case (5:8)
          ! THIS IS AN TOP SEGMENT.
          INDEOB= LEV+(IBVPT(ISGT)-1)*SGS%nLev
          NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEO0)
        case (9:18)
          ! THIS IS A MID-SEGMENT.
          do IP=LEV+1,SGS%nLev
            INDEOT = IP+(ITVPT(ISGT)-1)*SGS%nLev
            INDEOB = IP+(IBVPT(ISGT)-1)*SGS%nLev
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
          end do
        case (19:22)
          ! THIS IS A BOTTOM SEGMENT.
          do IP=LEV+1,SGS%nLev
            INDEOT = IP+(ITVPT(ISGT)-1)*SGS%nLev
            IPQ = (IP*(IP-1))/2+LEV
            INDEOB = 2*SGS%nLev + IPQ
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
          end do
        case default
          ! THIS IS A LOWER WALK.
          do INDEO=2*SGS%nLev+1,EXS%MxEO
            NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,INDEO)
          end do
      end select
    end do
  end do
end do

! store the accumulated number of upper walks
MXUP = 0
do MV=1,CIS%nMidV                  ! loop over midverticies
  IVLT = MV+SGS%MVSta-1            ! Get the absolute vertex index
  do LFTSYM=1,SGS%nSym             ! Loop over symmetries
    CIS%NOW(1,LFTSYM,MV) = NRL(LFTSYM,IVLT,0)
    MXUP = max(MXUP,CIS%NOW(1,LFTSYM,MV))
    do INDEO=1,EXS%MxEO
      EXS%NOCP(INDEO,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO)
    end do
  end do
end do

! For lower walks
NRL(:,SGS%MVSta:SGS%nVert,:) = 0
NRL(1,SGS%nVert,0) = 1

do IVLT=SGS%nVert-1,SGS%MVSta,-1
  LEV = SGS%DRT(IVLT,LTAB)
  do ISGT=1,nSeg
    IVLB = CIS%ISGM(IVLT,ISGT)
    if (IVLB == 0) cycle
    ICL = IC1(ISGT)
    ISYM = 1
    if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
    do ITSYM=1,SGS%nSym
      IBSYM = Mul(ITSYM,ISYM)
      select case (ISGT)
        case (23:)
          ! THIS IS A LOWER WALK.
          NRL(ITSYM,IVLT,0) = NRL(ITSYM,IVLT,0)+NRL(IBSYM,IVLB,0)
        case (19:22)
          ! THIS IS AN BOTTOM SEGMENT.
          INDEO = LEV+(ITVPT(ISGT)-1)*SGS%nLev
          NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,0)
        case (9:18)
          ! THIS IS A MID-SEGMENT.
          do IQ=1,LEV-1
            INDEOT = IQ+(ITVPT(ISGT)-1)*SGS%nLev
            INDEOB = IQ+(IBVPT(ISGT)-1)*SGS%nLev
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
          end do
        case (5:8)
          ! THIS IS AN TOP SEGMENT.
          do IQ=1,LEV-1
            INDEOB = IQ+(IBVPT(ISGT)-1)*SGS%nLev
            IPQ = (LEV*(LEV-1))/2+IQ
            INDEOT = IPQ+2*SGS%nLev
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
          end do
        case default
          ! THIS IS AN UPPER WALK.
          do IPQ=1,(LEV*(LEV-1))/2
            INDEO = 2*SGS%nLev+IPQ
            NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,INDEO)
          end do
      end select
    end do
  end do
end do

MXDWN = 0
do MV=1,CIS%nMidV
  IVLT = MV+SGS%MVSta-1
  do LFTSYM=1,SGS%nSym
    CIS%NOW(2,LFTSYM,MV) = NRL(LFTSYM,IVLT,0)
    MXDWN = max(MXDWN,CIS%NOW(2,LFTSYM,MV))
    do INDEO=1,EXS%MxEO
      N = NRL(LFTSYM,IVLT,INDEO)
      if (N /= 0) EXS%NOCP(INDEO,LFTSYM,MV) = N
    end do
  end do
end do

EXS%nICOup = 0
do INDEO=1,EXS%MxEO
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      EXS%IOCP(INDEO,LFTSYM,MV) = EXS%nICOup
      EXS%nICOup = EXS%nICOup +  EXS%NOCP(INDEO,LFTSYM,MV)
    end do
  end do
end do

Call Mk_IOW(CIS,SGS)
call CSFCOUNT(CIS,SGS)

!AR INSERT FOR US IN SIGMA ROUTINE

NSGMX = 1
NDWNS1 = 0  ! Dummy initialization
NSGTMP = max(MXUP,MXDWN)
do MV3=1,CIS%nMidV
  MV1 = EXS%MVL(MV3,2)
  MV2 = EXS%MVL(MV3,1)
  MV4 = EXS%MVR(MV3,1)
  MV5 = EXS%MVR(MV3,2)
  do ISYUS1=1,SGS%nSym
    NUPS1 = CIS%NOW(1,ISYUS1,MV3)
    do ISYDS1=1,SGS%nSym
      NDWNS1 = CIS%NOW(2,ISYDS1,MV3)
      NSGMX = max(NSGMX,CIS%NOCSF(ISYUS1,MV3,ISYDS1))

      if (MV1 /= 0) then
        NT4TMP = NUPS1*CIS%NOW(2,ISYDS1,MV1)
        NSGTMP = max(NSGTMP,NT4TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV1)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV2 /= 0) then
        NT3TMP = NUPS1*CIS%NOW(2,ISYDS1,MV2)
        NSGTMP = max(NSGTMP,NT3TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV2)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV4 /= 0) then
        NT1TMP = NUPS1*CIS%NOW(2,ISYDS1,MV4)
        NSGTMP = max(NSGTMP,NT1TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV4)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV5 /= 0) then
        NT2TMP = NUPS1*CIS%NOW(2,ISYDS1,MV5)
        NSGTMP = max(NSGTMP,NT2TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV5)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

    end do
  end do
end do

call mma_allocate(EXS%SGTMP,NSGTMP,Label='EXS%SGTMP')

#ifdef _DEBUGPRINT_
write(u6,600) MXUP,MXDWN,NSGMX,NSGMX,NSGTMP

!AR END OF INSERT
NUW=CIS%NUW
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*) ' TOTAL NR OF COUPL COEFFS ',EXS%nICOup
INDEO = 2*SGS%nLev+1
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
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
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

subroutine MKCOUP(SGS,CIS,EXS)
! Purpose: Compute and return the table ICOUP(1..3,ICOP).
! The number of coupling coeffs is obtained from NOCP, the offset to
! the ICOP numbering is given by IOCP. The numbers ICOUP(1..3,ICOP) are
! then the ket and bra half-walks, enumerated by the Lund scheme,
! and the index into the VTAB table (the pool of possible values of
! coupling coefficients).

! Any loop is regarded as a segment path from top to midlevel, or
! from midlevel to bottom.
! The segment path is described by the table ISGPTH. It is
! essentially a list of which one of segments nr 1..nSeg that are
! used at each level. The segments are of three types:
! Type 0: Upwalk segment or top loop segment.
! Type 1: Mid segment.
! Type 2: Bottom segment.
! Type 3: Downwalk segment.
! ISGPTH(IVLFT,LEV)=Left upper vertex.
! ISGPTH(ITYPE,LEV)=Type of segment, (0..3).
! ISGPTH(IAWSL,LEV)=Left arc weight sum (from top, or from bottom).
! ISGPTH(IAWSR,LEV)=Similar, right.
! ISGPTH(ILS  ,LEV)=Left symmetry label (accum from top or bottom).
! ISGPTH(ICS  ,LEV)=Left coupling case number (0..3).
! ISGPTH(ISEG ,LEV)=Segment type (1..nSeg).
! These indices are used to denote the columns of table ISGPTH.

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: i, i1, i2, IAWS, IC, ICL, ICOP, ICR, IHALF, iLnd, IndEO, iP, iPos, iQ, iSg, iSgt, iSym, iT, iTyp, &
                     iTypMx, iTypT, iVlb, iVlt, iVrt, iVrTop, iVTab, iVTEnd, iVTop, iVTSta, L, Lev, Lev1, Lev2, LftSym, LL, MV, &
                     nCheck, nVTab_Final
real(kind=wp) :: C
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I3, ICOP1, ICOP2, ICP1, ICP2, N, NRC, NRCPQ
real(kind=wp) :: CP
#endif
integer(kind=iwp), allocatable :: ILNDW(:), ISGPTH(:,:)
real(kind=wp), allocatable :: val(:), VTab(:)
integer(kind=iwp), parameter :: IVLFT = 1, ITYPE = 2, IAWSL = 3, IAWSR = 4, ILS = 5, ICS = 6, ISEG = 7, nVTab = 5000

call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')

!  Special case
If (SGS%nLev==1) Then
   nVTab_Final=0
   call mma_allocate(EXS%VTab,nVTab_Final,Label='EXS%VTab')
   Return
End If

! NOW IS ZEROED AND WILL BE USED AS AN ARRAY OF COUNTERS, BUT WILL BE RESTORED FINALLY.
CIS%NOW(:,:,:) = 0
! SIMILAR FOR THE COUPLING COEFFICIENT TABLE:
EXS%NOCP(:,:,:) = 0

call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
call mma_allocate(ISGPTH,[1,7],[0,SGS%nLev],Label='ISGPTH')
call mma_allocate(val,[0,SGS%nLev],Label='val')
call mma_allocate(VTab,nVTab,Label='VTab')

! COUPLING COEFFICIENT VALUE TABLE:
NVTAB_FINAL = 2
VTab(1) =  One
VTab(2) = -One

NCHECK = 0

do IHALF=1,2
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = SGS%nLev
    LEV2 = SGS%MidLev
    ITYPMX = 0
  else
    IVTSTA = SGS%MVSta
    IVTEND = SGS%MVEnd
    LEV1 = SGS%MidLev
    LEV2 = 0
    ITYPMX = 2
  end if

  do IVTOP=IVTSTA,IVTEND   ! just the head vertex, or the mid vertices
    do ITYP=0,ITYPMX       ! loop over all possible segments types:
                           ! for the top level we have only head segments.
      IVRTOP = IVTOP
      if (ITYP > 0) IVRTOP = CIS%IVR(IVTOP,ITYP)
      if (IVRTOP == 0) cycle
      LEV = LEV1

      ! Initialize the temporary table
      ISGPTH(IVLFT,LEV) = IVTOP    ! Left vertex index
      ISGPTH(ITYPE,LEV) = ITYP     ! segment type vertex
      ISGPTH(IAWSL,LEV) = 0
      ISGPTH(IAWSR,LEV) = 0
      ISGPTH(ILS,  LEV) = 1
      ISGPTH(ISEG, LEV) = 0
      Val(LEV)          = One

!     Loop from the head level to the mid level, or from the mid level to the root
      do while (LEV <= LEV1)
        ITYPT = ISGPTH(ITYPE,LEV)
        IVLT = ISGPTH(IVLFT,LEV)

        ! Find the next valid segment of this level and explore it.
        do ISGT=ISGPTH(ISEG,LEV)+1,nSeg
          IVLB = CIS%ISGM(IVLT,ISGT)
          if (IVLB == 0) cycle
          if (ITYPT == ITVPT(ISGT)) exit
        end do

        if (ISGT > nSeg) then   ! All possible segments investigated at this level
          ISGPTH(ISEG,LEV) = 0  ! Reset
          LEV = LEV+1           ! back up a level
        else
          ISGPTH(ISEG,LEV) = ISGT
          ICL = IC1(ISGT)
          ISYM = 1
          if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
          IVRT = IVLT
          if ((ITYPT == 1) .or. (ITYPT == 2)) IVRT = CIS%IVR(IVLT,ITYPT)
          ICR = IC2(ISGT)
          ISGPTH(ICS,LEV) = ICL
          LEV = LEV-1
          ISGPTH(IAWSL,LEV) = ISGPTH(IAWSL,LEV+1)+SGS%MAW(IVLT,ICL)
          ISGPTH(IAWSR,LEV) = ISGPTH(IAWSR,LEV+1)+SGS%MAW(IVRT,ICR)
          val(LEV) = val(LEV+1)*CIS%VSGM(IVLT,ISGT)
          ISGPTH(ILS,LEV) = Mul(ISYM,ISGPTH(ILS,LEV+1))
          ISGPTH(IVLFT,LEV) = IVLB
          ISGPTH(ITYPE,LEV) = IBVPT(ISGT)
          ISGPTH(ISEG,LEV) = 0
          if (LEV > LEV2) cycle

          MV = ISGPTH(IVLFT,SGS%MidLev)+1-SGS%MVSta
          LFTSYM = ISGPTH(ILS,LEV2)
          IT = ISGPTH(ITYPE,SGS%MidLev)
          if (IT == 0) IT = 3
          if (ISGPTH(ITYPE,LEV2) == 0) IT = 0

          if (IT == 0) then
            ILND = 1+CIS%NOW(IHALF,LFTSYM,MV)
            IAWS = ISGPTH(IAWSL,LEV2)
            ILNDW(IAWS) = ILND
            CIS%NOW(IHALF,LFTSYM,MV) = ILND
            IPOS = CIS%IOW(IHALF,LFTSYM,MV)+(ILND-1)*CIS%nIpWlk
            do LL=LEV2+1,LEV1,nPack
              IC = 0
              do L=min(LL+nPack-1,LEV1),LL,-1
                IC = 4*IC+ISGPTH(ICS,L)
              end do
              IPOS = IPOS+1
              CIS%ICase(IPOS) = IC
            end do
          else
            IP = 0
            IQ = 0
            do L=LEV2+1,LEV1
              ISG = ISGPTH(ISEG,L)
              if ((ISG >= 5) .and. (ISG <= 8)) IP = L
              if ((ISG >= 19) .and. (ISG <= 22)) IQ = L
            end do
            if (IP == 0) IP = IQ
            INDEO = SGS%nLev*(IT-1)+IP
            if (IT == 3) INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
            ICOP = 1+EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
            ICOP = EXS%IOCP(INDEO,LFTSYM,MV)+ICOP
            NCHECK = NCHECK+1
            if (ICOP > EXS%nICoup) then
              write(u6,*) ' ERROR: NICOUP=',EXS%nICoup
              write(u6,*) ' NR OF COUPS PRODUCED:',NCHECK
              write(u6,*) '           TYPE NR IT:',IT
              write(u6,*) '            IP,IQ    :',IP,IQ
              write(u6,*) '            INDEO    :',INDEO
              write(u6,*) '        MIDVERTEX MV :',MV
              write(u6,*) ' LEFT SYMMETRY LFTSYM:',LFTSYM
              write(u6,*) ' COUP OFFSET IOCP    :',EXS%IOCP(INDEO,LFTSYM,MV)
              write(u6,*) ' COUP SERIAL NR ICOP :',ICOP
              write(u6,*) ' D:O, WITHOUT OFFSET :',ICOP-EXS%IOCP(INDEO,LFTSYM,MV)
              write(u6,*) ' CURRENT NOCP NUMBER :',EXS%NOCP(INDEO,LFTSYM,MV)
              call ABEND()
            end if

            C = val(LEV2)
            do I=1,NVTAB_FINAL
              IVTAB = I
              if (abs(C-VTab(I)) < 1.0e-10_wp) exit
            end do
            if (I > NVTAB_FINAL) then
              NVTAB_FINAL = NVTAB_FINAL+1
              if (NVTAB_FINAL > nVTab) then
                write(u6,*) 'MKCOUP: NVTAB_FINAL=',NVTAB_FINAL
                write(u6,*) 'NVTAB_FINAL should not be allowed to grow'
                write(u6,*) 'beyond nVTab which was set provisionally'
                write(u6,*) 'in subroutine GINIT in file ginit.f.'
                write(u6,*) 'Now nVTab=',nVTab
                write(u6,*) 'This may indicate a problem with your input.'
                write(u6,*) 'If you do want to do this big calculation, try'
                write(u6,*) 'increasing nVTab in GINIT and recompile.'
                call ABEND()
              end if
              VTab(NVTAB_FINAL) = C
              IVTAB = NVTAB_FINAL
            end if
            EXS%ICoup(1,ICOP) = ISGPTH(IAWSL,LEV2)
            EXS%ICoup(2,ICOP) = ISGPTH(IAWSR,LEV2)
            EXS%ICoup(3,ICOP) = IVTAB
            if (ICOP > EXS%nICoup) then
              write(u6,*) 'MKCOUP: ICOP>NICOUP!'
              call ABEND()
            end if
          end if

          LEV = LEV+1
        end if
      end do

    end do
  end do
end do
! RENUMBER THE COUPLING COEFFICIENT INDICES BY LUND SCHEME:
do ICOP=1,EXS%nICoup
  I1 = EXS%ICoup(1,ICOP)
  I2 = EXS%ICoup(2,ICOP)
  EXS%ICoup(1,ICOP) = ILNDW(I1)
  EXS%ICoup(2,ICOP) = ILNDW(I2)
end do

call mma_deallocate(val)
call mma_deallocate(ISGPTH)
call mma_deallocate(ILNDW)

#ifdef _DEBUGPRINT_
ICOP1 = 0
ICOP2 = 0
write(u6,*) ' NR OF DIFFERENT VALUES OF COUP:',NVTAB_FINAL
do ICOP=1,EXS%nICoup
  I3 = EXS%ICoup(3,ICOP)
  if (I3 == 1) ICOP1 = ICOP1+1
  if (I3 == 2) ICOP2 = ICOP2+1
end do
write(u6,*)
write(u6,*) ' NR OF COUPS WITH VALUE  1.0:',ICOP1
write(u6,*) ' NR OF COUPS WITH VALUE -1.0:',ICOP2
write(u6,*)
write(u6,*) ' COUPLING COEFFICIENTS:'
write(u6,*) '    IP    IQ    MV LFTSYM NOCP'
write(u6,*) ' 1. OPEN LOOPS TYPE 1.'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(IP,LFTSYM,MV)
      write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
    end do
  end do
end do
write(u6,*) ' 2. OPEN LOOPS TYPE 2.'
do IP=1,SGS%nLev
  INDEO = SGS%nLev+IP
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(INDEO,LFTSYM,MV)
      write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
    end do
  end do
end do
write(u6,*) ' 3. CLOSED LOOPS.'
do IP=1,SGS%nLev
  do IQ=1,IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do LFTSYM=1,SGS%nSym
        N = EXS%NOCP(INDEO,LFTSYM,MV)
        write(u6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,N
      end do
    end do
  end do
end do
write(u6,*)
write(u6,*) ' COUPLING COEFFICIENTS:'
write(u6,*) '    IP    IQ    MV LFTSYM ICOP ICOUP1&2   COUP'
write(u6,*) ' 1. OPEN LOOPS TYPE 1.'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(IP,LFTSYM,MV)
      ICOP = EXS%IOCP(IP,LFTSYM,MV)
      do I=1,N
        ICOP = ICOP+1
        ICP1 = EXS%ICoup(1,ICOP)
        ICP2 = EXS%ICoup(2,ICOP)
        CP = VTab(EXS%ICoup(3,ICOP))
        write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
      end do
    end do
  end do
end do
write(u6,*) ' 2. OPEN LOOPS TYPE 2.'
do IP=1,SGS%nLev
  INDEO = SGS%nLev+IP
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(INDEO,LFTSYM,MV)
      ICOP = EXS%IOCP(INDEO,LFTSYM,MV)
      do I=1,N
        ICOP = ICOP+1
        ICP1 = EXS%ICoup(1,ICOP)
        ICP2 = EXS%ICoup(2,ICOP)
        CP = VTab(EXS%ICoup(3,ICOP))
        write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
      end do
    end do
  end do
end do
write(u6,*) ' 3. CLOSED LOOPS.'
do IP=1,SGS%nLev
  do IQ=1,IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do LFTSYM=1,SGS%nSym
        N = EXS%NOCP(INDEO,LFTSYM,MV)
        ICOP = EXS%IOCP(INDEO,LFTSYM,MV)
        do I=1,N
          ICOP = ICOP+1
          ICP1 = EXS%ICoup(1,ICOP)
          ICP2 = EXS%ICoup(2,ICOP)
          CP = VTab(EXS%ICoup(3,ICOP))
          write(u6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,ICOP,ICP1,ICP2,CP
        end do
      end do
    end do
  end do
end do
write(u6,*)
write(u6,*) ' CONVENTIONAL NR OF COUPLING COEFFS, BY PAIR:'
NRC = 0
do IP=2,SGS%MidLev
  do IQ=1,IP-1
    NRCPQ = 0
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(INDEO,LFTSYM,MV)*CIS%NOW(1,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
do IP=SGS%MidLev+1,SGS%nLev
  do IQ=1,SGS%MidLev
    NRCPQ = 0
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(IP,LFTSYM,MV)*EXS%NOCP(IQ,ISYM,MV)
        NRCPQ = NRCPQ+EXS%NOCP(SGS%nLev+IP,LFTSYM,MV)*EXS%NOCP(SGS%nLev+IQ,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
do IP=SGS%MidLev+2,SGS%nLev
  do IQ=SGS%MidLev+1,IP-1
    NRCPQ = 0
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(INDEO,LFTSYM,MV)*CIS%NOW(2,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
write(u6,*)
write(u6,*) ' TOTAL CONVENTIONAL COUPLING COEFFS:',NRC
#endif

call mma_allocate(EXS%VTab,nVTab_Final,Label='EXS%VTab')
EXS%VTab(1:nVTab_final) = VTab(1:nVTab_final)
call mma_deallocate(VTab)

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

  select case (ICLASS)
  case (1,2)
    NeedsRightPartner = .true.
  case default
    NeedsRightPartner = .false.
  end select
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

!subroutine MKTRANS(SGS,CIS,EXS,TRS)
subroutine MKTRANS(SGS,CIS,TRS)
  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
! type(EXStruct), intent(in)    :: EXS
  type(TRStruct), intent(inout) :: TRS

  integer(kind=iwp) :: IVLT, ISGT, ICLASS
  integer(kind=iwp) :: LEV, IVLB
  integer(kind=iwp) :: ICL, ICR, ITOP, IBOT
  integer(kind=iwp) :: IPRT, ISYM, IOBAND, IFLAG
  integer(kind=iwp) :: ITR, N
  integer(kind=iwp), allocatable :: IPOS(:,:)
  real(kind=wp) :: VSEG

  ! Initialize metadata
! call TRANS_Free(TRS) this should go out.

  TRS%nClass    = 4
  TRS%nOpenBand = 2

  ! Allocate bucket tables
  call mma_deallocate(TRS%NTR,safe='*')
  call mma_deallocate(TRS%ITR0,safe='*')
  Write (6,*) 'Allocate TRS%NTR'
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

!#ifdef _VERIFY_
  write(u6,*) 'MKTRANS: total transitions = ', TRS%nTrans
  write(u6,*) 'MKTRANS: sum NTR           = ', sum(TRS%NTR)
  ! Validation of the transition table
  If (TRS%nTrans/=sum(TRS%NTR)) Then
     Write (u6,*) 'TRS%nTrans/=sum(TRS%NTR)'
     Call Abend()
  End If
  call CheckTransitionBucketCounts(SGS,CIS,TRS)

  if (sum(TRS%NTR) /= CountAllTransitions(SGS,CIS,TRS%nClass)) then
    write(u6,*) 'MKTRANS ERROR: total transition count mismatch.'
    call Abend()
  end if
!#endif

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

      ITOP = ITVPT(ISGT)
      IBOT = IBVPT(ISGT)

      ICL = IC1(ISGT)
      ICR = IC2(ISGT)

      ISYM = 1
      if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)

      IPRT   = PartnerSlot(ITOP)
      IOBAND = OpenBand(ITOP)
      VSEG   = CIS%VSGM(IVLT,ISGT)

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

!#ifdef _VERIFY_

  call DumpAndCheckAllTransitionBuckets(SGS,CIS,TRS)

contains


subroutine DumpAndCheckAllTransitionBuckets(SGS,CIS,TRS)
  implicit none

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(TRStruct), intent(in) :: TRS

  integer(kind=iwp) :: IVLT, ICLASS, ISGT, LEV
  integer(kind=iwp) :: NDIR, NTRB, IT0, ITR, K
  integer(kind=iwp) :: IVLB, ITOP, IBOT, ICL, ICR
  integer(kind=iwp) :: IPRT, ISYM, IOBAND
  real(kind=wp)     :: VSEG
  logical           :: mismatch
  real(kind=wp), parameter :: TOL = 1.0e-12_wp

  ! Directly reconstructed bucket contents from old logic
  integer(kind=iwp), allocatable :: D_ISGT(:), D_IVLB(:)
  integer(kind=iwp), allocatable :: D_ITOP(:), D_IBOT(:)
  integer(kind=iwp), allocatable :: D_ICL(:), D_ICR(:)
  integer(kind=iwp), allocatable :: D_IPRT(:), D_ISYM(:), D_IOBAND(:)
  real(kind=wp),    allocatable :: D_VSEG(:)

  ! Allocate enough room for the worst case:
  ! at most nSeg entries in one bucket
  allocate(D_ISGT(nSeg), D_IVLB(nSeg))
  allocate(D_ITOP(nSeg), D_IBOT(nSeg))
  allocate(D_ICL(nSeg),  D_ICR(nSeg))
  allocate(D_IPRT(nSeg), D_ISYM(nSeg), D_IOBAND(nSeg))
  allocate(D_VSEG(nSeg))

  write(u6,*) ' '
  write(u6,*) '==========================================================='
  write(u6,*) 'MKTRANS: dump and full bucket comparison against old logic'
  write(u6,*) '==========================================================='
  write(u6,*) ' '

  do ICLASS = 0, TRS%nClass-1

    write(u6,'(A,I4)') '--- INPUT CLASS = ', ICLASS

    do IVLT = 1, SGS%nVert

      LEV = SGS%DRT(IVLT,LTAB)

      ! Reconstruct the old/direct bucket contents for (IVLT,ICLASS)
      NDIR = 0

      do ISGT = 1, nSeg

        IVLB = CIS%ISGM(IVLT,ISGT)
        if (IVLB == 0) cycle

        ITOP = ITVPT(ISGT)
        if (ITOP /= ICLASS) cycle

        IBOT = IBVPT(ISGT)
        ICL  = IC1(ISGT)
        ICR  = IC2(ISGT)

        ISYM = 1
        if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)

        IPRT   = PartnerSlot(ITOP)
        IOBAND = OpenBand(ITOP)
        VSEG   = CIS%VSGM(IVLT,ISGT)

        NDIR = NDIR + 1

        D_ISGT(NDIR)   = ISGT
        D_IVLB(NDIR)   = IVLB
        D_ITOP(NDIR)   = ITOP
        D_IBOT(NDIR)   = IBOT
        D_ICL(NDIR)    = ICL
        D_ICR(NDIR)    = ICR
        D_IPRT(NDIR)   = IPRT
        D_ISYM(NDIR)   = ISYM
        D_IOBAND(NDIR) = IOBAND
        D_VSEG(NDIR)   = VSEG

      end do

      ! Pull the new bucket from TRS
      NTRB = TRS%NTR(IVLT,ICLASS)
      IT0  = TRS%ITR0(IVLT,ICLASS)

      ! Dump the bucket contents
      write(u6,'(A,I6,A,I4,A,I4)') 'IVLT=',IVLT, ' direct=',NDIR, ' trs=',NTRB

      if (NDIR == 0 .and. NTRB == 0) cycle

      write(u6,*) '  Direct reconstruction:'
      if (NDIR == 0) then
        write(u6,*) '    [empty]'
      else
        do K = 1, NDIR
          write(u6,*) '    k=',K, &
                      ' ISGT=',D_ISGT(K), &
                      ' IVLB=',D_IVLB(K), &
                      ' ITOP=',D_ITOP(K), &
                      ' IBOT=',D_IBOT(K), &
                      ' ICL=', D_ICL(K), &
                      ' ICR=', D_ICR(K), &
                      ' IPRT=',D_IPRT(K), &
                      ' ISYM=',D_ISYM(K), &
                      ' IOB=', D_IOBAND(K), &
                      ' VSEG=',D_VSEG(K)
        end do
      end if

      write(u6,*) '  TRS bucket:'
      if (NTRB == 0) then
        write(u6,*) '    [empty]'
      else
        do K = 1, NTRB
          ITR = IT0 + K
          write(u6,*) '    k=',K, &
                      ' ISGT=',TRS%ISGT(ITR), &
                      ' IVLB=',TRS%IVLB(ITR), &
                      ' ITOP=',TRS%ITOP(ITR), &
                      ' IBOT=',TRS%IBOT(ITR), &
                      ' ICL=', TRS%ICL(ITR), &
                      ' ICR=', TRS%ICR(ITR), &
                      ' IPRT=',TRS%IPRT(ITR), &
                      ' ISYM=',TRS%ISYM(ITR), &
                      ' IOB=', TRS%IOBAND(ITR), &
                      ' VSEG=',TRS%VSEG(ITR)
        end do
      end if

      ! Compare bucket sizes
      mismatch = .false.

      if (NDIR /= NTRB) then
        write(u6,*) '  >>> MISMATCH: bucket size differs.'
        mismatch = .true.
      end if

      ! Compare entry by entry in order
      if (.not. mismatch) then
        do K = 1, NDIR
          ITR = IT0 + K

          if (TRS%ISGT(ITR)   /= D_ISGT(K))   mismatch = .true.
          if (TRS%IVLB(ITR)   /= D_IVLB(K))   mismatch = .true.
          if (TRS%ITOP(ITR)   /= D_ITOP(K))   mismatch = .true.
          if (TRS%IBOT(ITR)   /= D_IBOT(K))   mismatch = .true.
          if (TRS%ICL(ITR)    /= D_ICL(K))    mismatch = .true.
          if (TRS%ICR(ITR)    /= D_ICR(K))    mismatch = .true.
          if (TRS%IPRT(ITR)   /= D_IPRT(K))   mismatch = .true.
          if (TRS%ISYM(ITR)   /= D_ISYM(K))   mismatch = .true.
          if (TRS%IOBAND(ITR) /= D_IOBAND(K)) mismatch = .true.
          if (abs(TRS%VSEG(ITR) - D_VSEG(K)) > TOL) mismatch = .true.

          if (mismatch) then
            write(u6,*) '  >>> MISMATCH at entry k =', K
            write(u6,*) '      Direct: ', D_ISGT(K), D_IVLB(K), D_ITOP(K), D_IBOT(K), &
                         D_ICL(K), D_ICR(K), D_IPRT(K), D_ISYM(K), D_IOBAND(K), D_VSEG(K)
            write(u6,*) '      TRS   : ', TRS%ISGT(ITR), TRS%IVLB(ITR), TRS%ITOP(ITR), TRS%IBOT(ITR), &
                         TRS%ICL(ITR), TRS%ICR(ITR), TRS%IPRT(ITR), TRS%ISYM(ITR), &
                         TRS%IOBAND(ITR), TRS%VSEG(ITR)
            exit
          end if

        end do
      end if

      if (mismatch) then
        write(u6,*) 'MKTRANS ERROR: bucket dump/check failed.'
        write(u6,*) '  Failing vertex =', IVLT
        write(u6,*) '  Failing class  =', ICLASS
        call Abend()
      else
        write(u6,*) '  bucket OK'
      end if

    end do

    write(u6,*) ' '
  end do

  write(u6,*) 'MKTRANS: all buckets matched old logic.'
  write(u6,*) ' '

  deallocate(D_ISGT, D_IVLB)
  deallocate(D_ITOP, D_IBOT)
  deallocate(D_ICL, D_ICR)
  deallocate(D_IPRT, D_ISYM, D_IOBAND)
  deallocate(D_VSEG)

end subroutine DumpAndCheckAllTransitionBuckets


subroutine CheckTransitionBucketCounts(SGS,CIS,TRS)
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(TRStruct), intent(in) :: TRS

  integer(kind=iwp) :: ICLASS, IVLT, ISGT
  integer(kind=iwp) :: CountBuckets, CountDirect

  write(u6,*) ' '
  write(u6,*) 'MKTRANS: checking transition bucket counts by class'

  do ICLASS = 0, TRS%nClass-1

    ! Count transitions from the bucket table
    CountBuckets = 0
    do IVLT = 1, SGS%nVert
      CountBuckets = CountBuckets + TRS%NTR(IVLT,ICLASS)
    end do

    ! Count transitions directly from ISGM / ITVPT
    CountDirect = 0
    do IVLT = 1, SGS%nVert
      do ISGT = 1, nSeg
        if (CIS%ISGM(IVLT,ISGT) == 0) cycle
        if (ITVPT(ISGT) /= ICLASS) cycle
        CountDirect = CountDirect + 1
      end do
    end do

    write(u6,'(A,I3,A,I8,A,I8)') '  class=',ICLASS, &
         '  bucket sum=',CountBuckets, '  direct count=',CountDirect

    if (CountBuckets /= CountDirect) then
      write(u6,*) 'MKTRANS ERROR: transition bucket mismatch.'
      write(u6,*) '  class       =', ICLASS
      write(u6,*) '  bucket sum  =', CountBuckets
      write(u6,*) '  direct count=', CountDirect
      call Abend()
    end if

  end do

  write(u6,*) 'MKTRANS: bucket-count check passed.'
  write(u6,*) ' '

end subroutine CheckTransitionBucketCounts


integer(kind=iwp) function CountAllTransitions(SGS,CIS,nClass)
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: nClass

  integer(kind=iwp) :: IVLT, ISGT

  CountAllTransitions = 0
  do IVLT = 1, SGS%nVert
    do ISGT = 1, nSeg
      if (CIS%ISGM(IVLT,ISGT) /= 0) then
        if (ITVPT(ISGT) >= 0 .and. ITVPT(ISGT) < nClass) then
          CountAllTransitions = CountAllTransitions + 1
        end if
      end if
    end do
  end do
end function CountAllTransitions

!#endif

end subroutine MKTRANS

end module SGUGA
