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
! Construction of the split-graph tables: segment values (MKSEG),
! coupling tables (MkNRCOUP, MKCOUP), walk numbering (MKSGNUM),
! transition tables (MKTRANS).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Tables

use SGUGA_Core

#include "intent.fh"

implicit none

contains

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
      IAL = SGS%DRT(IVL,CSFG_ATab)! Pick up the a- and b-value of the left vertex
      IBL = SGS%DRT(IVL,CSFG_BTab)
      INL = 2*IAL+IBL
      do IV=IVL+1,IV2         ! loop over all vertices of level LEV, right to left, but the first vertex
        IA = SGS%DRT(IV,CSFG_ATab)! Pick up the a-value of the right vertex
        IB = SGS%DRT(IV,CSFG_BTab)
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
      IAL = SGS%DRT(IVL,CSFG_ATab)
      IBL = SGS%DRT(IVL,CSFG_BTab)
      do IV=IV1,IV2
        IA = SGS%DRT(IV,CSFG_ATab)
        IB = SGS%DRT(IV,CSFG_BTab)
        if ((IA == IAL+1) .and. (IB == IBL-2)) CIS%IVR(IVL,3) = IV
        if ((IA == IAL-1) .and. (IB == IBL+2)) CIS%IVR(IVL,4) = IV
      end do
    end do
  end do

! CONSTRUCT THE MVL AND MVR TABLES:
  EXS%MVR(:,:) = 0
  do IVL=SGS%MVSta,SGS%MVEnd   ! Loop over vertices of the MIDLEV (absolute index)
    MVLL = IVL-SGS%MVSta+1     ! Compute relative index
    ! If there is a valid node for delta(b)=-1, store relative index.
    if (CIS%IVR(IVL,1) /= 0) EXS%MVR(MVLL,1) = CIS%IVR(IVL,1)-SGS%MVSta+1
    if (CIS%IVR(IVL,2) /= 0) EXS%MVR(MVLL,2) = CIS%IVR(IVL,2)-SGS%MVSta+1 ! Dito delta(b)=+1
  end do

! constructe the recipical version MVL
  EXS%MVL(:,:) = 0
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
  end do
  write(u6,*) ' MVR TABLE:'
  do MV=1,CIS%nMidV
    write(u6,1234) MV,EXS%MVR(MV,1),EXS%MVR(MV,2)
  end do
1234 format('  MV=',I2,'    UPPER WALKS:',2I6)
  write(u6,*)
  write(u6,*) ' VERTEX PAIR TABLE IVR IN MKSEG:'
  do IVL=1,SGS%nVert
    write(u6,1235) IVL,CIS%IVR(IVL,1),CIS%IVR(IVL,2)
  end do
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

      IB = SGS%DRT(IVLT,CSFG_BTab)    ! Pick up the b-value from the DRT
      if (ISGT >= iSegTHeadBeg) then
        ! triplet operator: the ASTA values, referred to the b of the right
        ! (ket) top vertex; zero segments (b = 0 with db = 0) are not valid
        IBR = SGS%DRT(IVRT,CSFG_BTab)
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
          V = sqrt(real(IB,kind=wp)/real(1+IB,kind=wp))
        case (6)
          V = (-One)**IB*sqrt(real(1+IB,kind=wp)/real(2+IB,kind=wp))
        case (7)
          V = sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
        case (8)
          V = (-One)**IB*sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
        case (9)
          V = (-One)**IB*One/sqrt(real(IB,kind=wp)*real(1+IB,kind=wp))
        case (10)
          V = sqrt(real(IB,kind=wp)*real(2+IB,kind=wp))/real(1+IB,kind=wp)
        case (11)
          V = (-One)**IB*One/sqrt(real(1+IB,kind=wp)*real(2+IB,kind=wp))
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
  type(TRStruct), intent(in) :: TRS

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

  NRL_OpenBlock = nOpenBands*SGS%nLev
  EXS_OpenBlock = 0
  do band=1,nOpenBands
    if (ActiveBand(band)) EXS_OpenBlock = EXS_OpenBlock+SGS%nLev
  end do

  EXS%MxEO = NRL_OpenBlock+(SGS%nLev*(SGS%nLev+1))/2
  Memory = EXS_OpenBlock+(SGS%nLev*(SGS%nLev+1))/2

  call mma_allocate(EXS%NOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
  call mma_allocate(EXS%IOCP,Memory,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
  EXS%NOCP(:,:,:) = 0
  EXS%IOCP(:,:,:) = 0

! NRL(sym,vertex,indeo): number of valid segment paths of a given symmetry and operators class arriving at a given vertex
  call mma_allocate(NRL,[1,SGS%nSym],[1,SGS%nVert],[0,EXS%MxEO],Label='NRL')
! indeo=0 denotes a ordinary half-walk with no open- or closed-loop attached
  INDEO0 = 0

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
  do IVLT=1,SGS%MVSta-1    ! This formally denotes the upper left node of an segment
    LEV = SGS%DRT(IVLT,CSFG_LTab)

    ! loop over all possible transition classes: upper walk, delta(b)=-1/+1, or lower walk
    ! A segment is classified by its top and bottom transition class
    do ICLASS=0,TRS%nClass-1

      IT0 = TRS%ITR0(IVLT,ICLASS)   ! Offset to list of valid segments
      NT = TRS%NTR(IVLT,ICLASS)    ! # of valid segments for this vertex

      if (NT == 0) cycle

      do K=1,NT                  ! loop over all segment that got us here!
        ITR = IT0+K              ! counter
        ISGT = TRS%ISGT(ITR)        ! segment index
        IVLB = TRS%IVLB(ITR)        ! vertex index of lower left node of the segments
        ITOP = TRS%ITOP(ITR)        ! index of top connection type (delta(b)=+1 or -1), ITOP=ICLASS
        IBOT = TRS%IBOT(ITR)        ! dito lower connection type
        ISYM = TRS%ISYM(ITR)        ! the symmetry index of the segment

        do ITSYM=1,SGS%nSym      ! loop over the symmetry index of the top left  vertex
          IBSYM = Mul(ITSYM,ISYM)   ! get the symmetry index of left bottom vertex given the symmetry index of the segment

          select case (TRS%IFLAG(ITR))

            case (TR_WALK)
              ! ordinary upper walk
              NRL(IBSYM,IVLB,INDEO0) = NRL(IBSYM,IVLB,INDEO0)+NRL(ITSYM,IVLT,INDEO0)

            case (TR_OPEN)
              ! loop opening
              IP = LEV
              INDEO = IP+(OpenBand(IBOT)-1)*SGS%nLev
              NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,INDEO0)

            case (TR_MID)
              ! intermediate open-loop propagation

              ! Include propagations from all active open loops IP>LEV
              do IP=LEV+1,SGS%nLev
                INDEOT = IP+(OpenBand(ITOP)-1)*SGS%nLev
                INDEOB = IP+(OpenBand(IBOT)-1)*SGS%nLev
                NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
              end do

            case (TR_CLOSE)
              ! loop closing
              IQ = LEV
              ! Include propagations from all active loops IP>LEV
              do IP=LEV+1,SGS%nLev
                INDEOT = IP+(OpenBand(ITOP)-1)*SGS%nLev
                IPQ = (IP*(IP-1))/2+IQ
                INDEOB = NRL_OpenBlock+IPQ
                NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
              end do

            case (TR_TAIL)
              ! tail/downwalk propagation of closed loops
              do IPQ=1,SGS%nLEV*(SGS%nLev+1)/2
                INDEO = NRL_OpenBlock+IPQ
                NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,INDEO)
              end do

            case default
              write(u6,*) 'MkNRCOUP(upper): unexpected TRS%IFLAG = ',TRS%IFLAG(ITR)
              write(u6,*) '  ITR  = ',ITR
              write(u6,*) '  ISGT = ',ISGT
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
      do band=1,nOpenBands
        if (.not. ActiveBand(band)) cycle

        do IP=1,SGS%nLev
          INDEO = IP+(band-1)*SGS%nLev
          EXS%NOCP(INDEO,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO)
        end do

      end do

      ! ---- closed loops ----
      do IP=1,SGS%nLev
        do IQ=1,IP
          IPQ = IP*(IP-1)/2+IQ
          INDEO_EXS = IPQ+EXS_OpenBlock
          INDEO_NRL = IPQ+NRL_OpenBlock
          EXS%NOCP(INDEO_EXS,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO_NRL)
        end do
      end do

    end do
  end do

! For lower walks
  NRL(:,SGS%MVSta:SGS%nVert,:) = 0
  NRL(1,SGS%nVert,0) = 1

  do IVLT=SGS%nVert-1,SGS%MVSta,-1
    LEV = SGS%DRT(IVLT,CSFG_LTab)

    do ICLASS=0,TRS%nClass-1

      IT0 = TRS%ITR0(IVLT,ICLASS)
      NT = TRS%NTR(IVLT,ICLASS)

      if (NT == 0) cycle

      do K=1,NT
        ITR = IT0+K
        ISGT = TRS%ISGT(ITR)
        IVLB = TRS%IVLB(ITR)
        ITOP = TRS%ITOP(ITR)
        IBOT = TRS%IBOT(ITR)
        ISYM = TRS%ISYM(ITR)

        do ITSYM=1,SGS%nSym
          IBSYM = Mul(ITSYM,ISYM)

          select case (TRS%IFLAG(ITR))

            case (TR_TAIL)
              ! ordinary lower walk
              NRL(ITSYM,IVLT,INDEO0) = NRL(ITSYM,IVLT,INDEO0)+NRL(IBSYM,IVLB,INDEO0)

            case (TR_CLOSE)
              ! create open-loop class from below
              IQ = LEV
              INDEO = IQ+(OpenBand(ITOP)-1)*SGS%nLev
              NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,INDEO0)

            case (TR_MID)
              ! propagate open-loop class upward
              do IQ=1,LEV-1
                INDEOT = IQ+(OpenBand(ITOP)-1)*SGS%nLev
                INDEOB = IQ+(OpenBand(IBOT)-1)*SGS%nLev
                NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
              end do

            case (TR_OPEN)
              ! close open loop into closed-loop class
              IP = LEV
              do IQ=1,LEV-1
                INDEOB = IQ+(OpenBand(IBOT)-1)*SGS%nLev
                IPQ = (IP*(IP-1))/2+IQ
                INDEOT = NRL_OpenBlock+IPQ
                NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
              end do

            case (TR_WALK)
              ! propagate closed loops upward
              do IPQ=1,(LEV*(LEV-1))/2
                INDEO = NRL_OpenBlock+IPQ
                NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,INDEO)
              end do

            case default
              write(u6,*) 'MkNRCOUP(lower): unexpected TRS%IFLAG = ',TRS%IFLAG(ITR)
              write(u6,*) '  ITR  = ',ITR
              write(u6,*) '  ISGT = ',ISGT
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
      do band=1,nOpenBands
        if (.not. ActiveBand(band)) cycle

        do IQ=1,SGS%nLev
          INDEO = IQ+(band-1)*SGS%nLev
          EXS%NOCP(INDEO,LFTSYM,MV) = max(EXS%NOCP(INDEO,LFTSYM,MV),NRL(LFTSYM,IVLT,INDEO))
        end do

      end do

      ! ---- closed loops ----
      do IP=1,SGS%nLev
        do IQ=1,IP
          IPQ = IP*(IP-1)/2+IQ
          INDEO_EXS = IPQ+EXS_OpenBlock
          INDEO_NRL = IPQ+NRL_OpenBlock
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
        EXS%nICOup = EXS%nICOup+EXS%NOCP(INDEO,LFTSYM,MV)
      end do
    end do
  end do

  call CSFG_Mk_IOW(CIS,SGS)

  call CSFG_Count(CIS,SGS)

! INSERT FOR USE IN SIGMA ROUTINE

  NSGMX = 1
  NDWNS1 = 0  ! Dummy initialization
  NSGTMP = MXUP*MXDWN   ! Max size of intermediate sigma block

  do MVSGM=1,CIS%nMidV
    MV1 = EXS%MVL(MVSGM,2)
    MV2 = EXS%MVL(MVSGM,1)

    MV4 = EXS%MVR(MVSGM,1)
    MV5 = EXS%MVR(MVSGM,2)

    do ISYUS1=1,SGS%nSym
      NUPS1 = CIS%NOW(1,ISYUS1,MVSGM)
      do ISYDS1=1,SGS%nSym
        NDWNS1 = CIS%NOW(2,ISYDS1,MVSGM)

        NSGMX = max(NSGMX,CIS%NOCSF(ISYUS1,MVSGM,ISYDS1))

        if (MV1 /= 0) then
          NT4TMP = NUPS1*CIS%NOW(2,ISYDS1,MV1)
          NT5TMP = CIS%NOW(1,ISYUS1,MV1)*NDWNS1
          NSGTMP = max(NSGTMP,NT4TMP,NT5TMP)
        end if

        if (MV2 /= 0) then
          NT3TMP = NUPS1*CIS%NOW(2,ISYDS1,MV2)
          NT5TMP = CIS%NOW(1,ISYUS1,MV2)*NDWNS1
          NSGTMP = max(NSGTMP,NT3TMP,NT5TMP)
        end if

        if (MV4 /= 0) then
          NT1TMP = NUPS1*CIS%NOW(2,ISYDS1,MV4)
          NT5TMP = CIS%NOW(1,ISYUS1,MV4)*NDWNS1
          NSGTMP = max(NSGTMP,NT1TMP,NT5TMP)
        end if

        if (MV5 /= 0) then
          NT2TMP = NUPS1*CIS%NOW(2,ISYDS1,MV5)
          NT5TMP = CIS%NOW(1,ISYUS1,MV5)*NDWNS1
          NSGTMP = max(NSGTMP,NT2TMP,NT5TMP)
        end if

      end do
    end do
  end do

! Set up the code to wheater or not intermediate sigma vectors are to be reused.
  call mma_maxDBLE(NT1TMP)
  if (NSGTMP*2*SGS%nSym*SGS%MidLev < NT1TMP/4) then
    call mma_allocate(EXS%SGTMP,NSGTMP*2*SGS%nSym*SGS%MidLev,Label='EXS%SGTMP')
    EXS%Reuse_SGTMP = .true.
!  EXS%Reuse_SGTMP=.false.
  else
    call mma_allocate(EXS%SGTMP,NSGTMP,Label='EXS%SGTMP')
    EXS%Reuse_SGTMP = .false.
  end if
  EXS%SGTMP(:) = Zero
  EXS%Save_Valid = .false.

#ifdef _DEBUGPRINT_
  write(u6,600) MXUP,MXDWN,NSGMX,NSGMX,NSGTMP

  NUW = CIS%NUW
  write(u6,*)
  write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
  write(u6,*) '                    LOWER ',CIS%nWalk-NUW
  write(u6,*) '                     SUM  ',CIS%nWalk
  write(u6,*) ' TOTAL NR OF COUPL COEFFS ',EXS%nICOup
  INDEO = EXS_OpenBlock+1
  write(u6,*) '         OF TYPE 1&2 ONLY:',EXS%IOCP(INDEO,1,1)
  write(u6,*)
  write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
  write(u6,'(8(1X,I8))')(CIS%NCSF(IS),IS=1,SGS%nSym)
  write(u6,*)

  write(u6,*)
  write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN MkNRCOUP'
  write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
  do MV=1,CIS%nMidV
    write(u6,*)
    write(u6,1234) MV,(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
    write(u6,1235)(CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
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

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(inout) :: CIS
  type(EXStruct), intent(inout) :: EXS
  type(TRStruct), intent(in) :: TRS

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
  real(kind=wp), allocatable :: val(:), VTab(:)

  ! Rows of ISGPTH
  integer(kind=iwp), parameter :: IVLFT = 1
  integer(kind=iwp), parameter :: ITYPE = 2
  integer(kind=iwp), parameter :: IAWSL = 3
  integer(kind=iwp), parameter :: IAWSR = 4
  integer(kind=iwp), parameter :: ILS = 5
  integer(kind=iwp), parameter :: ICS = 6
  integer(kind=iwp), parameter :: ISEG = 7
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

  MxEO_Block = nOpenBands*SGS%nLev
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
  hsize = 2*nVTab+1
  call mma_allocate(VHashKey,[0,hsize-1],Label='VHashKey')
  call mma_allocate(VHashVal,[0,hsize-1],Label='VHashVal')
  VHashKey(:) = -huge(VHashKey)
  VHashVal(:) = 0

  ! Coupling coefficient value table
  NVTAB_FINAL = 2
  VTab(1) = One
  VTab(2) = -One

  NCHECK = 0

  do IHALF=1,2  ! Upper and lower SG factors

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

    ! For HALF==1 the level is 1, for HALF==2 the level is SGS%MidLEV
    do IVTOP=IVTSTA,IVTEND ! Loop over all vertex on this level, left to right

      do ITYP=0,ITYPMX

        IVRTOP = IVTOP   !
        if (ITYP > 0) IVRTOP = CIS%IVR(IVTOP,ITYP)
        if (IVRTOP == 0) cycle

        LEV = LEV1

        ! Initialize path table at starting level
        ISGPTH(IVLFT,LEV) = IVTOP
        ISGPTH(ITYPE,LEV) = ITYP
        ISGPTH(IAWSL,LEV) = 0
        ISGPTH(IAWSR,LEV) = 0
        ISGPTH(ILS,LEV) = 1
        ISGPTH(ICS,LEV) = 0
        ISGPTH(ISEG,LEV) = 0   ! bucket cursor
        ISGPTH(IRSEG,LEV) = 0   ! raw segment number

        val(LEV) = One

        ! Walk from top to midlevel, or from midlevel to root
        do while (LEV <= LEV1)

          ITYPT = ISGPTH(ITYPE,LEV)
          IVLT = ISGPTH(IVLFT,LEV)

          IT0 = TRS%ITR0(IVLT,ITYPT)
          NT = TRS%NTR(IVLT,ITYPT)

          K = ISGPTH(ISEG,LEV)+1

          if (K > NT) then
            ! No more transitions left in this bucket: reset and backtrack
            ISGPTH(ISEG,LEV) = 0
            ISGPTH(IRSEG,LEV) = 0
            LEV = LEV+1
            cycle
          end if

          ITR = IT0+K
          ISGT = TRS%ISGT(ITR)
          IVLB = TRS%IVLB(ITR)

          ICL = TRS%ICL(ITR)
          ICR = TRS%ICR(ITR)
          ISYM = TRS%ISYM(ITR)

          ! Right upper vertex
          IVRT = IVLT
          if (TRS%IPRT(ITR) /= 0) IVRT = CIS%IVR(IVLT,TRS%IPRT(ITR))

          ! Store current bucket cursor and raw segment number
          ISGPTH(ISEG,LEV) = K
          ISGPTH(IRSEG,LEV) = ISGT
          ISGPTH(ICS,LEV) = ICL

          ! Descend one level
          LEV = LEV-1

          ! Keep old MAW lookups in this first transition-based version
          ISGPTH(IAWSL,LEV) = ISGPTH(IAWSL,LEV+1)+SGS%MAW(IVLT,ICL)
          ISGPTH(IAWSR,LEV) = ISGPTH(IAWSR,LEV+1)+SGS%MAW(IVRT,ICR)

          val(LEV) = val(LEV+1)*TRS%VSEG(ITR)

          ISGPTH(ILS,LEV) = Mul(ISYM,ISGPTH(ILS,LEV+1))
          ISGPTH(IVLFT,LEV) = IVLB
          ISGPTH(ITYPE,LEV) = TRS%IBOT(ITR)
          ISGPTH(ISEG,LEV) = 0
          ISGPTH(IRSEG,LEV) = 0
          ISGPTH(ICS,LEV) = 0

          if (LEV > LEV2) cycle

          ! ------------------------------------------------------
          ! Bottom of current half-path reached
          ! ------------------------------------------------------
          MV = ISGPTH(IVLFT,SGS%MidLev)+1-SGS%MVSta
          LFTSYM = ISGPTH(ILS,LEV2)

          IT = ISGPTH(ITYPE,SGS%MidLev)
          if (IT == 0) IT = 3
          if (ISGPTH(ITYPE,LEV2) == 0) IT = 0

          if (IT == 0) then

            ! Ordinary walk
            ILND = 1+CIS%NOW(IHALF,LFTSYM,MV)
            IAWS = ISGPTH(IAWSL,LEV2)
            ILNDW(IAWS) = ILND
            CIS%NOW(IHALF,LFTSYM,MV) = ILND

            IPOS = CIS%IOW(IHALF,LFTSYM,MV)+(ILND-1)*CIS%nIpWlk

            do LL=LEV2+1,LEV1,CSFG_nPack
              IC = 0
              do L=min(LL+CSFG_nPack-1,LEV1),LL,-1
                IC = 4*IC+ISGPTH(ICS,L)
              end do
              IPOS = IPOS+1
              CIS%ICase(IPOS) = IC
            end do

          else

            ! Open or closed loop
            IP = 0
            IQ = 0

            do L=LEV2+1,LEV1
              ISG = ISGPTH(IRSEG,L)

              if ((ISG >= 5) .and. (ISG <= 8)) IP = L
              if ((ISG >= 19) .and. (ISG <= 22)) IQ = L
            end do

            if (IP == 0) IP = IQ

            if (IT == 3) then
              INDEO = MxEO_Block+(IP*(IP-1))/2+IQ
            else
              INDEO = SGS%nLev*(IT-1)+IP
            end if

            ! Increment the local index and store
            ICOP = 1+EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
            ! Add the offset and get the global index
            ICOP = EXS%IOCP(INDEO,LFTSYM,MV)+ICOP

            NCHECK = NCHECK+1

            if (ICOP > EXS%nICoup) then
              write(u6,*) 'ERROR in MKCOUP: ICOP > EXS%nICoup'
              write(u6,*) ' ICOP      = ',ICOP
              write(u6,*) ' nICoup    = ',EXS%nICoup
              write(u6,*) ' INDEO     = ',INDEO
              write(u6,*) ' MV        = ',MV
              write(u6,*) ' LFTSYM    = ',LFTSYM
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
            do i=1,NVTAB_FINAL
              IVTAB = i
              if (abs(C-VTab(i)) < 1.0e-10_wp) exit
            end do

            ! If not in the table add it at the end
            if (i > NVTAB_FINAL) then
              NVTAB_FINAL = NVTAB_FINAL+1
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
          LEV = LEV+1

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
  do ICOP=1,EXS%nICoup
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

  write(u6,*) 'NR OF DIFFERENT COUPLING VALUES: ',NVTAB_FINAL
  ! The value  1.0 is stored at i1=1
  ! The value -1.0 is stored at i1=2
  do ICOP=1,EXS%nICoup
    i1 = EXS%ICoup(3,ICOP)
    if (i1 == 1) ICOP1 = ICOP1+1
    if (i1 == 2) ICOP2 = ICOP2+1
  end do

  write(u6,*) 'NR OF COUPS WITH VALUE  1.0: ',ICOP1
  write(u6,*) 'NR OF COUPS WITH VALUE -1.0: ',ICOP2
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
          if (JPOS == CSFG_nPack+1) then
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
          if (JPOS == CSFG_nPack+1) then
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

  call mma_deallocate(TRS%NTR,safe='*')
  call mma_deallocate(TRS%ITR0,safe='*')

  call mma_deallocate(TRS%ISGT,safe='*')
  call mma_deallocate(TRS%IVLT,safe='*')
  call mma_deallocate(TRS%IVLB,safe='*')

  call mma_deallocate(TRS%ITOP,safe='*')
  call mma_deallocate(TRS%IBOT,safe='*')

  call mma_deallocate(TRS%ICL,safe='*')
  call mma_deallocate(TRS%ICR,safe='*')

  call mma_deallocate(TRS%IPRT,safe='*')
  call mma_deallocate(TRS%ISYM,safe='*')
  call mma_deallocate(TRS%IOBAND,safe='*')
  call mma_deallocate(TRS%IFLAG,safe='*')

  call mma_deallocate(TRS%MAWL,safe='*')
  call mma_deallocate(TRS%MAWR,safe='*')

  call mma_deallocate(TRS%VSEG,safe='*')

end subroutine TRANS_Free

subroutine MKTRANS(SGS,CIS,TRS)

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(TRStruct), intent(inout) :: TRS
  integer(kind=iwp) :: IBOT, ICL, ICLASS, ICR, IFLAG, IOBAND, IPRT, ISGT, ISYM, ITOP, ITR, IVLB, IVLT, LEV, N
  real(kind=wp) :: VSEG
  integer(kind=iwp), allocatable :: IPOS(:,:)

  ! Initialize metadata
  call TRANS_Free(TRS)

  TRS%nClass = 4
  TRS%nOpenBand = 2

  ! Allocate bucket tables
  call mma_deallocate(TRS%NTR,safe='*')
  call mma_deallocate(TRS%ITR0,safe='*')
  call mma_allocate(TRS%NTR,[1,SGS%nVert],[0,TRS%nClass-1],Label='TRS%NTR')
  call mma_allocate(TRS%ITR0,[1,SGS%nVert],[0,TRS%nClass-1],Label='TRS%ITR0')

  TRS%NTR(:,:) = 0
  TRS%ITR0(:,:) = 0

  ! First pass: couny transitions per bucket
  ! For every source vertex IVLT, count how many valid segments require each top class ITOP
  do IVLT=1,SGS%nVert
    do ISGT=1,nSeg
      IVLB = CIS%ISGM(IVLT,ISGT)
      if (IVLB == 0) cycle
      ITOP = ITVPT(ISGT)
      TRS%NTR(IVLT,ITOP) = TRS%NTR(IVLT,ITOP)+1
    end do
  end do

  ! Build offsets
  N = 0
  do ICLASS=0,TRS%nClass-1
    do IVLT=1,SGS%nVert
      TRS%ITR0(IVLT,ICLASS) = N
      N = N+TRS%NTR(IVLT,ICLASS)
    end do
  end do

  TRS%nTrans = N
  ! Now the bucket layout is fixed.

  ! Allocate the flat transition arrays.

  call mma_allocate(TRS%ISGT,TRS%nTrans,Label='TRS%ISGT')
  call mma_allocate(TRS%IVLT,TRS%nTrans,Label='TRS%IVLT')
  call mma_allocate(TRS%IVLB,TRS%nTrans,Label='TRS%IVLB')

  call mma_allocate(TRS%ITOP,TRS%nTrans,Label='TRS%ITOP')
  call mma_allocate(TRS%IBOT,TRS%nTrans,Label='TRS%IBOT')

  call mma_allocate(TRS%ICL,TRS%nTrans,Label='TRS%ICL')
  call mma_allocate(TRS%ICR,TRS%nTrans,Label='TRS%ICR')

  call mma_allocate(TRS%IPRT,TRS%nTrans,Label='TRS%IPRT')
  call mma_allocate(TRS%ISYM,TRS%nTrans,Label='TRS%ISYM')
  call mma_allocate(TRS%IOBAND,TRS%nTrans,Label='TRS%IOBAND')
  call mma_allocate(TRS%IFLAG,TRS%nTrans,Label='TRS%IFLAG')

  call mma_allocate(TRS%MAWL,TRS%nTrans,Label='TRS%MAWL')
  call mma_allocate(TRS%MAWR,TRS%nTrans,Label='TRS%MAWR')

  call mma_allocate(TRS%VSEG,TRS%nTrans,Label='TRS%VSEG')

  ! Allocate a fill-position array
  call mma_allocate(IPOS,[1,SGS%nVert],[0,TRS%nClass-1],Label='IPOS')
  IPOS(:,:) = TRS%ITR0(:,:)   ! IPOS is a write pointer per bucket.

  ! Second pass: fill transition arrays

  do IVLT=1,SGS%nVert      ! Loop over all vertices
    LEV = SGS%DRT(IVLT,CSFG_LTab)  ! Pick up the level index of the vertex

    ! loop over all possible segment types allowed for this vertex.
    do ISGT=1,nSeg

      IVLB = CIS%ISGM(IVLT,ISGT) ! Get the index of Vertex Left Bottom
      if (IVLB == 0) cycle

      ITOP = ITVPT(ISGT)   ! top connection class
      IBOT = IBVPT(ISGT)   ! bottom connection class

      ICL = IC1(ISGT)      ! left case
      ICR = IC2(ISGT)      ! right case

      ! symmetry index of the vertex pair (iVLT-iVLB)
      ISYM = 1
      if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)

      IPRT = PartnerSlot(ITOP)
      IOBAND = OpenBand(ITOP)
      VSEG = CIS%VSGM(IVLT,ISGT)

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

      IPOS(IVLT,ITOP) = IPOS(IVLT,ITOP)+1
      ITR = IPOS(IVLT,ITOP)

      TRS%ISGT(ITR) = ISGT
      TRS%IVLT(ITR) = IVLT
      TRS%IVLB(ITR) = IVLB

      TRS%ITOP(ITR) = ITOP
      TRS%IBOT(ITR) = IBOT

      TRS%ICL(ITR) = ICL
      TRS%ICR(ITR) = ICR

      TRS%IPRT(ITR) = IPRT
      TRS%ISYM(ITR) = ISYM
      TRS%IOBAND(ITR) = IOBAND
      TRS%IFLAG(ITR) = IFLAG

      TRS%VSEG(ITR) = VSEG

      ! Leave these as placeholders in the first version
      TRS%MAWL(ITR) = 0
      TRS%MAWR(ITR) = 0

    end do
  end do

  call mma_deallocate(IPOS)

end subroutine MKTRANS

!> The walk tables of state iState for the CSFs of symmetry LSM: ICase
!> (CSFG_MkCOT) and the walk numbering USGN/LSGN (MkSgNum), which depends on the
!> symmetry and is rebuilt when another symmetry is asked for (a state can
!> serve several symmetries, e.g. MCLR's state for the whole run).
subroutine SG_Num_Setup(iState,LSM)
  integer(kind=iwp), intent(in) :: iState, LSM
  if (.not. allocated(CIS(iState)%ICase)) call CSFG_MkCOT(SGS(iState),CIS(iState))
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

end module SGUGA_Tables
