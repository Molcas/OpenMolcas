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

use Molcas, only: MxLev
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
use Constants, only: Zero, One

implicit none
private

! Split-Graph descriptor, sizes, addresses...
type SGStruct
  integer(kind=iwp) :: NSym = 0, nActEl = 0, IFRAS = 0
  integer(kind=iwp) :: IA0, IB0, IC0, iSpin, nLev, nVert, nVert0, MidLev, MVSta, MVEnd, MXUP, MXDWN, LV1RAS, LM1RAS, LV3RAS, LM3RAS
  integer(kind=iwp), allocatable :: ISm(:), DRT(:,:), DRT0(:,:), Down(:,:), Down0(:,:), Up(:,:), Ver(:), MAW(:,:), LTV(:), &
                                    DAW(:,:), RAW(:,:), SCR(:,:)
  integer(kind=iwp), pointer :: DRTP(:,:), DOWNP(:,:)
end type SGStruct

! CI Structures, addresses,..
type CIStruct
  integer(kind=iwp) :: nMidV, nIpWlk, nWalk
  integer(kind=iwp), allocatable :: NOW(:,:,:), IOW(:,:,:), NCSF(:), NOCSF(:,:,:), IOCSF(:,:,:), ICase(:), IVR(:,:), ISGM(:,:)
  real(kind=wp), allocatable :: VSGM(:,:)
end type CIStruct

! Excitation operators, coupling coefficients,...
type EXStruct
  integer(kind=iwp) :: MxEO, nICoup
  integer(kind=iwp), allocatable :: NOCP(:,:,:), IOCP(:,:,:), ICoup(:,:), MVL(:,:), MVR(:,:), USGN(:,:), LSGN(:,:)
  real(kind=wp), allocatable :: VTab(:), SGTMP(:)
end type EXStruct

type(SGStruct), target :: SGS
type(CIStruct), target :: CIS
type(EXStruct), target :: EXS


integer(kind=iwp) :: iq
integer(kind=iwp), protected :: L2ACT(MXLEV)=[(iq,iq=1,MXLEV)]
integer(kind=iwp), protected :: LEVEL(MXLEV)=[(iq,iq=1,MXLEV)]

public :: CIS, CIStruct, EXS, EXStruct, L2ACT, LEVEL, SGS, SGStruct

public :: SG_Init, MKSGUGA, MkCOT, MkCList, MkMAW, MkSeg, NrCOUP, MkCoup, MkISM_rasscf, MkSgNum, SG_Free
public :: SG_Init_Simple, MKISM_Raw

integer(kind=iwp), parameter :: IBVPT(26) = [0,0,0,0,1,1,2,2,1,1,2,1,1,2,2,1,2,2,3,3,3,3,3,3,3,3], &
                                IC1(26)   = [0,1,2,3,0,2,0,1,0,1,1,2,3,0,1,2,2,3,1,3,2,3,0,1,2,3], &
                                IC2(26)   = [0,1,2,3,1,3,2,3,0,1,2,2,3,0,1,1,2,3,0,2,0,1,0,1,2,3], &
                                ISVC(26)  = [1,1,1,1,1,6,1,5,1,2,4,7,2,1,7,3,2,2,1,5,1,6,1,1,1,1], &
                                ITVPT(26) = [0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,1,1,2,2,3,3,3,3]

contains

subroutine MKSGUGA(SGS,CIS)
! PURPOSE: MAKE THE GUGA TABLES
! NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!          THE START ADDRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!          THREE USER DEFINED TYPES. Consult the gugx module for the details.

  type(SGStruct), target, intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2, ATAB = 3, BTAB = 4, CTAB = 5
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
    call mkRAS(SGS)

    ! REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)

    call mkDRT(SGS)

    ! IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED DRT TABLE

  end if

  nullify(SGS%DOWNP,SGS%DRTP)

  ! CALCULATE ARC WEIGHT.

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
    IAC = min(SGS%IA0,SGS%IC0)
    SGS%nVert0 = ((SGS%IA0+1)*(SGS%IC0+1)*(2*SGS%IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

  end subroutine mknVert0

  subroutine mkDRT0(SGS)
  ! PURPOSE: CONSTRUCT THE UNRESTRICTED GUGA TABLE

    use Index_Functions, only: nTri_Elem1

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: ADDR, ADWN, AUP, BC, BDWN, BUP, CDWN, CUP, DWN, LEV, MXADDR, NACTEL, nTmp, STEP, VDWN, VEND, VSTA, VUP, &
                         VUPS
    integer(kind=iwp), allocatable :: TMP(:)
    integer(kind=iwp), parameter :: DA(0:3) = [0,0,1,1], DB(0:3) = [0,1,-1,0], DC(0:3) = [1,0,1,0]
#   ifdef _DEBUGPRINT_
    integer(kind=iwp) :: VERT
#   endif

    NTMP = nTri_Elem1(SGS%nLev+1)
    call mma_allocate(TMP,NTMP,Label='TMP')

    ! SET UP TOP ROW

    NACTEL = 2*SGS%IA0+SGS%IB0
    SGS%nLev = SGS%IA0+SGS%IB0+SGS%IC0
    SGS%DRTP(1,LTAB) = SGS%nLev
    SGS%DRTP(1,NTAB) = NACTEL
    SGS%DRTP(1,ATAB) = SGS%IA0
    SGS%DRTP(1,BTAB) = SGS%IB0
    SGS%DRTP(1,CTAB) = SGS%IC0
    VSTA = 1
    VEND = 1
#   ifdef _DEBUGPRINT_
    write(u6,*) 'A0,B0,C0,NVERT=',SGS%IA0,SGS%IB0,SGS%IC0,SGS%nVert
#   endif

    ! LOOP OVER ALL LEVELS

    do LEV=SGS%nLev,1,-1
      MXADDR = ((LEV+1)*(LEV+2))/2
      TMP(1:MXADDR) = 0

      ! LOOP OVER VERTICES

      do VUP=VSTA,VEND
        AUP = SGS%DRTP(VUP,ATAB)
        BUP = SGS%DRTP(VUP,BTAB)
        CUP = SGS%DRTP(VUP,CTAB)

        ! LOOP OVER CASES
        ! AND STORE ONLY VALID CASE NUMBERS WITH ADDRESSES

        do STEP=0,3
          SGS%DownP(VUP,STEP) = 0
          ADWN = AUP-DA(STEP)
          if (ADWN < 0) cycle
          BDWN = BUP-DB(STEP)
          if (BDWN < 0) cycle
          CDWN = CUP-DC(STEP)
          if (CDWN < 0) cycle
          BC = BDWN+CDWN
          ADDR = 1+(BC*(BC+1))/2+CDWN
          TMP(ADDR) = 4*VUP+STEP
          SGS%DownP(VUP,STEP) = ADDR
        end do
      end do
      VDWN = VEND

      ! NOW INSERT VALID CASES INTO DRT TABLE

      do ADDR=1,MXADDR
        VUPS = TMP(ADDR)
        if (VUPS == 0) cycle
        VUP = VUPS/4
        STEP = mod(VUPS,4)
        VDWN = VDWN+1
        SGS%DRTP(VDWN,ATAB) = SGS%DRTP(VUP,ATAB)-DA(STEP)
        SGS%DRTP(VDWN,BTAB) = SGS%DRTP(VUP,BTAB)-DB(STEP)
        SGS%DRTP(VDWN,CTAB) = SGS%DRTP(VUP,CTAB)-DC(STEP)
        TMP(ADDR) = VDWN
      end do

      ! CREATE DOWN CHAIN TABLE

      do VUP=VSTA,VEND
        do STEP=0,3
          DWN = SGS%DownP(VUP,STEP)
          if (DWN /= 0) SGS%DownP(VUP,STEP) = TMP(DWN)
        end do
      end do
      VSTA = VEND+1
      VEND = VDWN
    end do
    ! End of loop over levels.

    ! ADDING THE ZERO LEVEL TO DRT AND DOWNCHAIN TABLE

    SGS%DRTP(VEND,1:5) = 0
    SGS%DownP(VEND,0:3) = 0

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

    call mma_deallocate(TMP)

  end subroutine mkDRT0

  subroutine mkRAS(SGS)

    use UnixInfo, only: ProgName
    type(SGStruct), target, intent(inout) :: SGS

    if (ProgName(1:5) == 'rassi') then
      call rmvert(SGS)
    else
      call RESTR(SGS)
    end if

  end subroutine mkRAS

  subroutine mkDRT(SGS)
  ! PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
  !          THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
  !          VERTICES FROM THE OLD TO THE NEW DRT TABLE

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, ID, IDNEW, IV, IVNEW

    call mma_allocate(SGS%DRT,SGS%nVert,5,Label='DRT')
    call mma_allocate(SGS%Down,[1,SGS%nVert],[0,3],Label='SGS%DOWN')

    do IV=1,SGS%nVert0
      IVNEW = SGS%Ver(IV)
      if (IVNEW == 0) cycle
      SGS%DRT(IVNEW,:) = SGS%DRT0(IV,1:5)
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        IDNEW = 0
        if (ID /= 0) IDNEW = SGS%Ver(ID)
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
    integer(kind=iwp) :: IC, IDWN, ISUM, IV

    call mma_allocate(SGS%DAW,[1,SGS%nVert],[0,4],Label='SGS%DAW')

    ! BEGIN TO CONSTRUCT DOWN CHAIN TABLE

    SGS%DAW(SGS%nVert,0:3) = 0
    SGS%DAW(SGS%nVert,4) = 1
    do IV=SGS%nVert-1,1,-1
      ISUM = 0
      do IC=0,3
        SGS%DAW(IV,IC) = 0
        IDWN = SGS%Down(IV,IC)
        if (IDWN == 0) cycle
        SGS%DAW(IV,IC) = ISUM
        ISUM = ISUM+SGS%DAW(IDWN,4)
      end do
      SGS%DAW(IV,4) = ISUM
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
    integer(kind=iwp) :: IC, IDWN, ISUM, IU, IV

    call mma_allocate(SGS%UP,[1,SGS%nVert],[0,3],Label='SGS%UP')
    call mma_allocate(SGS%RAW,[1,SGS%nVert],[0,4],Label='SGS%RAW')

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

    SGS%RAW(:,0:3) = 0
    SGS%RAW(1,4) = 1
    do IV=2,SGS%nVert
      ISUM = 0
      do IC=0,3
        IU = SGS%Up(IV,IC)
        if (IU == 0) cycle
        SGS%RAW(IV,IC) = ISUM
        ISUM = ISUM+SGS%RAW(IU,4)
      end do
      SGS%RAW(IV,4) = ISUM
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
  ! PURPOSE: FIND THE MIDLEVEL

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IV, LEV

    call mma_allocate(SGS%LTV,[-1,SGS%nLev],Label='LTV')

    ! SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:

    SGS%LTV(:) = 0

    do IV=1,SGS%nVert
      LEV = SGS%DRT(IV,LTAB)
      SGS%LTV(LEV) = SGS%LTV(LEV)+1
    end do

    do LEV=SGS%nLev,0,-1
      SGS%LTV(LEV-1) = SGS%LTV(LEV-1)+SGS%LTV(LEV)
    end do

    do LEV=-1,SGS%nLev-1
      SGS%LTV(LEV) = 1+SGS%LTV(LEV+1)
    end do

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

    if (SGS%nLev == 0) then
      SGS%MidLev = 0
    else
      SGS%MidLev = 1
    end if
    MINW = 1000000
    do IL=1,SGS%nLev-1
      NW = 0
      do IV=SGS%LTV(IL),SGS%LTV(IL-1)-1
        NW = NW+SGS%RAW(IV,4)-SGS%DAW(IV,4)
      end do
      NW = abs(NW)
      if (NW >= MINW) cycle
      SGS%MidLev = IL
      MINW = NW
    end do
    SGS%MVSta = SGS%LTV(SGS%MidLev)
    SGS%MVEnd = SGS%LTV(SGS%MidLev-1)-1
    CIS%nMidV = SGS%MVEnd-SGS%MVSta+1

    ! NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
    ! (DISREGARDING SYMMETRY)

    SGS%MxUp = 0
    SGS%MxDwn = 0
    do MV=SGS%MVSta,SGS%MVEnd
      if (SGS%MxUp < SGS%RAW(MV,4)) SGS%MxUp = SGS%RAW(MV,4)
      if (SGS%MxDwn < SGS%DAW(MV,4)) SGS%MxDwn = SGS%DAW(MV,4)
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

use RasDef, only: nRas, nRsPrt, nRasEl

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, iRO, iSy, IV, L, Lev, N, NCHANGES, NLD, NV
logical(kind=iwp) :: Test
integer(kind=iwp), allocatable :: CONN(:), Lim(:)
integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2

! Construct a restricted graph.
call mma_allocate(Lim,SGS%nLev,Label='Lim')
Lim(:) = 0
! Fill in the occupation limit table:
Lev = 0
do iRO=1,nRsPrt
  do iSy=1,SGS%nSym
    Lev = Lev+nRas(iSy,iRO)
  end do
  if (Lev > 0) Lim(Lev) = nRasEl(iRO)
end do

call mma_allocate(SGS%Ver,SGS%nVert0,Label='SGS%Ver')
call mma_allocate(CONN,SGS%nVert,Label='CONN')

! KILL VERTICES THAT DO NOT OBEY RESTRICTIONS.
do IV=1,SGS%nVert-1
  SGS%Ver(IV) = 1
  L = SGS%DRT0(IV,LTAB)
  N = SGS%DRT0(IV,NTAB)
  if (N < Lim(L)) SGS%Ver(IV) = 0
end do
SGS%Ver(SGS%nVert) = 1

NCHANGES = 1 ! Initiate first loop
do while (NCHANGES > 0)
  ! REMOVE ARCS HAVING A DEAD UPPER OR LOWER VERTEX.
  ! COUNT THE NUMBER OF ARCS REMOVED OR VERTICES KILLED.
  NCHANGES = 0
  do IV=1,SGS%nVert-1
    if (SGS%Ver(IV) == 0) then
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        if (ID > 0) then
          SGS%Down0(IV,IC) = 0
          NCHANGES = NCHANGES+1
        end if
      end do
    else
      NLD = 0
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        if (ID > 0) then
          if (SGS%Ver(ID) == 0) then
            SGS%Down0(IV,IC) = 0
            NCHANGES = NCHANGES+1
          else
            NLD = NLD+1
          end if
        end if
      end do
      if (NLD == 0) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES+1
      end if
    end if
  end do
  ! ALSO CHECK ON CONNECTIONS FROM ABOVE:
  CONN(:) = 0
  CONN(1) = SGS%Ver(1)
  do IV=1,SGS%nVert-1
    if (SGS%Ver(IV) == 1) then
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        Test = ID > 0
        if (Test) Test = SGS%Ver(ID) == 1
        if (Test) CONN(ID) = 1
      end do
    end if
  end do
  do IV=1,SGS%nVert
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

NV = 0
do IV=1,SGS%nVert
  if (SGS%Ver(IV) == 1) then
    NV = NV+1
    SGS%Ver(IV) = NV
  end if
end do
SGS%nVert = NV

call mma_deallocate(CONN)
call mma_deallocate(Lim)

end subroutine RMVERT

subroutine RESTR(SGS)
! PURPOSE: PUT THE RAS CONSTRAINT TO THE DRT TABLE BY CREATING A MASK

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, IV, IVD, IVV, LEV, MASK, N
integer(kind=iwp), parameter :: i_and(0:3,0:3) = reshape([0,0,0,0,0,1,0,1,0,0,2,2,0,1,2,3],[4,4]), &
                                i_or(0:3,0:3) = reshape([0,1,2,3,1,1,3,3,2,3,2,3,3,3,3,3],[4,4]), &
                                LTAB = 1, NTAB = 2

call mma_allocate(SGS%Ver,SGS%nVert0,Label='V11')

! LOOP OVER ALL VERTICES AND CHECK ON RAS CONDITIONS
! CREATE MASK

do IV=1,SGS%nVert0
  LEV = SGS%DRT0(IV,LTAB)
  N = SGS%DRT0(IV,NTAB)
  SGS%Ver(IV) = 0
  if ((LEV == SGS%LV1RAS) .and. (N >= SGS%LM1RAS)) SGS%Ver(IV) = 1
  if ((LEV == SGS%LV3RAS) .and. (N >= SGS%LM3RAS)) SGS%Ver(IV) = SGS%Ver(IV)+2
end do

! NOW LOOP FORWARDS, MARKING THOSE VERTICES CONNECTED FROM ABOVE.
! SINCE VER WAS INITIALIZED TO ZERO, NO CHECKING IS NEEDED.

do IV=1,SGS%nVert0-1
  IVV = SGS%Ver(IV)
  do IC=0,3
    ID = SGS%Down0(IV,IC)
    if (ID /= 0) SGS%Ver(ID) = i_or(SGS%Ver(ID),IVV)
  end do
end do

! THEN LOOP BACKWARDS. SAME RULES, EXCEPT THAT CONNECTIVITY
! SHOULD BE PROPAGATED ONLY ABOVE THE RESTRICTION LEVELS.

do IV=SGS%nVert0-1,1,-1
  LEV = SGS%DRT0(IV,LTAB)
  MASK = 0
  if (LEV > SGS%LV1RAS) MASK = 1
  if (LEV > SGS%LV3RAS) MASK = MASK+2
  IVV = SGS%Ver(IV)
  do IC=0,3
    ID = SGS%Down0(IV,IC)
    if (ID /= 0) then
      IVD = SGS%Ver(ID)
      IVV = i_or(IVV,i_and(MASK,IVD))
    end if
  end do
  SGS%Ver(IV) = IVV
end do

! WE ARE NOW INTERESTED ONLY IN VERTICES CONNECTED BOTH TO
! ALLOWED VERTICES FOR RAS-SPACE 1 AND RAS-SPACE 3.
! THOSE ARE NUMBERED IN ASCENDING ORDER, THE REST ARE ZEROED.

SGS%nVert = 0
do IV=1,SGS%nVert0
  if (SGS%Ver(IV) == 3) then
    SGS%nVert = SGS%nVert+1
    SGS%Ver(IV) = SGS%nVert
  else
    SGS%Ver(IV) = 0
  end if
end do
if (SGS%nVert == 0) call SysAbendMsg('Restr','No configuration was found\n','Check NACTEL, RAS1, RAS2, RAS3 values')
#ifdef _DEBUGPRINT_
write(u6,*) 'RESTR:'
write(u6,*) 'LV1RAS, LV3RAS, LM1RAS, LM3RAS=',SGS%LV1RAS,SGS%LV3RAS,SGS%LM1RAS,SGS%LM3RAS
do IV=1,SGS%nVert0
  write(u6,*) 'VER(:)=',SGS%Ver(IV)
end do
#endif

end subroutine RESTR

end subroutine MKSGUGA


SUBROUTINE SG_Init_Simple(nSym,nActEl,iSpin,                   &
                          SGS,CIS,EXS,                         &
                          nHole1,nEle3,nRs1,nRs2,nRs3,         &
                          xLevel,xL2Act,xNLEV,xNSM)
IMPLICIT None
Integer(kind=iwp), intent(in):: nSym,nActEl,iSpin
Type(SGStruct), intent(inout):: SGS
Type(CIStruct), intent(inout):: CIS
Type(EXStruct),  optional, intent(inout):: EXS
Integer(kind=iwp), optional, intent(in):: nHole1,nEle3,nRs1(nSym),nRs2(nSym),nRs3(nSym)
Integer(kind=iwp), optional, intent(in):: xLevel(MxLev), xL2Act(MxLev)
Integer(kind=iwp), optional, intent(in):: xNLEV, xNSM(MxLev)

Integer(kind=iwp) :: nRas1T,nRas2T,nRas3T,IS


! Make sure that we start from a clean slate.
If (Present(EXS)) THEN
   ! Here if the extended parameter list was used.
   Call SG_Free(SGS,CIS,EXS)
Else
   ! Here if the terse parameter list was used.
   Call SG_Free(SGS,CIS)
End If

If (nSym<1 .or. nSym>8) Then
   Write (u6,*) ' SG_Init_Simple: illegal nSym value:', nSym
   Call Abend()
End If
If (iSpin<1) Then
   Write (u6,*) ' SG_Init_Simple: illegal iSpin value:', iSpin
   Call Abend()
End If
If (nActEl<0) Then
   Write (u6,*) ' SG_Init_Simple: illegal nActEl value:', nActEl
   Call Abend()
End If

SGS%nSym=nSym
SGS%iSpin=iSpin
SGS%nActEl=nActEl

If (Present(EXS)) THEN

   nRAS1T = sum(nRs1(1:nSym))
   nRAS2T = sum(nRs2(1:nSym))
   nRAS3T = sum(nRs3(1:nSym))

   SGS%LV1RAS=NRAS1T
   SGS%LV3RAS=nRAS1T+NRAS2T
   SGS%LM1RAS=2*nRas1T-NHOLE1
   SGS%LM3RAS=NACTEL-NELE3
   IF ((NRAS1T+NRAS3T)/=0) Then
      SGS%IFRAS=1
      do IS=1,NSYM
         if (nRs1(IS)+nRs2(IS)+nRs3(IS) /= 0) SGS%IFRAS = SGS%IFRAS+1
      end do

   Else
      SGS%IFRAS=0
   End If
Else
   SGS%LV1RAS=0
   SGS%LV3RAS=0
   SGS%LM1RAS=0
   SGS%LM3RAS=0
End IF

If (Present(xLevel)) Level(:)=xLevel(:)
If (Present(xL2Act)) L2Act(:)=xL2Act(:)

! CREATE THE SYMMETRY INDEX VECTOR

If (Present(xnLev)) Then
   Call MKISM(SGS,xnLev,xNSM)
Else
   Call MKISM(SGS)
End If

Call MkSGUGA(SGS,CIS)

END SUBROUTINE SG_Init_Simple


SUBROUTINE SG_Init(nSym,nActEl,iSpin,                   &
                  SGS,CIS,EXS,                          &
                  nHole1,nEle3,nRs1,nRs2,nRs3,          &
                  xLevel,xL2Act)
IMPLICIT None
Integer(kind=iwp), intent(in):: nSym,nActEl,iSpin
Type(SGStruct), intent(inout):: SGS
Type(CIStruct), intent(inout):: CIS
Integer(kind=iwp), optional, intent(in):: nHole1,nEle3,nRs1(nSym),nRs2(nSym),nRs3(nSym)
Type(EXStruct),  optional, intent(inout):: EXS
Integer(kind=iwp), optional, intent(in):: xLevel(MxLev), xL2Act(MxLev)

Call SG_Init_Simple(nSym,nActEl,iSpin,                   &
                    SGS,CIS,EXS,                         &
                    nHole1,nEle3,nRs1,nRs2,nRs3,         &
                    xLevel,xL2Act)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

CALL MKMAW(SGS)

If (Present(EXS)) THEN
!     FORM VARIOUS OFFSET TABLES:
!     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
!
   CALL MKCOT(SGS,CIS)
!
!     CONSTRUCT THE CASE LIST
!
   Call MKCLIST(SGS,CIS)

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

   CALL MKSEG(SGS,CIS,EXS)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.

   CALL NRCOUP(SGS,CIS,EXS)

   CALL MKCOUP(SGS,CIS,EXS)
End If

END SUBROUTINE SG_Init

subroutine MKCOT(SGS,CIS)
! PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
! NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
!          THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
!          WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
!          SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
!          THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp) :: IHALF, ILND, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, LEV, LEV1, LEV2, MV, NUW
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

CIS%nIpWlk = 1+(SGS%MidLev-1)/15
CIS%nIpWlk = max(CIS%nIpWlk,1+(SGS%nLev-SGS%MidLev-1)/15)
call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF')
call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF')
call mma_allocate(SGS%Scr,[1,3],[0,SGS%nLev],Label='SGS%Scr')

! CLEAR ARRAYS IOW AND NOW

CIS%NOW(:,:,:) = 0
CIS%IOW(:,:,:) = 0

! CLEAR ARRAYS IOCSF AND NOCSF

CIS%IOCSF(:,:,:) = 0
CIS%NOCSF(:,:,:) = 0

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
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

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
    LEV = LEV1
    SGS%Scr(IVERT,LEV) = IVTOP
    SGS%Scr(ISYM,LEV) = 1
    SGS%Scr(ISTEP,LEV) = -1
    do while (LEV <= LEV1)
      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
      IVT = SGS%Scr(IVERT,LEV)
      do ISTP=SGS%Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) exit
      end do
      ! NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
      if (ISTP > 3) then
        SGS%Scr(ISTEP,LEV) = -1
        LEV = LEV+1
        cycle
      end if
      ! SUCH AN ARC WAS FOUND. WALK DOWN:
      SGS%Scr(ISTEP,LEV) = ISTP
      ISML = 1
      if ((ISTP == 1) .or. (ISTP == 2)) ISML = SGS%ISm(LEV)
      LEV = LEV-1
      SGS%Scr(ISYM,LEV) = Mul(ISML,SGS%Scr(ISYM,LEV+1))
      SGS%Scr(IVERT,LEV) = IVB
      SGS%Scr(ISTEP,LEV) = -1
      if (LEV > LEV2) cycle
      ! WE HAVE REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
      ! FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
      MV = SGS%Scr(IVERT,SGS%MidLev)+1-SGS%MVSta
      IWSYM = SGS%Scr(ISYM,LEV2)
      ILND = 1+CIS%NOW(IHALF,IWSYM,MV)
      ! SAVE THE MAX WALK NUMBER FOR GIVEN SYMMETRY AND MIDVERTEX
      CIS%NOW(IHALF,IWSYM,MV) = ILND
      ! BACK UP ONE LEVEL AND TRY AGAIN:
      LEV = LEV+1
    end do
  end do
end do

call CSFCOUNT(CIS,SGS%nSym,NUW)

#ifdef _DEBUGPRINT_
Block
integer(kind=iwp) :: IS
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)
write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:',(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,'(A,8I6)') '           LOWER WALKS:',(CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do ISTP=1,SGS%nSym
    write(u6,'(A,I2,A,8I6)') ' ISTP=',ISTP,'  CONFIGURATIONS:',(CIS%NOCSF(IS,MV,ISTP),IS=1,SGS%nSym)
  end do
end do
End Block
#endif

end subroutine MKCOT

subroutine MKCLIST(SGS,CIS)
! PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
!          STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
!          IN THE ARRAY ICASE. GROUPS OF 15 CASES ARE PACKED
!          INTO ONE INTEGER WORD.

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp) :: IC, IHALF, ILND, IPOS, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, L, LEV, LEV1, LEV2, LL, MV
logical(kind=iwp) :: Found
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')
call mma_allocate(SGS%Scr,[1,3],[0,SGS%nLev],Label='SGS%Scr',safe='*')

! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY

CIS%NOW(:,:,:) = 0

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
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

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH:
    LEV = LEV1
    SGS%Scr(IVERT,LEV) = IVTOP
    SGS%Scr(ISYM,LEV) = 1
    SGS%Scr(ISTEP,LEV) = -1
    do while (LEV <= LEV1)
      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX:
      IVT = SGS%Scr(IVERT,LEV)
      Found = .false.
      do ISTP=SGS%Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) then
          Found = .true.
          exit
        end if
      end do
      if (Found) then
        ! ALT A -- SUCH AN ARC WAS FOUND. WALK DOWN:
        SGS%Scr(ISTEP,LEV) = ISTP
        ISML = 1
        if ((ISTP == 1) .or. (ISTP == 2)) ISML = SGS%ISm(LEV)
        LEV = LEV-1
        SGS%Scr(ISYM,LEV) = Mul(ISML,SGS%Scr(ISYM,LEV+1))
        SGS%Scr(IVERT,LEV) = IVB
        SGS%Scr(ISTEP,LEV) = -1
        if (LEV > LEV2) cycle
        ! WE HAVE REACHED THE LOWER LEVEL. THE WALK IS COMPLETE.
        ! MIDVERTEX NUMBER:
        MV = SGS%Scr(IVERT,SGS%MidLev)+1-SGS%MVSta
        ! SYMMETRY LABEL OF THIS WALK:
        IWSYM = SGS%Scr(ISYM,LEV2)
        ! ITS ORDERING NUMBER WITHIN THE SAME BATCH OF (IHALF,IWSYM,MV):
        ILND = 1+CIS%NOW(IHALF,IWSYM,MV)
        CIS%NOW(IHALF,IWSYM,MV) = ILND
        ! CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
        IPOS = CIS%IOW(IHALF,IWSYM,MV)+(ILND-1)*CIS%nIpWlk
        ! PACK THE STEPS IN GROUPS OF 15 LEVELS PER INTEGER:
        do LL=LEV2+1,LEV1,15
          IC = 0
          do L=min(LL+14,LEV1),LL,-1
            IC = 4*IC+SGS%Scr(ISTEP,L)
          end do
          IPOS = IPOS+1
          CIS%ICase(IPOS) = IC
        end do
        ! FINISHED WITH THIS WALK. BACK UP ONE LEVEL AND TRY AGAIN:
        LEV = LEV+1
      else
        ! ALT B -- NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
        SGS%Scr(ISTEP,LEV) = -1
        LEV = LEV+1
      end if
    end do
  end do
end do

call mma_deallocate(SGS%Scr)

end subroutine MKCLIST

subroutine MKMAW(SGS)

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, ISUM, IU, IV

call mma_allocate(SGS%MAW,[1,SGS%nVert],[0,3],Label='SGS%MAW')

! COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO MAW:
SGS%MAW(SGS%MVSta:SGS%nVert,0:3) = SGS%DAW(SGS%MVSta:SGS%nVert,0:3)
! COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO MAW. HOWEVER,
!    NOTE THAT THE MAW TABLE IS ACCESSED BY THE UPPER VERTEX.
SGS%MAW(1:SGS%MVSta-1,0:3) = 0
do IU=1,SGS%MVSta-1
  do IC=0,3
    ID = SGS%Down(IU,IC)
    if (ID /= 0) SGS%MAW(IU,IC) = SGS%RAW(ID,IC)
  end do
end do
! FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
ISUM = 1
do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    IU = SGS%Up(IV,IC)
    if (IU == 0) cycle
    SGS%MAW(IU,IC) = ISUM+SGS%MAW(IU,IC)
  end do
  ISUM = ISUM+SGS%RAW(IV,4)
end do
do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    if (SGS%Down(IV,IC) == 0) cycle
    SGS%MAW(IV,IC) = ISUM+SGS%MAW(IV,IC)
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
!    ISGT=1,..,26, WHOSE TOP LEFT VERTEX IS IVLT. ISGM GIVES
!    ZERO IF THE SEGMENT IS IMPOSSIBLE IN THE GRAPH DEFINED BY
!    THE PALDUS TABLE DRT, ELSE IT IS THE BOTTOM LEFT VERTEX
!    NUMBER OF THE SEGMENT. THE SEGMENT VALUE IS THEN VSGM.

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IA, IAL, IB, IBL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB, IVLT, IVRB, IVRT, LEV, MV, MVLL, MVRR
real(kind=wp) :: V
integer(kind=iwp), parameter :: IATAB = 3, IBTAB = 4
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ICL, ICR, ID
character(len=20) :: TEXT
character(len=*), parameter :: FRML(7) = ['        1.0         ','       -1.0         ','        1/(B+1)     ', &
                                          '       -1/(B+1)     ','  SQRT(   B /(B+1)) ','  SQRT((B+2)/(B+1)) ', &
                                          '  SQRT(B(B+2))/(B+1)']
#endif

call mma_allocate(CIS%IVR,SGS%nVert,2,Label='CIS%IVR')
call mma_allocate(CIS%ISGM,SGS%nVert,26,Label='CIS%ISGM')
call mma_allocate(CIS%VSGM,SGS%nVert,26,Label='CIS%VSGM')
call mma_allocate(EXS%MVL,CIS%nMidV,2,Label='EXS%MVL')
call mma_allocate(EXS%MVR,CIS%nMidV,2,Label='EXS%MVR')

CIS%IVR(:,:) = 0
do LEV=1,SGS%nLev
  IV1 = SGS%LTV(LEV)
  IV2 = SGS%LTV(LEV-1)-1
  do IVL=IV1,IV2
    IAL = SGS%DRT(IVL,IATAB)
    IBL = SGS%DRT(IVL,IBTAB)
    do IV=IVL+1,IV2
      IA = SGS%DRT(IV,IATAB)
      if (IA == IAL) then
        IB = SGS%DRT(IV,IBTAB)
        if (IB == (IBL-1)) CIS%IVR(IVL,1) = IV
      else if (IA == (IAL-1)) then
        IB = SGS%DRT(IV,IBTAB)
        if (IB == (IBL+1)) CIS%IVR(IVL,2) = IV
      end if
    end do
  end do
end do
! CONSTRUCT THE MVL AND MVR TABLES:
do IVL=SGS%MVSta,SGS%MVEnd
  MVLL = IVL-SGS%MVSta+1
  MVRR = 0
  if (CIS%IVR(IVL,1) /= 0) MVRR = CIS%IVR(IVL,1)-SGS%MVSta+1
  EXS%MVR(MVLL,1) = MVRR
  MVRR = 0
  if (CIS%IVR(IVL,2) /= 0) MVRR = CIS%IVR(IVL,2)-SGS%MVSta+1
  EXS%MVR(MVLL,2) = MVRR
  EXS%MVL(MVLL,1) = 0
  EXS%MVL(MVLL,2) = 0
end do
do MV=1,CIS%nMidV
  if (EXS%MVR(MV,1) /= 0) EXS%MVL(EXS%MVR(MV,1),1) = MV
  if (EXS%MVR(MV,2) /= 0) EXS%MVL(EXS%MVR(MV,2),2) = MV
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' MIDVERT PAIR TABLES MVL,MVR IN MKSEG:'
write(u6,*) ' MVL TABLE:'
write(u6,1234) (MV,EXS%MVL(MV,1),EXS%MVL(MV,2),MV=1,CIS%nMidV)
write(u6,*) ' MVR TABLE:'
write(u6,1234) (MV,EXS%MVR(MV,1),EXS%MVR(MV,2),MV=1,CIS%nMidV)
write(u6,*)
write(u6,*) ' VERTEX PAIR TABLE IVR IN MKSEG:'
write(u6,1234) (IVL,CIS%IVR(IVL,1),CIS%IVR(IVL,2),IVL=1,SGS%nVert)
#endif

! INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
CIS%ISGM(:,:) = 0
CIS%VSGM(:,:) = Zero

do IVLT=1,SGS%nVert
  do ISGT=1,26
    ITT = ITVPT(ISGT)
    IVRT = IVLT
    if ((ITT == 1) .or. (ITT == 2)) IVRT = CIS%IVR(IVLT,ITT)
    if (IVRT == 0) cycle
    IVLB = SGS%Down(IVLT,IC1(ISGT))
    if (IVLB == 0) cycle
    IVRB = SGS%Down(IVRT,IC2(ISGT))
    if (IVRB == 0) cycle
    ! SEGMENT IS NOW ACCEPTED AS POSSIBLE.

    CIS%ISGM(IVLT,ISGT) = IVLB
    IB = SGS%DRT(IVLT,IBTAB)
    select case (ISVC(ISGT))
      case (1)
        V = One
      case (2)
        V = -One
      case (3)
        V = One/real(1+IB,kind=wp)
      case (4)
        V = -One/real(1+IB,kind=wp)
      case (5)
        V = sqrt(real(IB,kind=wp)/real(1+IB,kind=wp))
      case (6)
        V = sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
      case (7)
        V = sqrt(real(IB*(2+IB),kind=wp))/real(1+IB,kind=wp)
      case default
        V = Zero ! Dummy assignment
        call Abend()
    end select
    CIS%VSGM(IVLT,ISGT) = V
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' SEGMENT TABLE IN MKSEG.'
write(u6,*) ' VLT SGT ICL ICR VLB       SEGMENT TYPE         FORMULA'
do IV=1,SGS%nVert
  do ISGT=1,26
    ID = CIS%ISGM(IV,ISGT)
    if (ID == 0) cycle
    ICL = IC1(ISGT)
    ICR = IC2(ISGT)
    if (ISGT <= 4) TEXT = '  WALK SECTION.'
    if ((ISGT >= 5) .and. (ISGT <= 8)) TEXT = ' TOP SEGMENT.'
    if ((ISGT >= 9) .and. (ISGT <= 18)) TEXT = ' MID-SEGMENT.'
    if ((ISGT >= 19) .and. (ISGT <= 22)) TEXT = ' BOTTOM SEGMENT.'
    if (ISGT > 22) TEXT = ' DOWN-WALK SECTION.'
    write(u6,2345) IV,ISGT,ICL,ICR,ID,TEXT,FRML(ISVC(ISGT))
  end do
end do
1234 format(3(3(1X,I4),4X))
2345 format(1X,5I4,5X,A20,5X,A20)
# endif

end subroutine MKSEG

subroutine NRCOUP(SGS,CIS,EXS)

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IBSYM, ICL, INDEO, INDEOB, INDEOT, IP, IPQ, IQ, ISGT, ISYDS1, ISYM, ISYUS1, ITSYM, IVLB, IVLT, LEV, LFTSYM, &
                     MV, MV1, MV2, MV3, MV4, MV5, MXDWN, MXUP, N, NDWNS1, NSGMX, NSGTMP, NT1TMP, NT2TMP, NT3TMP, NT4TMP, NT5TMP, &
                     NUPS1, NUW
integer(kind=iwp), allocatable :: NRL(:,:,:)
integer(kind=iwp), parameter :: LTAB = 1
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IS, IST, NCP
#endif

call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW',safe='*')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW',safe='*')
call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF',safe='*')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF',safe='*')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF',safe='*')

EXS%MxEO = (SGS%nLev*(SGS%nLev+5))/2
call mma_allocate(EXS%NOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
call mma_allocate(EXS%IOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
call mma_allocate(NRL,[1,SGS%nSym],[1,SGS%nVert],[0,EXS%MxEO],Label='NRL')

NRL(:,1:SGS%MVEnd,:) = 0
NRL(1,1,0) = 1

do IVLT=1,SGS%MVSta-1
  LEV = SGS%DRT(IVLT,LTAB)
  do ISGT=1,26
    IVLB = CIS%ISGM(IVLT,ISGT)
    if (IVLB == 0) cycle
    ICL = IC1(ISGT)
    ISYM = 1
    if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
    do ITSYM=1,SGS%nSym
      IBSYM = Mul(ITSYM,ISYM)

      select case (ISGT)
        case (:4)
          ! THIS IS AN UPPER WALK.
          NRL(IBSYM,IVLB,0) = NRL(IBSYM,IVLB,0)+NRL(ITSYM,IVLT,0)
        case (5:8)
          ! THIS IS AN TOP SEGMENT.
          INDEO = LEV+(IBVPT(ISGT)-1)*SGS%nLev
          NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,0)
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
            INDEOB = IPQ+2*SGS%nLev
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

MXUP = 0
do MV=1,CIS%nMidV
  IVLT = MV+SGS%MVSta-1
  do LFTSYM=1,SGS%nSym
    CIS%NOW(1,LFTSYM,MV) = NRL(LFTSYM,IVLT,0)
    MXUP = max(MXUP,CIS%NOW(1,LFTSYM,MV))
    do INDEO=1,EXS%MxEO
      EXS%NOCP(INDEO,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO)
    end do
  end do
end do

NRL(:,SGS%MVSta:,:) = 0
NRL(1,SGS%nVert,0) = 1
do IVLT=SGS%nVert-1,SGS%MVSta,-1
  LEV = SGS%DRT(IVLT,LTAB)
  do ISGT=1,26
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
      EXS%nICOup = EXS%nICOup+EXS%NOCP(INDEO,LFTSYM,MV)
    end do
  end do
end do

call CSFCOUNT(CIS,SGS%nSym,NUW)

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
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
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

end subroutine NRCOUP


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
! essentially a list of which one of segments nr 1..26 that are
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
! ISGPTH(ISEG ,LEV)=Segment type (1..26).
! These indices are used to denote the columns of table ISGPTH.

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: i, i1, i2, IAWS, IC, ICL, ICOP, ICR, IHALF, iLnd, IndEO, iP, iPos, iQ, iS, iSg, iSgt, iSym, iT, iTyp, iTypMx, &
                     iTypT, iVlb, iVlt, iVrt, iVrTop, iVTab, iVTEnd, iVTop, iVTSta, L, Lev, Lev1, Lev2, LftSym, LL, MV, nCheck, &
                     nVTab_Final
real(kind=wp) :: C
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I3, ICOP1, ICOP2, ICP1, ICP2, N, NRC, NRCPQ
real(kind=wp) :: CP
#endif
integer(kind=iwp), allocatable :: ILNDW(:), ISGPTH(:,:)
real(kind=wp), allocatable :: val(:), VTab(:)
integer(kind=iwp), parameter :: IVLFT = 1, ITYPE = 2, IAWSL = 3, IAWSR = 4, ILS = 5, ICS = 6, ISEG = 7, nVTab = 5000

call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
call mma_allocate(ISGPTH,[1,7],[0,SGS%nLev],Label='ISGPTH')
call mma_allocate(val,[0,SGS%nLev],Label='val')
call mma_allocate(VTab,nVTab,Label='VTab')

call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')

!CIS%nIpWlk = 1+(SGS%MidLev-1)/15
!CIS%nIpWlk = max(CIS%nIpWlk,1+(nLev-SGS%MidLev-1)/15)
! NOW IS ZEROED AND WILL BE USED AS AN ARRAY OF COUNTERS, BUT WILL BE RESTORED FINALLY.
do IHALF=1,2
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      CIS%NOW(IHALF,IS,MV) = 0
    end do
  end do
end do
! SIMILAR FOR THE COUPLING COEFFICIENT TABLE:
do INDEO=1,EXS%MxEO
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      EXS%NOCP(INDEO,IS,MV) = 0
    end do
  end do
end do

! COUPLING COEFFICIENT VALUE TABLE:
NVTAB_FINAL = 2
VTab(1) = One
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
  do IVTOP=IVTSTA,IVTEND
    do ITYP=0,ITYPMX
      IVRTOP = IVTOP
      if (ITYP > 0) IVRTOP = CIS%IVR(IVTOP,ITYP)
      if (IVRTOP == 0) cycle
      LEV = LEV1
      ISGPTH(IVLFT,LEV) = IVTOP
      ISGPTH(ITYPE,LEV) = ITYP
      ISGPTH(IAWSL,LEV) = 0
      ISGPTH(IAWSR,LEV) = 0
      ISGPTH(ILS,LEV) = 1
      ISGPTH(ISEG,LEV) = 0
      val(LEV) = One
      do while (LEV <= LEV1)
        ITYPT = ISGPTH(ITYPE,LEV)
        IVLT = ISGPTH(IVLFT,LEV)
        do ISGT=ISGPTH(ISEG,LEV)+1,26
          IVLB = CIS%ISGM(IVLT,ISGT)
          if (IVLB == 0) cycle
          if (ITYPT == ITVPT(ISGT)) exit
        end do
        if (ISGT > 26) then
          ISGPTH(ISEG,LEV) = 0
          LEV = LEV+1
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
            do LL=LEV2+1,LEV1,15
              IC = 0
              do L=min(LL+14,LEV1),LL,-1
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
        if (JPOS == 16) then
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
        if (JPOS == 16) then
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

subroutine SG_Free(SGS,CIS,EXS)
! PURPOSE: FREE THE SGUGA TABLES

#include "intent.fh"

implicit none
type(SGStruct), intent(_OUT_) :: SGS
type(CIStruct), intent(_OUT_) :: CIS
type(EXStruct), optional, intent(_OUT_) :: EXS

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
call mma_deallocate(SGS%SCR,safe='*')
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

If (Present(EXS)) THEN
   call mma_deallocate(EXS%NOCP,safe='*')
   call mma_deallocate(EXS%IOCP,safe='*')
   call mma_deallocate(EXS%ICoup,safe='*')
   call mma_deallocate(EXS%VTab,safe='*')
   call mma_deallocate(EXS%SGTMP,safe='*')
   call mma_deallocate(EXS%MVL,safe='*')
   call mma_deallocate(EXS%MVR,safe='*')
   call mma_deallocate(EXS%USGN,safe='*')
   call mma_deallocate(EXS%LSGN,safe='*')
End If

end subroutine SG_Free

subroutine CSFCOUNT(CIS,NSYM,NUW)

implicit none
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp), intent(in) :: NSYM
integer(kind=iwp), intent(out) :: NUW
integer(kind=iwp) :: ISYDWN, ISYM, ISYTOT, ISYUP, MV, N

! CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
! SEPARATED FOR EACH MIDVERTEX AND SYMMETRY

NUW = 0
do MV=1,CIS%nMidV
  do ISYM=1,NSYM
    CIS%IOW(1,ISYM,MV) = NUW*CIS%nIpWlk
    NUW = NUW+CIS%NOW(1,ISYM,MV)
  end do
end do
CIS%nWalk = NUW
do MV=1,CIS%nMidV
  do ISYM=1,NSYM
    CIS%IOW(2,ISYM,MV) = CIS%nWalk*CIS%nIpWlk
    CIS%nWalk = CIS%nWalk+CIS%NOW(2,ISYM,MV)
  end do
end do

! CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
! SEPARATED BY MIDVERTICES AND SYMMETRY.
! FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.

CIS%NCSF(:) = 0
do ISYTOT=1,NSYM
  do MV=1,CIS%nMidV
    do ISYUP=1,NSYM
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

subroutine MkISM_RAW(SGS,nLev,XLevel,XL2Act)

type(SGStruct), target, intent(inout) :: SGS
integer(kind=iwp), intent(in):: nLev
integer(kind=iwp), optional, intent(in):: xLevel(MxLev), xL2Act(MxLev)
integer(kind=iwp) iq

Write (6,*) 'MkISM_RAW'
Write (6,*) 'NLEV=',NLEV
Write (6,*) ALLOCATED(SGS%ISM)
SGS%NLEV = nLEV
! Allocate Level to Symmetry table ISm:
call mma_allocate(SGS%ISM,SGS%nLev,Label='SGS%ISM')


! Initiate if not already set externally.
If (Present(XLevel).and.Present(XL2Act)) Then
   Level(1:MxLev)=xLevel(1:MxLev)
   L2Act(1:MxLev)=xL2Act(1:MxLev)
Else If (Present(XLevel)) Then
   Level(1:MxLev)=xLevel(1:MxLev)
Else If (Present(XL2Act)) Then
   L2Act(1:MxLev)=xL2Act(1:MxLev)
Else
   If (LEVEL(1)==0) LEVEL(1:SGS%nLev)=[(iq,iq=1,SGS%nLev)]
   If (L2Act(1)==0) L2Act(1:SGS%nLev)=[(iq,iq=1,SGS%nLev)]
End If

! Default to incremental index if not properly set
If (Level(1)==0) Level(1:SGS%nLev)=[(iq,iq=1,SGS%nLev)]
If (L2Act(1)==0) L2Act(1:SGS%nLev)=[(iq,iq=1,SGS%nLev)]

End subroutine MkISM_RAW

  subroutine MKISM_RASSCF(SGS,xnLev,xNSM)
  ! PURPOSE: CREATE THE SYMMETRY INDEX VECTOR

    use gas_data, only: NGAS, NGSSH
    use rasscf_global, only: NSM

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp), optional, intent(in) :: xnLev, xNSM(MxLev)

    integer(kind=iwp) :: IGAS, ISYM, NLEV, NSTA

!#define _OLD_
#ifdef _OLD_
    NLEV = 0
    do IGAS=1,NGAS
       do ISYM=1,SGS%NSYM
          NSTA = NLEV+1
          NLEV = NLEV+NGSSH(IGAS,ISYM)
          NSM(NSTA:NLEV) = ISYM
       end do
    end do

    Call MkISM_RAW(SGS,nLev)

    If (Present(xnLev)) Then
       SGS%ISM(1:SGS%nLev) = xNSM(1:SGS%nLev)
    Else
       SGS%ISM(1:SGS%nLev) = NSM(1:SGS%nLev)
    End If
#else
    If (Present(xnLev)) Then
       Call MkISM_RAW(SGS,xnLev)
       SGS%ISM(1:SGS%nLev) = xNSM(1:SGS%nLev)
    Else
       NLEV = 0
       do IGAS=1,NGAS
         do ISYM=1,SGS%NSYM
           NSTA = NLEV+1
           NLEV = NLEV+NGSSH(IGAS,ISYM)
           NSM(NSTA:NLEV) = ISYM
         end do
       end do
       Call MkISM_RAW(SGS,nLev)
       SGS%ISM(1:SGS%nLev) = NSM(1:SGS%nLev)
    End If
#endif

  end subroutine MKISM_RASSCF

  subroutine MKISM(SGS,xnLev,xNSM)

  use UnixInfo, only: ProgName

  type(SGStruct), target, intent(inout) :: SGS
  integer(iwp), optional, intent(in):: xnLev, xNSM(MxLev)

    Select Case(ProgName(1:6))
    Case ('rassi ')
      call mkISM_Rassi(SGS)
    Case ('mclr  ')
      If (Present(xnLev)) Then
         Write (6,*) 'xNLEV=',xNLEV
         Write (6,*) 'xNSM=',xNSM
      Else
         Write (6,*) 'Something is missing'
      End If
      call mkISM_mclr(SGS)
    Case ('caspt2')
      call mkISM_cp2(SGS)
    Case ('rasscf','casvb ')
      If (Present(xnLev)) Then
         Write (6,*) 'MKISM: xnLev=',xnLev
         Write (6,*) 'MKISM: xNSM(1:5)=',xNSM(1:5)
         call mkISM_rasscf(SGS,xnLev,xNSM)
      Else
         call mkISM_rasscf(SGS)
      End If
    Case Default
       Write (u6,*) 'MkISM: not setup for program:', ProgName
       Call Abend()
    End Select

  end subroutine MKISM

  subroutine MKISM_MCLR(SGS)

    use input_mclr, only: NRS1, NRS2, NRS3, NLEV=>NTASH

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: iBas, iOrb, iSym

    Call MkISM_RAW(SGS,nLev)

    iOrb = 0
    do iSym=1,SGS%nSym
      do iBas=1,nRs1(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do
    do iSym=1,SGS%nSym
      do iBas=1,nRs2(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do
    do iSym=1,SGS%nSym
      do iBas=1,nRs3(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do

    Write (6,*) 'MKISM_MCLR:', SGS%ISM(:)

  end subroutine MKISM_MCLR

  subroutine MKISM_RASSI(SGS)

   use rassi_data, only: NASH, NLEV=>NASHT

   type(SGStruct), target, intent(inout) :: SGS
   integer(kind=iwp) :: iOrb, ISYM, IT, ILEV

    Call MkISM_RAW(SGS,nLev)

    iOrb = 0
    do ISYM=1,SGS%NSYM
      do IT=1,NASH(ISYM)
        iOrb = iOrb+1
        ILEV = LEVEL(iOrb)
        SGS%ISM(ILEV) = ISYM
      end do
    end do

  end subroutine MKISM_RASSI

  subroutine mkism_cp2(SGS)

    use caspt2_module, only: nAsh, NLEV=>nAshT

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: ILEV, ISYM, IT, iOrb

    Call MkISM_RAW(SGS,nLev)

    ! PAM060612: With true RAS space, the orbitals must be ordered
    ! first by RAS type, then by symmetry.

    iOrb = 0
    do ISYM=1,SGS%NSYM
      do IT=1,NASH(ISYM)
        iOrb = iOrb+1
        ILEV = LEVEL(iOrb)
        SGS%ISM(ILEV) = ISYM
      end do
    end do

  end subroutine mkism_cp2

end module SGUGA
