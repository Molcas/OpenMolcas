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

integer(kind=iwp) :: L2ACT(MXLEV), LEVEL(MXLEV)

public :: CIS, CIStruct, EXS, EXStruct, L2ACT, LEVEL, SGS, SGStruct


public :: SGINIT, MKSGUGA, MkCOT, MkCList

contains

subroutine MKSGUGA(SGS,CIS)
! PURPOSE: MAKE THE GUGA TABLES
! NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!          THE START ADDRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!          THREE USER DEFINED TYPES. Consult the gugx module for the details.

  type(SGStruct), target, intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2, ATAB = 3, BTAB = 4, CTAB = 5

  ! CREATE THE SYMMETRY INDEX VECTOR

  call MKISM(SGS)

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

  subroutine MKISM(SGS)

  use UnixInfo, only: ProgName

  type(SGStruct), target, intent(inout) :: SGS

    if (ProgName(1:6) == 'rassi') then
      call mkISm_Rassi(SGS)
    else if (ProgName(1:4) == 'mclr') then
      call mkISm_mclr(SGS)
    else if (ProgName(1:6) == 'caspt2') then
      call mkISm_cp2(SGS)
    else
      call mkNSM(SGS)
    end if

  end subroutine MKISM

  subroutine MKISM_MCLR(SGS)

    use input_mclr, only: NRS1, NRS2, NRS3, NSYM, NTASH

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: iBas, iOrb, iSym

    SGS%NLEV = ntASh

    call mma_allocate(SGS%ISM,SGS%nLev,Label='SGS%ISM')

    iOrb = 0
    do iSym=1,nSym
      do iBas=1,nRs1(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do
    do iSym=1,nSym
      do iBas=1,nRs2(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do
    do iSym=1,nSym
      do iBas=1,nRs3(iSym)
        iOrb = iOrb+1
        SGS%ISM(iOrb) = iSym
      end do
    end do

  end subroutine MKISM_MCLR

  subroutine MKISM_RASSI(SGS)

   use rassi_data, only: NASH, NASHT

   type(SGStruct), target, intent(inout) :: SGS
   integer(kind=iwp) :: ITABS, ISYM, IT, ILEV, nSym

    nSym = SGS%nSym
    SGS%NLEV = NASHT ! Total number of active orbitals
    ! Allocate Level to Symmetry table ISm:
    call mma_allocate(SGS%ISm,SGS%nLev,Label='SGS%ISm')
    ITABS = 0
    do ISYM=1,NSYM
      do IT=1,NASH(ISYM)
        ITABS = ITABS+1
        ILEV = LEVEL(ITABS)
        SGS%ISM(ILEV) = ISYM
      end do
    end do

  end subroutine MKISM_RASSI

  subroutine mkism_cp2(SGS)

    use fciqmc_interface, only: DoFCIQMC
    use caspt2_module, only: DoCumulant, nAsh, nAshT, nSym

    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: ILEV, iq, ISYM, IT, ITABS, nLev

    NLEV = NASHT
    SGS%nLev = NLEV
    call mma_allocate(SGS%ISM,NLEV,Label='ISM')
    ! ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
    ! PAM060612: With true RAS space, the orbitals must be ordered
    ! first by RAS type, then by symmetry.
    ITABS = 0
    do ISYM=1,NSYM
      do IT=1,NASH(ISYM)
        ITABS = ITABS+1
        ! Quan: Bug in LEVEL(ITABS) and L2ACT
        if (DoCumulant .or. DoFCIQMC) then
          do iq=1,NLEV
            LEVEL(iq) = iq
            L2ACT(iq) = iq
          end do
        end if
        ILEV = LEVEL(ITABS)
        SGS%ISM(ILEV) = ISYM
      end do
    end do

  end subroutine mkism_cp2

  subroutine MKNSM(SGS)
  ! PURPOSE: CREATE THE SYMMETRY INDEX VECTOR

    use gas_data, only: NGAS, NGSSH
    use rasscf_global, only: NSM
    use general_data, only: NSYM

    ! to get some dimensions
    type(SGStruct), target, intent(inout) :: SGS
    integer(kind=iwp) :: IGAS, ISYM, NLEV, NSTA

    NLEV = 0
    do IGAS=1,NGAS
      do ISYM=1,NSYM
        NSTA = NLEV+1
        NLEV = NLEV+NGSSH(IGAS,ISYM)
        NSM(NSTA:NLEV) = ISYM
      end do
    end do

    if (SGS%nSym /= 0) then
      SGS%nLev = nLev
      call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')
      SGS%ISM(1:nLev) = NSM(1:nLev)
    end if

  end subroutine MKNSM

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

SUBROUTINE SGINIT(nSym,nActEl,iSpin,                    &
                  SGS,CIS,EXS,                          &
                  nHole1,nEle3,nRas1T,nRas2T,nRas3T)
IMPLICIT None
Integer(kind=iwp), intent(in):: nSym,nActEl,iSpin
Type(SGStruct), intent(inout):: SGS
Type(CIStruct), intent(inout):: CIS
Integer(kind=iwp), optional, intent(in):: nHole1,nEle3,nRas1T,nRas2T,nRas3T
Type(EXStruct),  optional, intent(inout):: EXS

SGS%nSym=nSym
SGS%iSpin=iSpin
SGS%nActEl=nActEl

If (Present(EXS)) THEN
   SGS%LV1RAS=NRAS1T
   SGS%LV3RAS=nRas1T+NRAS2T
   SGS%LM1RAS=2*nRas1T-NHOLE1
   SGS%LM3RAS=NACTEL-NELE3
   IF ((NRAS1T+NRAS3T)/=0) Then
      SGS%IFRAS=1
   Else
      SGS%IFRAS=0
   End If
End IF

Call MkSGUGA(SGS,CIS)

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
End IF

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

CALL MKMAW(SGS)

If (Present(EXS)) THEN
! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

CALL MKSEG(SGS,CIS,EXS)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.

CALL NRCOUP(SGS,CIS,EXS)

CALL MKCOUP(SGS,CIS,EXS)
End If

CALL mma_deallocate(SGS%DAW)
CALL mma_deallocate(SGS%RAW)

If (Present(EXS)) THEN
   Call mma_deallocate(CIS%ISGM)
   Call mma_deallocate(CIS%VSGM)
   Call mma_deallocate(CIS%IVR)

   Call mma_deallocate(SGS%MAW)

   CALL mma_deallocate(SGS%DRT)
   Call mma_deallocate(SGS%DOWN)
   CALL mma_deallocate(SGS%UP)
   Call mma_deallocate(SGS%LTV)
End If

END SUBROUTINE SGINIT

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

end module SGUGA
