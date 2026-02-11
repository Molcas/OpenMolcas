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

module MkGUGA_mod

! This module contains a single function, to avoid explicit interfaces

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
private

public :: MKGUGA

contains

subroutine MKGUGA(SGS,CIS)
! PURPOSE: MAKE THE GUGA TABLES
! NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!          THE START ADDRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!          THREE USER DEFINED TYPES. Consult the gugx module for the details.

  use gugx, only: CIStruct, SGStruct

  type(SGStruct), target, intent(inout) :: SGS
  type(CIStruct), intent(inout) :: CIS
  integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2, ATAB = 3, BTAB = 4, CTAB = 5

  ! CREATE THE SYMMETRY INDEX VECTOR

  call MKISM()

  ! COMPUTE TOP ROW OF THE GUGA TABLE

  call mknVert0()

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

  call mkDRT0()

  ! IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
  ! VERTICES WHICH VIOLATE THE FORMER.

  if (SGS%IFRAS /= 0) then
    call mkRAS()

    ! REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)

    call mkDRT()

    ! IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED DRT TABLE

  end if

  nullify(SGS%DOWNP,SGS%DRTP)

  ! CALCULATE ARC WEIGHT.

  call MKDAW()

  ! COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS

  call MKRAW()

  ! COMPUTE LTV TABLES.

  call MKLTV()

  ! COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.

  call MKMID()

contains

  subroutine MKISM()

    use UnixInfo, only: ProgName

    if (ProgName(1:6) == 'rassi') then
      call mkISm_Rassi()
    else if (ProgName(1:4) == 'mclr') then
      call mkISm_mclr()
    else if (ProgName(1:6) == 'caspt2') then
      call mkISm_cp2()
    else
      call mkNSM()
    end if

  end subroutine MKISM

  subroutine MKISM_MCLR()

    use input_mclr, only: NRS1, NRS2, NRS3, NSYM, NTASH

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

  subroutine MKISM_RASSI()

    use gugx, only: LEVEL
    use rassi_data, only: NASH, NASHT

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

  subroutine mkism_cp2()

    use fciqmc_interface, only: DoFCIQMC
    use gugx, only: L2ACT, LEVEL
    use caspt2_module, only: DoCumulant, nAsh, nAshT, nSym

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

  subroutine MKNSM()
  ! PURPOSE: CREATE THE SYMMETRY INDEX VECTOR

    use gugx, only: SGS
    use gas_data, only: NGAS, NGSSH
    use rasscf_global, only: NSM
    use general_data, only: NSYM

    ! to get some dimensions
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

  subroutine mknVert0()

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

  subroutine mkDRT0()
  ! PURPOSE: CONSTRUCT THE UNRESTRICTED GUGA TABLE

    use Index_Functions, only: nTri_Elem1

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

  subroutine mkRAS()

    use UnixInfo, only: ProgName

    if (ProgName(1:5) == 'rassi') then
      call rmvert(SGS)
    else
      call RESTR(SGS)
    end if

  end subroutine mkRAS

  subroutine mkDRT()
  ! PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
  !          THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
  !          VERTICES FROM THE OLD TO THE NEW DRT TABLE

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

  subroutine MKDAW()
  ! PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE

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

  subroutine MKRAW()
  ! PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS

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

  subroutine MKLTV()
  ! PURPOSE: FIND THE MIDLEVEL

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

  subroutine MKMID()
  ! PURPOSE: FIND THE MIDLEVEL

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

end subroutine MKGUGA

end module MkGUGA_mod
