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

subroutine Input_Grid_It(iRun,INPORB,iReturn)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
!                                                                      *
! Object: input module for grid code                                   *
!                                                                      *
!***********************************************************************

use grid_it_globals, only: CutOff, GridAxis1, GridAxis2, GridAxis3, GridDense, GridNormal, GridOrigin, GridSparse, iGauss, &
                           iGridNpt, iMaxDown, iMaxUp, imoPack, ipGrid, iReq, isAll, isAtom, isAuMO, isBinary, isColor, isCurDens, &
                           isCutOff, isDebug, isDensity, isDerivative, isLine, isLuscus, isSphere, isTheOne, isTotal, isUserGrid, &
                           isVirt, isXField, itRange, MAXGRID, nBytesPackedVal, nGridPoints, NoOrb, NoSort, nReq, OneCoor, Region, &
                           TheGap, TheName, Title1, Virt
use Constants, only: Zero, Four, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iRun, iReturn
character(len=*), intent(out) :: INPORB
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i_bits, iCustOrig, iD, iErr, iKey, inUnit, isAuto, isFileOrb, isNet, isReadNet, iSubBlock, iTemp(3), j, &
                     magicValue, n
real(kind=wp) :: dTemp(4), rD, rSubBlock, SubBlock(3), XFminCh, xHiErr, xLeft, xLoErr, xRight
character(len=265) :: AllKeys
character(len=256) :: FileStr, FileIn
character(len=120) :: SelectStr
character(len=80) :: Key, MULLprt
integer(kind=iwp), external :: MyGetKey

AllKeys = 'PRIN BINA ASCI NPOI DENS SPAR ORBI REGI ONE  TITL GAP  END  NODE TOTA NAME VB   ALL  ATOM CUBE GRID '// &
          'PACK PKLI PKBI NOOR LINE ORAN ERAN DEBU CUTO NOPA GORI SELE NOSO FILE SPHR COLO VIRT MULL SUBB XDER '// &
          'YDER ZDER GDER CURD CRXJ UMAX NOLU XFIE LUS1 LUS2 PLUS MINU XFMI'

!do i=1,nRout
!  nPrint(i) = 5
!end do
isBinary = 3
isAuto = 1
isNet = 0
isReadNet = 0
isTheOne = 0
Title1 = ' '
TheName = ' '
TheGap = Four
!nMOmin = 0
!nGrid = 5
nReq = -1
isDensity = 1
isDerivative = 0
isCurDens = 0
!isRxJ = 0
!iuseMaxes = 0
isAuMO = -1
isAtom = 0
isTotal = 0
!isVB = 0
isAll = 0
isUserGrid = 0
iGauss = 0
NoOrb = 0
iMaxUp = 7
iMaxDown = 7
isLine = 0
itRange = 1
isDebug = 0
isCutOff = 0
CutOff = 2.5_wp
iCustOrig = 0
NoSort = 0
isFileOrb = 0
isSphere = 0
isColor = 0
isVirt = 0
!isMULL = 0
!isLONGPRT = 0
!isWDW = 0
iSubBlock = 0
isLuscus = 1
!isLusMath = 0
!aLusMath = -1
! not preparing the GRIDCHARGE file as external source for XFIELD input
isXField = 0
XFminCh = Zero
! Default values for packing
imoPack = 0
!isBinPack = 0
xLeft = 0.005_wp
xRight = 0.7_wp
! (really, half range:)
!iyRange = 128
nBytesPackedVal = 1
!iMinYLeft = 4
xLoErr = 0.10_wp
xHiErr = Quart
INPORB = 'INPORB'
if (iRun == 0) then
  ! make defaults for a fast run via call from other module
  isNet = 1  ! set sparse
  !isBinary = 0 ! temporary set Ascii output
  !iMaxUp = 1
  !iMaxDown = 5
  !imoPack = 0 ! packed grids not working with 64bit (?)
  !isCutOff = 1
  !goto 500 (?)
end if

! KeyWord directed input

InUnit = 5
! Function MyGetKey (InUnit, What, IValue, RValue, SValue, N, IArray, RArray)

call RdNLst(InUnit,'GRID_IT')
do
  if (MyGetKey(InUnit,'S',iD,rD,Key,iD,[iD],[rD]) /= 0) exit
  iKey = index(AllKeys,Key(1:4))
  if ((iKey == 0) .or. ((iKey-1)/5*5 /= (iKey-1))) then
    write(u6,'(a,a)') 'Unrecognized keyword in input file:',Key(1:4)
    call Quit_OnUserError()
  end if
  iKey = (iKey-1)/5+1

  select case (iKey)
    case (1)
      ! PRIN
      if (MyGetKey(InUnit,'I',n,rD,Key,iD,[iD],[rD]) /= 0) call error()
      do j=1,n
        if (MyGetKey(InUnit,'A',iD,rD,Key,2,iTemp,[rD]) /= 0) call error()
        !write(u6,*) 'debug'
        !nPrint(iTemp(1)) = iTemp(2)
      end do
    case (2)
      ! BINARY = default
      isBinary = 1
    case (3)
      ! ASCII = for debug
      !write(u6,*) ' Keyword ASCII is obsolete'
      !write(u6,*) ' It can be used only for debugging purpose'
      !write(u6,*) ' Note that .lus files produced with this option '
      !write(u6,*) '      can not be visualised'
      isBinary = 0
    case (4)
      ! NPOI
      if (MyGetKey(InUnit,'A',iD,rD,Key,3,iGridNpt,[rD]) /= 0) call error()
      isNet = -1
      isReadNet = isReadNet+1
    case (5)
      ! DENSE - dense grid network..
      isNet = 2
      isReadNet = isReadNet+1
    case (6)
      ! SPARSE - rare grid network..
      isNet = 1
      isReadNet = isReadNet+1
    case (7)
      ! ORBI Orbitals
      if (nReq > 0) then
        write(u6,*) 'ORBI keyword can not be used together with SELEct'
        call Quit_OnUserError()
      end if
      if (MyGetKey(InUnit,'I',nReq,rD,Key,iD,[iD],[rD]) /= 0) call error()

      if (nReq > MAXGRID) then
        write(u6,'(a,i5,a,i5)') 'Too many requested orbitals ',nReq,'>',MAXGRID
        call Quit_OnUserError()
      end if
      read(inUnit,*,iostat=iErr) (iReq(i),i=1,nReq*2)
      if (iErr /= 0) call error()
      !if (MyGetKey(InUnit,'A',iD,rD,Key,nReq*2,iReq,[rD]) /= 0) call error()
      isAuMO = 0
    case (8)
      ! REGION
      if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) call error()
      itRange = 1
      isAuMO = 1
      write(u6,*) ' *** Warning keyword REGION is obsolete'
      write(u6,*) ' ***         assumimg Energy range '
    case (9)
      ! ONE - debug option
      if (MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor) /= 0) call error()
      isTheOne = 1
      isBinary = 0
    case (10)
      ! TITLE
      ! NOTE: Title can be only ONE line here!!!
      if (MyGetKey(InUnit,'S',iD,rD,Title1,iD,[iD],[rD]) /= 0) call error()
    case (11)
      ! GAP
      if (MyGetKey(InUnit,'R',iD,TheGap,Key,iD,[iD],[rD]) /= 0) call error()
    case (12)
      ! END
      exit
    case (13)
      ! NODENSITY
      isDensity = 0
    case (14)
      ! TOTAL
      isTotal = 1
    case (15)
      ! NAME
      read(InUnit,'(a)') TheName
      !if (MyGetKey(InUnit,'S',iD,rD,TheName,iD,[iD],[rD]) /= 0) call error()
      ! unfortunately MyGetKey uppercases strings!
    case (16)
      ! VB
      !isVB = 1
      call Quit_OnUserError()
    case (17)
      ! All
      isAll = 1
    case (18)
      ! Atom
      isAtom = 1
      isNet = -1
      iGridNpt(1) = 0
      iGridNpt(2) = 0
      iGridNpt(3) = 0
    case (19)
      ! CUBE
      iGauss = 1
      isBinary = 0
      write(u6,*) 'Cube option is moved to grid2cube'
      call Quit_OnUserError()
    case (20)
      ! Grid
      isUserGrid = 1
      isBinary = 0
      isNet = -1
      iGridNpt(1) = 0
      iGridNpt(2) = 0
      iGridNpt(3) = 0
      if (MyGetKey(InUnit,'I',nGridPoints,rD,Key,iD,[iD],[rD]) /= 0) call error()
      call GetMem('Grid','ALLO','REAL',ipGrid,nGridPoints*3)
      read(InUnit,*,iostat=iErr) (Work(ipGrid+i-1),i=1,nGridPoints*3)
      if (iErr /= 0) call error()
    case (21)
      ! Pack
      !imoPack = 1
    case (22)
      ! PkLims
      if (MyGetKey(InUnit,'D',iD,rD,Key,4,[iD],dTemp) /= 0) call error()
      xLeft = dTemp(1)
      xRight = dTemp(2)
      xLoErr = dTemp(3)
      xHiErr = dTemp(4)
    case (23)
      ! PkBits
      if (MyGetKey(InUnit,'I',i_bits,rD,Key,iD,[iD],[rD]) /= 0) call error()
      if (i_bits == 16) then
        !iyRange = 32768
        nBytesPackedVal = 2
      end if
    case (24)
      ! NoOrbitals
      NoOrb = 1
    case (25)
      ! LINE - density on line
      if (MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor) /= 0) call error()
      isTheOne = 1
      isTotal = 1
      isBinary = 0
      isLine = 1
    case (26)
      ! ORANGE
      if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) call error()
      itRange = 0
      isAuMO = 1
      NoSort = 1
    case (27)
      ! ERANGE
      if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) call error()
      itRange = 1
      isAuMO = 1
    case (28)
      ! DEBUG
      isBinary = 0
      isDebug = 1
    case (29)
      ! CUTOFF
      if (MyGetKey(InUnit,'R',iD,CutOff,Key,iD,[iD],[rD]) /= 0) call error()
      isCutOff = 1
    case (30)
      ! NOPACK
      !imoPack = 0
    case (31)
      ! GORI
      iCustOrig = 1
      if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridOrigin) /= 0) call error()
      if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis1) /= 0) call error()
      if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis2) /= 0) call error()
      if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis3) /= 0) call error()
    case (32)
      ! SELEct
      if (MyGetKey(InUnit,'S',iD,rD,SelectStr,iD,[iD],[rD]) /= 0) call error()
      if (nReq > 0) then
        write(u6,*) 'SELEct keyword can not be used together with ORBItals'
        call Quit_OnUserError()
      end if
      call gridExpandSelect(SelectStr)
      isAuMO = 0
    case (33)
      ! NOSOrt
      NoSort = 1
    case (34)
      ! FILE
      read(InUnit,'(A)') FileIn
      isFileOrb = 1
      call fileorb(FileIn,FileStr)
      write(u6,*) 'INPORB file: ',trim(FileStr)
    case (35)
      ! SPHR
      isSphere = 1
    case (36)
      ! COLOr
      isColor = 1
    case (37)
      ! VIRT
      isVirt = 1
      if (MyGetKey(InUnit,'R',iD,Virt,Key,iD,[iD],[rD]) /= 0) call error()
    case (38)
      ! MULLiken charges per MO
      !isMULL = 1
      read(InUnit,'(A)') MULLPRT
      call upCASE(MULLPRT)
      call LeftAd(MULLPRT)
      !if (MULLPRT(1:4) == 'LONG') isLONGPRT = 1
    case (39)
      ! SUBBLOCK
      if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],SubBlock) /= 0) call error()
      if (MyGetKey(InUnit,'R',iD,rSubBlock,Key,iD,[iD],[rD]) /= 0) call error()
      iSubBlock = 1
    case (40)
      ! XDER
      isDerivative = 1
      call Quit_OnUserError()
    case (41)
      ! YDER
      isDerivative = 2
      call Quit_OnUserError()
    case (42)
      ! ZDER
      isDerivative = 3
      call Quit_OnUserError()
    case (43)
      isDerivative = 4
      call Quit_OnUserError()
    case (44)
      ! CURD (current density)
      isCurDens = 1
      call Quit_OnUserError()
    case (45)
      ! CRXJ (current density, rxj)
      isCurDens = 1
      !isRxJ = 1
      call Quit_OnUserError()
    case (46)
      ! UMAX (use magnetic axes)
      !iuseMaxes = 1
    case (47)
      ! NOLUSCUS
      isLuscus = 0
      isBinary = 0
    case (48)
      ! XFIEld - ask Grid_It to compute electronic density on a DFT integration grid
      isXField = 1
      isReadNet = isReadNet+1 !make the grid definition exclusive
    case (49)
      ! LUS1
      write(u6,*) 'Not implemented'
      !isLusMath = 1
      !read(InUnit,'(a)') LUS1
    case (50)
      ! LUS2
      write(u6,*) 'Not implemented'
      !isLusMath = 1
      !read(InUnit,'(a)') LUS2
    case (51)
      ! PLUS
      write(u6,*) 'Not implemented'
      !aLusMath = 1
    case (52)
      ! MINUS
      write(u6,*) 'Not implemented'
      !aLusMath = -1
    case (53)
      ! XFMI xfield minimum charge of each grid point to be stored
      if (MyGetKey(InUnit,'R',iD,XFminCh,Key,iD,[iD],[rD]) /= 0) call error()
    case default
  end select
end do

!***********************************************************************
!                                                                      *
!                       End of input section.                          *
!                                                                      *
!***********************************************************************
!if (isLusMath == 1) return
if ((isLuscus == 1) .and. (isBinary == 0)) then
  write(u6,*) 'ASCII keyword is set, but NoLUSCUS is not'
  write(u6,*) 'calling abend as the best option available'
  call Quit_OnUserError()
end if
close(InUnit)
if (isLuscus == 1) then
  if (isLine /= 0) then
    write(u6,*) 'LUSCUS and LINE options are not compatible'
    call Quit_OnUserError()
  end if
end if

if (isReadNet > 1) write(u6,'(a)') 'Warning: Double definition of GRID net'

! Well, there is something to do!

if (isFileOrb == 1) then
  INPORB = FileStr
end if
call OpenGrid(INPORB)

! try to generate grid position automatically

magicValue = 0
if (isNet >= 0) then
  magicValue = GridNormal
  if (isNet == 1) magicValue = GridSparse
  if (isNet == 2) magicValue = GridDense
end if
if (iCustOrig == 1) then
  if (isNet /= -1) then
    write(u6,*) 'GORI can be used only with NPOI'
    call Quit_OnUserError()
  end if
end if
call MyCoor(isAuto,GridOrigin(1),GridOrigin(2),GridOrigin(3),GridAxis1(1),GridAxis2(2),GridAxis3(3),iGridNpt(1),iGridNpt(2), &
            iGridNpt(3),magicValue,iCustOrig)

if (iSubBlock == 1) then
  do i=1,3
    GridOrigin(i) = SubBlock(i)-rSubBlock
  end do
  GridAxis1(1) = rSubBlock*2
  GridAxis2(2) = rSubBlock*2
  GridAxis3(3) = rSubBlock*2
  iGridNpt(1) = 40
  iGridNpt(2) = 40
  iGridNpt(3) = 40
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iReturn)

contains

subroutine error()
  write(u6,'(a,a,a)') 'Error during reading ',Key(1:20),'section in input file'
  call Quit_OnUserError()
end subroutine error

end subroutine Input_Grid_It
