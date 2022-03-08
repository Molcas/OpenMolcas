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

subroutine Input_Grid_It(iRun,INPORB)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
!                                                                      *
! Object: input module for grid code                                   *
!                                                                      *
!***********************************************************************

use grid_it_globals, only: CutOff, Grid, GridAxis1, GridAxis2, GridAxis3, GridDense, GridNormal, GridOrigin, GridSparse, iGauss, &
                           iGridNpt, iMaxDown, iMaxUp, iReq, isAll, isAtom, iAuMO, iBinary, isColor, isCurDens, isCutOff, isDebug, &
                           isDensity, iDerivative, isLine, isLuscus, isMOPack, isSphere, isTheOne, isTotal, isUserGrid, isVirt, &
                           itRange, MAXGRID, nBytesPackedVal, nGridPoints, NoOrb, NoSort, nReq, OneCoor, Region, TheGap, TheName, &
                           Title1, Virt
use stdalloc, only: mma_allocate
use Constants, only: Zero, Four
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), intent(in) :: iRun
character(len=*), intent(out) :: INPORB
integer(kind=iwp) :: i, i_bits, iD, iDum(1), iErr, iKey, inUnit, iNet, iReadNet, iTemp(3), j, magicValue, n
logical(kind=iwp) :: isAuto, isCustOrig, isFileOrb, isSubBlock
real(kind=wp) :: dTemp(4), rD, rDum(1), rSubBlock, SubBlock(3), XFminCh
character(len=256) :: FileStr, FileIn
character(len=120) :: SelectStr
character(len=80) :: Key, MULLprt
character :: What
character(len=*), parameter :: AllKeys = 'PRIN BINA ASCI NPOI DENS SPAR ORBI REGI ONE  TITL &
                                         &GAP  END  NODE TOTA NAME VB   ALL  ATOM CUBE GRID &
                                         &PACK PKLI PKBI NOOR LINE ORAN ERAN DEBU CUTO NOPA &
                                         &GORI SELE NOSO FILE SPHR COLO VIRT MULL SUBB XDER &
                                         &YDER ZDER GDER CURD CRXJ UMAX NOLU XFIE LUS1 LUS2 &
                                         &PLUS MINU XFMI'
integer(kind=iwp), external :: MyGetKey

!do i=1,nRout
!  nPrint(i) = 5
!end do
iBinary = 3
isAuto = .true.
iNet = 0
iReadNet = 0
isTheOne = .false.
Title1 = ' '
TheName = ' '
TheGap = Four
!nMOmin = 0
!nGrid = 5
nReq = -1
isDensity = .true.
iDerivative = 0
isCurDens = .false.
!isRxJ = 0
!iuseMaxes = 0
iAuMO = -1
isAtom = .false.
isTotal = .false.
!isVB = 0
isAll = .false.
isUserGrid = .false.
iGauss = 0
NoOrb = .false.
iMaxUp = 7
iMaxDown = 7
isLine = .false.
itRange = 1
isDebug = .false.
isCutOff = .false.
CutOff = 2.5_wp
isCustOrig = .false.
NoSort = .false.
isFileOrb = .false.
isSphere = .false.
isColor = .false.
isVirt = .false.
!isMULL = 0
!isLONGPRT = 0
!isWDW = 0
isSubBlock = .false.
isLuscus = .true.
!isLusMath = 0
!aLusMath = -1
! not preparing the GRIDCHARGE file as external source for XFIELD input
!isXField = 0
XFminCh = Zero
! Default values for packing
isMOPack = .false.
!isBinPack = 0
!xLeft = 0.005_wp
!xRight = 0.7_wp
! (really, half range:)
!iyRange = 128
nBytesPackedVal = 1
!iMinYLeft = 4
!xLoErr = 0.10_wp
!xHiErr = Quart
INPORB = 'INPORB'
if (iRun == 0) then
  ! make defaults for a fast run via call from other module
  iNet = 1  ! set sparse
  !iBinary = 0 ! temporary set Ascii output
  !iMaxUp = 1
  !iMaxDown = 5
  !isMOPack = .false. ! packed grids not working with 64bit (?)
  !isCutOff = .true.
  !goto 500 (?)
end if
GridOrigin(:) = Zero
GridAxis1(:) = Zero
GridAxis2(:) = Zero
GridAxis3(:) = Zero

! KeyWord directed input

InUnit = u5
! Function MyGetKey (InUnit, What, IValue, RValue, SValue, N, IArray, RArray)

call RdNLst(InUnit,'GRID_IT')
do
  What = 'S'
  if (MyGetKey(InUnit,What,iD,rD,Key,iD,iDum,rDum) /= 0) exit
  iKey = index(AllKeys,Key(1:4))
  if ((iKey == 0) .or. ((iKey-1)/5*5 /= (iKey-1))) then
    write(u6,'(a,a)') 'Unrecognized keyword in input file:',Key(1:4)
    call Quit_OnUserError()
  end if
  iKey = (iKey-1)/5+1

  select case (iKey)
    case (1)
      ! PRIN
      What = 'I'
      if (MyGetKey(InUnit,What,n,rD,Key,iD,iDum,rDum) /= 0) call error()
      do j=1,n
        What = 'A'
        if (MyGetKey(InUnit,What,iD,rD,Key,2,iTemp,rDum) /= 0) call error()
        !write(u6,*) 'debug'
        !nPrint(iTemp(1)) = iTemp(2)
      end do
    case (2)
      ! BINARY = default
      iBinary = 1
    case (3)
      ! ASCII = for debug
      !write(u6,*) ' Keyword ASCII is obsolete'
      !write(u6,*) ' It can be used only for debugging purpose'
      !write(u6,*) ' Note that .lus files produced with this option '
      !write(u6,*) '      can not be visualised'
      iBinary = 0
    case (4)
      ! NPOI
      What = 'A'
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iGridNpt,rDum) /= 0) call error()
      iNet = -1
      iReadNet = iReadNet+1
    case (5)
      ! DENSE - dense grid network..
      iNet = 2
      iReadNet = iReadNet+1
    case (6)
      ! SPARSE - rare grid network..
      iNet = 1
      iReadNet = iReadNet+1
    case (7)
      ! ORBI Orbitals
      if (nReq > 0) then
        write(u6,*) 'ORBI keyword can not be used together with SELEct'
        call Quit_OnUserError()
      end if
      What = 'I'
      if (MyGetKey(InUnit,What,nReq,rD,Key,iD,iDum,rDum) /= 0) call error()

      if (nReq > MAXGRID) then
        write(u6,'(a,i5,a,i5)') 'Too many requested orbitals ',nReq,'>',MAXGRID
        call Quit_OnUserError()
      end if
      read(inUnit,*,iostat=iErr) (iReq(i),i=1,nReq*2)
      if (iErr /= 0) call error()
      !What = 'A'
      !if (MyGetKey(InUnit,What,iD,rD,Key,nReq*2,iReq,rDum) /= 0) call error()
      iAuMO = 0
    case (8)
      ! REGION
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,2,iDum,Region) /= 0) call error()
      itRange = 1
      iAuMO = 1
      write(u6,*) ' *** Warning keyword REGION is obsolete'
      write(u6,*) ' ***         assumimg Energy range '
    case (9)
      ! ONE - debug option
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,7,iDum,OneCoor) /= 0) call error()
      isTheOne = .true.
      iBinary = 0
    case (10)
      ! TITLE
      ! NOTE: Title can be only ONE line here!!!
      What = 'S'
      if (MyGetKey(InUnit,What,iD,rD,Title1,iD,iDum,rDum) /= 0) call error()
    case (11)
      ! GAP
      What = 'R'
      if (MyGetKey(InUnit,What,iD,TheGap,Key,iD,iDum,rDum) /= 0) call error()
    case (12)
      ! END
      exit
    case (13)
      ! NODENSITY
      isDensity = .false.
    case (14)
      ! TOTAL
      isTotal = .true.
    case (15)
      ! NAME
      read(InUnit,'(a)') TheName
      !What = 'S'
      !if (MyGetKey(InUnit,What,iD,rD,TheName,iD,iDum,rDum) /= 0) call error()
      ! unfortunately MyGetKey uppercases strings!
    case (16)
      ! VB
      !isVB = 1
      call Quit_OnUserError()
    case (17)
      ! All
      isAll = .true.
    case (18)
      ! Atom
      isAtom = .true.
      iNet = -1
      iGridNpt(1) = 0
      iGridNpt(2) = 0
      iGridNpt(3) = 0
    case (19)
      ! CUBE
      iGauss = 1
      iBinary = 0
      write(u6,*) 'Cube option is moved to grid2cube'
      call Quit_OnUserError()
    case (20)
      ! Grid
      isUserGrid = .true.
      iBinary = 0
      iNet = -1
      iGridNpt(1) = 0
      iGridNpt(2) = 0
      iGridNpt(3) = 0
      What = 'I'
      if (MyGetKey(InUnit,What,nGridPoints,rD,Key,iD,iDum,rDum) /= 0) call error()
      call mma_allocate(Grid,3,nGridPoints,label='Grid')
      read(InUnit,*,iostat=iErr) Grid(:,:)
      if (iErr /= 0) call error()
    case (21)
      ! Pack
      !isMOPack = .true.
    case (22)
      ! PkLims
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,4,iDum,dTemp) /= 0) call error()
      !xLeft = dTemp(1)
      !xRight = dTemp(2)
      !xLoErr = dTemp(3)
      !xHiErr = dTemp(4)
    case (23)
      ! PkBits
      What = 'I'
      if (MyGetKey(InUnit,What,i_bits,rD,Key,iD,iDum,rDum) /= 0) call error()
      if (i_bits == 16) then
        !iyRange = 32768
        nBytesPackedVal = 2
      end if
    case (24)
      ! NoOrbitals
      NoOrb = .true.
    case (25)
      ! LINE - density on line
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,7,iDum,OneCoor) /= 0) call error()
      isTheOne = .true.
      isTotal = .true.
      iBinary = 0
      isLine = .true.
    case (26)
      ! ORANGE
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,2,iDum,Region) /= 0) call error()
      itRange = 0
      iAuMO = 1
      NoSort = .true.
    case (27)
      ! ERANGE
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,2,iDum,Region) /= 0) call error()
      itRange = 1
      iAuMO = 1
    case (28)
      ! DEBUG
      iBinary = 0
      isDebug = .true.
    case (29)
      ! CUTOFF
      What = 'R'
      if (MyGetKey(InUnit,What,iD,CutOff,Key,iD,iDum,rDum) /= 0) call error()
      isCutOff = .true.
    case (30)
      ! NOPACK
      !isMOPack = .false.
    case (31)
      ! GORI
      isCustOrig = .true.
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iDum,GridOrigin) /= 0) call error()
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iDum,GridAxis1) /= 0) call error()
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iDum,GridAxis2) /= 0) call error()
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iDum,GridAxis3) /= 0) call error()
    case (32)
      ! SELEct
      What = 'S'
      if (MyGetKey(InUnit,What,iD,rD,SelectStr,iD,iDum,rDum) /= 0) call error()
      if (nReq > 0) then
        write(u6,*) 'SELEct keyword can not be used together with ORBItals'
        call Quit_OnUserError()
      end if
      call gridExpandSelect(SelectStr)
      iAuMO = 0
    case (33)
      ! NOSOrt
      NoSort = .true.
    case (34)
      ! FILE
      read(InUnit,'(A)') FileIn
      isFileOrb = .true.
      call fileorb(FileIn,FileStr)
      write(u6,*) 'INPORB file: ',trim(FileStr)
    case (35)
      ! SPHR
      isSphere = .true.
    case (36)
      ! COLOr
      isColor = .true.
    case (37)
      ! VIRT
      isVirt = .true.
      What = 'R'
      if (MyGetKey(InUnit,What,iD,Virt,Key,iD,iDum,rDum) /= 0) call error()
    case (38)
      ! MULLiken charges per MO
      !isMULL = 1
      read(InUnit,'(A)') MULLPRT
      call upCASE(MULLPRT)
      MULLPRT = adjustl(MULLPRT)
      !if (MULLPRT(1:4) == 'LONG') isLONGPRT = 1
    case (39)
      ! SUBBLOCK
      What = 'D'
      if (MyGetKey(InUnit,What,iD,rD,Key,3,iDum,SubBlock) /= 0) call error()
      What = 'R'
      if (MyGetKey(InUnit,What,iD,rSubBlock,Key,iD,iDum,rDum) /= 0) call error()
      isSubBlock = .true.
    case (40)
      ! XDER
      iDerivative = 1
      call Quit_OnUserError()
    case (41)
      ! YDER
      iDerivative = 2
      call Quit_OnUserError()
    case (42)
      ! ZDER
      iDerivative = 3
      call Quit_OnUserError()
    case (43)
      iDerivative = 4
      call Quit_OnUserError()
    case (44)
      ! CURD (current density)
      isCurDens = .true.
      call Quit_OnUserError()
    case (45)
      ! CRXJ (current density, rxj)
      isCurDens = .true.
      !isRxJ = 1
      call Quit_OnUserError()
    case (46)
      ! UMAX (use magnetic axes)
      !iuseMaxes = 1
    case (47)
      ! NOLUSCUS
      isLuscus = .false.
      iBinary = 0
    case (48)
      ! XFIEld - ask Grid_It to compute electronic density on a DFT integration grid
      !isXField = 1
      iReadNet = iReadNet+1 !make the grid definition exclusive
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
      What = 'R'
      if (MyGetKey(InUnit,What,iD,XFminCh,Key,iD,iDum,rDum) /= 0) call error()
    case default
  end select
end do

!***********************************************************************
!                                                                      *
!                       End of input section.                          *
!                                                                      *
!***********************************************************************
!if (isLusMath == 1) return
if (isLuscus .and. (iBinary == 0)) then
  write(u6,*) 'ASCII keyword is set, but NoLUSCUS is not'
  write(u6,*) 'calling abend as the best option available'
  call Quit_OnUserError()
end if
close(InUnit)
if (isLuscus .and. isLine) then
  write(u6,*) 'LUSCUS and LINE options are not compatible'
  call Quit_OnUserError()
end if

if (iReadNet > 1) write(u6,'(a)') 'Warning: Double definition of GRID net'

! Well, there is something to do!

if (isFileOrb) then
  INPORB = FileStr
end if
call OpenGrid(INPORB)

! try to generate grid position automatically

magicValue = 0
if (iNet >= 0) then
  magicValue = GridNormal
  if (iNet == 1) magicValue = GridSparse
  if (iNet == 2) magicValue = GridDense
end if
if (isCustOrig .and. (iNet /= -1)) then
  write(u6,*) 'GORI can be used only with NPOI'
  call Quit_OnUserError()
end if
call MyCoor(isAuto,GridOrigin(1),GridOrigin(2),GridOrigin(3),GridAxis1(1),GridAxis2(2),GridAxis3(3),iGridNpt(1),iGridNpt(2), &
            iGridNpt(3),magicValue,isCustOrig)

if (isSubBlock) then
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

contains

subroutine error()
  write(u6,'(a,a,a)') 'Error during reading ',Key(1:20),'section in input file'
  call Quit_OnUserError()
end subroutine error

end subroutine Input_Grid_It
