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

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"

character Key*80
character INPORB*(*)
integer iTemp(3)
dimension dTemp(4)
character AllKeys*265
character SelectStr*120
character FileStr*256, FileIn*256
character MULLprt*80
!
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
TheGap = 4.0d0
nMOmin = 0
nGrid = 5
nReq = -1
isDensity = 1
isDerivative = 0
isCurDens = 0
isRxJ = 0
iuseMaxes = 0
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
CutOff = 2.5d0
iCustOrig = 0
NoSort = 0
isFileOrb = 0
isSphere = 0
isColor = 0
isVirt = 0
isMULL = 0
isLONGPRT = 0
isWDW = 0
iSubBlock = 0
isLuscus = 1
!isLusMath = 0
!aLusMath = -1
! not preparing the GRIDCHARGE file as external source for XFIELD input
isXField = 0
XFminCh = 0.d0
! Default values for packing
imoPack = 0
isBinPack = 0
xLeft = 0.005d0
xRight = 0.7d0
! (really, half range:)
iyRange = 128
nBytesPackedVal = 1
iMinYLeft = 4
xLoErr = 0.10d0
xHiErr = 0.25d0
INPORB = 'INPORB'
if (iRun == 0) then
  ! make defaults for a fast run via call from other module
  isNet = 1  ! set sparse
  !isBinary = 0 ! temporary set Ascii output
  !iMaxUp = 1
  !iMaxDown = 5
  !imoPack = 0 ! packed grids not working with 64bit (?)
  !isCutOff = 1
  !goto 500
end if

! KeyWord directed input

InUnit = 5
! Function MyGetKey (InUnit, What, IValue, RValue, SValue, N, IArray, RArray)

call RdNLst(InUnit,'GRID_IT')
998 if (MyGetKey(InUnit,'S',iD,rD,Key,iD,[iD],[rD]) /= 0) goto 997
iKey = index(AllKeys,Key(1:4))
if (iKey == 0 .or. (iKey-1)/5*5 /= (iKey-1)) then
  write(6,'(a,a)') 'Unrecognized keyword in input file:',Key(1:4)
  call Quit_OnUserError()
end if
iKey = (iKey-1)/5+1

if (iKey == 1) then
  ! PRIN
  if (MyGetKey(InUnit,'I',n,rD,Key,iD,[iD],[rD]) /= 0) goto 666
  do j=1,n
    if (MyGetKey(InUnit,'A',iD,rD,Key,2,iTemp,[rD]) /= 0) goto 666
    !write(6,*) 'debug'
    !nPrint(iTemp(1)) = iTemp(2)
  end do
end if
if (iKey == 2) then
  ! BINARY = default
  isBinary = 1
end if
if (iKey == 3) then
  ! ASCII = for debug
  !write(6,*) ' Keyword ASCII is obsolete'
  !write(6,*) ' It can be used only for debugging purpose'
  !write(6,*) ' Note that .lus files produced with this option '
  !write(6,*) '      can not be visualised'
  isBinary = 0
end if
if (iKey == 4) then
  ! NPOI
  if (MyGetKey(InUnit,'A',iD,rD,Key,3,iGridNpt,[rD]) /= 0) goto 666
  isNet = -1
  isReadNet = isReadNet+1
end if
if (iKey == 5) then
  ! DENSE - dense grid network..
  isNet = 2
  isReadNet = isReadNet+1
end if
if (iKey == 6) then
  ! SPARSE - rare grid network..
  isNet = 1
  isReadNet = isReadNet+1
end if
if (iKey == 7) then
  ! ORBI Orbitals
  if (nReq > 0) then
    write(6,*) 'ORBI keyword can not be used together with SELEct'
    call Quit_OnUserError()
  end if
  if (MyGetKey(InUnit,'I',nReq,rD,Key,iD,[iD],[rD]) /= 0) goto 666

  if (nReq > MAXGRID) then
    write(6,'(a,i5,a,i5)') 'Too many requested orbitals ',nReq,'>',MAXGRID
    call Quit_OnUserError()
  end if
  read(inUnit,*,err=666,end=666) (iReq(i),i=1,nReq*2)
  !if (MyGetKey(InUnit,'A',iD,rD,Key,nReq*2,iReq,[rD]) /= 0) goto 666
  isAuMO = 0
end if
if (iKey == 8) then
  ! REGION
  if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) goto 666
  itRange = 1
  isAuMO = 1
  write(6,*) ' *** Warning keyword REGION is obsolete'
  write(6,*) ' ***         assumimg Energy range '
end if
if (iKey == 9) then
  ! ONE - debug option
  if (MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor) /= 0) goto 666
  isTheOne = 1
  isBinary = 0
end if
if (iKey == 10) then
  ! TITLE
  ! NOTE: Title can be only ONE line here!!!
  if (MyGetKey(InUnit,'S',iD,rD,Title1,iD,[iD],[rD]) /= 0) goto 666
end if
if (iKey == 11) then
  ! GAP
  if (MyGetKey(InUnit,'R',iD,TheGap,Key,iD,[iD],[rD]) /= 0) goto 666
end if

if (iKey == 12) then
  ! END
  goto 997
end if
if (iKey == 13) then
  ! NODENSITY
  isDensity = 0
end if
if (iKey == 14) then
  ! TOTAL
  isTotal = 1
end if
if (iKey == 15) then
  ! NAME
  read(InUnit,'(a)') TheName
  !if (MyGetKey(InUnit,'S',iD,rD,TheName,iD,[iD],[rD]) /= 0) goto 666
  ! unfortunately MyGetKey uppercases strings!
end if
if (iKey == 16) then
  ! VB
  !isVB = 1
  call Quit_OnUserError()
end if
if (iKey == 17) then
  ! All
  isAll = 1
end if

if (iKey == 18) then
  ! Atom
  isAtom = 1
  isNet = -1
  iGridNpt(1) = 0
  iGridNpt(2) = 0
  iGridNpt(3) = 0
end if
if (iKey == 19) then
  ! CUBE
  iGauss = 1
  isBinary = 0
  write(6,*) 'Cube option is moved to grid2cube'
  call Quit_OnUserError()
end if
if (iKey == 20) then
  ! Grid
  isUserGrid = 1
  isBinary = 0
  isNet = -1
  iGridNpt(1) = 0
  iGridNpt(2) = 0
  iGridNpt(3) = 0
  if (MyGetKey(InUnit,'I',nGridPoints,rD,Key,iD,[iD],[rD]) /= 0) goto 666
  call GetMem('Grid','ALLO','REAL',ipGrid,nGridPoints*3)
  read(InUnit,*,Err=666,end=666) (Work(ipGrid+i-1),i=1,nGridPoints*3)
end if
if (iKey == 21) then
  ! Pack
  !imoPack = 1
end if
if (iKey == 22) then
  ! PkLims
  if (MyGetKey(InUnit,'D',iD,rD,Key,4,[iD],dTemp) /= 0) goto 666
  xLeft = dTemp(1)
  xRight = dTemp(2)
  xLoErr = dTemp(3)
  xHiErr = dTemp(4)
end if
if (iKey == 23) then
  ! PkBits
  if (MyGetKey(InUnit,'I',ibits,rD,Key,iD,[iD],[rD]) /= 0) goto 666
  if (ibits == 16) then
    iyRange = 32768
    nBytesPackedVal = 2
  end if
end if
if (iKey == 24) then
  ! NoOrbitals
  NoOrb = 1
end if
if (iKey == 25) then
  ! LINE - density on line
  if (MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor) /= 0) goto 666
  isTheOne = 1
  isTotal = 1
  isBinary = 0
  isLine = 1
end if
if (iKey == 26) then
  ! ORANGE
  if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) goto 666
  itRange = 0
  isAuMO = 1
  NoSort = 1
end if
if (iKey == 27) then
  ! ERANGE
  if (MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region) /= 0) goto 666
  itRange = 1
  isAuMO = 1
end if
if (iKey == 28) then
  ! DEBUG
  isBinary = 0
  isDebug = 1
end if
if (iKey == 29) then
  ! CUTOFF
  if (MyGetKey(InUnit,'R',iD,CutOff,Key,iD,[iD],[rD]) /= 0) goto 666
  isCutOff = 1
end if
if (iKey == 30) then
  ! NOPACK
  !imoPack = 0
end if
if (iKey == 31) then
  ! GORI
  iCustOrig = 1
  if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridOrigin) /= 0) goto 666
  if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis1) /= 0) goto 666
  if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis2) /= 0) goto 666
  if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis3) /= 0) goto 666
end if
if (iKey == 32) then
  ! SELEct
  if (MyGetKey(InUnit,'S',iD,rD,SelectStr,iD,[iD],[rD]) /= 0) then
    goto 666
  end if
  if (nReq > 0) then
    write(6,*) 'SELEct keyword can not be used together with ORBItals'
    call Quit_OnUserError()
  end if
  call gridExpandSelect(SelectStr)
  isAuMO = 0
end if
if (iKey == 33) then
  ! NOSOrt
  NoSort = 1
end if
if (iKey == 34) then
  ! FILE
  read(InUnit,'(A)') FileIn
  isFileOrb = 1
  call fileorb(FileIn,FileStr)
  write(6,*) 'INPORB file: ',FileStr(:mylen(FileStr))
end if
if (iKey == 35) then
  ! SPHR
  isSphere = 1
end if
if (iKey == 36) then
  ! COLOr
  isColor = 1
end if
if (iKey == 37) then
  ! VIRT
  isVirt = 1
  if (MyGetKey(InUnit,'R',iD,Virt,Key,iD,[iD],[rD]) /= 0) goto 666
end if
if (iKey == 38) then
  ! MULLiken charges per MO
  isMULL = 1
  read(InUnit,'(A)') MULLPRT
  call upCASE(MULLPRT)
  call LeftAd(MULLPRT)
  if (MULLPRT(1:4) == 'LONG') isLONGPRT = 1
end if
if (iKey == 39) then
  ! SUBBLOCK
  if (MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],SubBlock) /= 0) goto 666
  if (MyGetKey(InUnit,'R',iD,rSubBlock,Key,iD,[iD],[rD]) /= 0) goto 666
  iSubBlock = 1
end if
! XDER,YDER,ZDER,GDER
if (iKey == 40) then
  isDerivative = 1
  call Quit_OnUserError()
end if
if (iKey == 41) then
  isDerivative = 2
  call Quit_OnUserError()
end if
if (iKey == 42) then
  isDerivative = 3
  call Quit_OnUserError()
end if
if (iKey == 43) then
  isDerivative = 4
  call Quit_OnUserError()
end if

if (iKey == 44) then
  ! CURD (current density)
  isCurDens = 1
  call Quit_OnUserError()
end if

if (iKey == 45) then
  ! CRXJ (current density, rxj)
  isCurDens = 1
  isRxJ = 1
  call Quit_OnUserError()
end if
if (iKey == 46) then
  ! UMAX (use magnetic axes)
  iuseMaxes = 1
end if
if (iKey == 47) then
  ! NOLUSCUS
  isLuscus = 0
  isBinary = 0
end if
if (iKey == 48) then
  ! XFIEld - ask Grid_It to compute electronic density on a DFT integration grid
  isXField = 1
  isReadNet = isReadNet+1 !make the grid definition exclusive
end if
if (iKey == 49) then
  ! LUS1
  write(6,*) 'Not implemented'
  !isLusMath = 1
  !read(InUnit,'(a)') LUS1
end if
if (iKey == 50) then
  ! LUS2
  write(6,*) 'Not implemented'
  !isLusMath = 1
  !read(InUnit,'(a)') LUS2
end if
if (iKey == 51) then
  ! PLUS
  write(6,*) 'Not implemented'
  !aLusMath = 1
end if
if (iKey == 52) then
  ! MINUS
  write(6,*) 'Not implemented'
  !aLusMath = -1
end if
if (iKey == 53) then
  ! XFMI xfield minimum charge of each grid point to be stored
  if (MyGetKey(InUnit,'R',iD,XFminCh,Key,iD,[iD],[rD]) /= 0) goto 666
end if
goto 998

666 write(6,'(a,a,a)') 'Error during reading ',Key(1:20),'section in input file'
call Quit_OnUserError()

!***********************************************************************
!                                                                      *
!                       End of input section.                          *
!                                                                      *
!***********************************************************************
997 continue
!if (isLusMath == 1) return
if (isLuscus == 1 .and. isBinary == 0) then
  write(6,*) 'ASCII keyword is set, but NoLUSCUS is not'
  write(6,*) 'calling abend as the best option available'
  call Quit_OnUserError()
end if
close(InUnit)
if (isLuscus == 1) then
  if (isLine /= 0) then
    write(6,*) 'LUSCUS and LINE options are not compatible'
    call Quit_OnUserError()
  end if
end if

if (isReadNet > 1) write(6,'(a)') 'Warning: Double definition of GRID net'

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
    write(6,*) 'GORI can be used only with NPOI'
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

end subroutine Input_Grid_It
