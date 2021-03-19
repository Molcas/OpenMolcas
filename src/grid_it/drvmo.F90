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
! Copyright (C) Roland Lindh                                           *
!               Valera Veryazov                                        *
!***********************************************************************

subroutine DrvMO(iRun,INPORB)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             Valera Veryazov, Dept. Theoretical Chemistry             *
!***********************************************************************

use Symmetry_Info, only: nIrrep
use Basis_Info, only: nBas
implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "grid.fh"

!logical Debug
character str*128, Crypt*7
character INPORB*(*)
character*80 myTitle
character*128 line
real*8 pp(3)
integer nTypes(7)
logical is_error
data Crypt/'fi123sd'/
!
!---- Set size of batches
!
parameter(nIncPack=18*1024)

character cMoBlock(nIncPack*2)

nInc = nIncPack
!                                                                      *
!***********************************************************************
!                                                                      *
!... Prologue
!Debug = .false.
isEner = 1
ipCutOff = ip_iDummy

dNorm = 0
!ddNorm = 0
if (iRun == 1 .and. levelprint >= 3) then
  write(6,*)
  write(6,'(A,8I5)') 'Irreps  : ',(i,i=1,nIrrep)
  write(6,'(A,8I5)') 'Basis   : ',(nBas(i),i=0,nIrrep-1)
  write(6,'(A,3I5)') 'Grid Net: ',(iGridNpt(i),i=1,3)
  write(6,*)
end if

!... Compute the size of the densities

nCMO = 0
nMOs = 0
do iIrrep=0,nIrrep-1
  nCMO = nCMO+nBas(iIrrep)**2
  nMOs = nMOs+nBas(iIrrep)
end do

if (isUHF /= 0 .and. (isDerivative /= 0 .or. isCurDens /= 0)) then
  write(6,*) 'ERROR - Current density or derivatives not implemented for UHF!!!'
  call Abend()
end if
if (NoOrb /= 0 .and. (isDerivative /= 0 .or. isCurDens /= 0)) then
  write(6,*) 'ERROR - Current density or derivatives not implemented with the NOORB keyword!!!'
  call Abend()
end if
if (isAtom /= 0 .and. (isDerivative /= 0 .or. isCurDens /= 0)) then
  write(6,*) 'ERROR - Current density or derivatives not implemented with the ATOM keyword!!!'
  call Abend()
end if

call GetMem('CMO','ALLO','REAL',ipCMO,nCMO)

call GetMem('Ener','ALLO','REAL',ipE,nMOs)
call GetMem('Occu','ALLO','REAL',ipOcc,nMOs)
call GetMem('Occ2','ALLO','REAL',ipOoo,nMOs)
if (isVirt == 1) then
  call GetMem('ddNo','ALLO','REAL',ipdd,nMOs*nMOs)
  do i=1,nMOs*nMOs
    Work(ipdd+i-1) = Zero
  end do
end if
call GetMem('iTyp','ALLO','INTE',ipType,nMOs)
call GetMem('Vol','ALLO','REAL',ipVol,nMOs)
call GetMem('Sort','ALLO','INTE',ipSort,nMOs)
call GetMem('Nzer','ALLO','INTE',ipNZ,nMOs*2)
call GetMem('NRef','ALLO','INTE',ipGRef,nMOs)
call GetMem('DoIt','ALLO','INTE',ipDoIt,nMOs)
call GetMem('Pab','ALLO','REAL',ipPab,nMOs)
do i=0,nMOs-1
  iWork(ipGRef+i) = -1
end do

! Read information from INPORB file

LuOrb = isFreeUnit(46)

if (isUHF == 0) then

  call RdVec(INPORB,LuOrb,'COE',nIrrep,NBAS,NBAS,Work(ipCMO),Work(ipOcc),Work(ipE),iWork(ip_iDummy),myTitle,0,iErr)

  ! construct Pab
  if (NoOrb == 1) then
    !write(6,*) 'nCMO,nMOs', nCMO,nMOs
    call makePab(Work(ipCMO),Work(ipOcc),Work(ipPab),nMOs,nMOs,nIrrep,nBas)
    !write(6,*) 'Pab=', (Work(ipPab+i),i=0,nMOs-1)
  end if
  if (iErr == 1) then
    do j=0,nMOs-1
      Work(ipE+j) = Zero
    end do
  end if
  call RdVec(INPORB,LuOrb,'I',nIrrep,NBAS,NBAS,Work(ip_Dummy),Work(ip_Dummy),Work(ip_Dummy),iWork(ipType),myTitle,0,iErr)
  if (iErr == 1) then
    do j=0,nMOs-1
      iWork(ipType+j) = 0
    end do
  end if
else  ! UHF case
  ! allocate memory for extra arrays.
  call GetMem('CMO_ab','ALLO','REAL',ipCMO_ab,nCMO)
  call GetMem('Ener_ab','ALLO','REAL',ipE_ab,nMOs)
  call GetMem('Occu_ab','ALLO','REAL',ipOcc_ab,nMOs)
  call GetMem('Sort_ab','ALLO','INTE',ipSort_ab,nMOs)
  call GetMem('NRef_ab','ALLO','INTE',ipGRef_ab,nMOs)
  call GetMem('DoIt_ab','ALLO','INTE',ipDoIt_ab,nMOs)

  call RdVec_(INPORB,LuOrb,'COE',1,nIrrep,NBAS,NBAS,Work(ipCMO),Work(ipCMO_ab),Work(ipOcc),Work(ipOcc_ab),Work(ipE),Work(ipE_ab), &
              iWork(ip_iDummy),myTitle,0,iErr,iWFtype)
  ! it can be only after SCF, so we do not need TypeIndex info
  do j=0,nMOs-1
    iWork(ipType+j) = 0
  end do

end if
do j=1,7
  nTypes(j) = 0
end do
do j=0,nMOs-1
  jj = iWork(ipType+j)
  if (jj > 0) nTypes(jj) = nTypes(jj)+1
end do

close(LuOrb)

! Calculate net.

if (iRun == 1) then
  write(6,'(A)') '   Input vectors read from INPORB'
  write(6,'(A,A)') '   Orbital file label: ',myTitle(:mylen(myTitle))
end if
do j=1,3
  iGridNpt(j) = iGridNpt(j)+1
end do

nCoor = iGridNpt(1)*iGridNpt(2)*iGridNpt(3)
iiCoord = nCoor

!---- if isCutOff is in used - recalculate nCoor

if (isCutOff == 1) then
  if (isUserGrid == 1) then
    write(6,*) 'Not implemented'
    call Quit_OnUserError()
  end if

  call GetMem('CUTFL','ALLO','INTE',ipCutOff,nCoor)
  ie1 = max(iGridNpt(1)-1,1)
  ie2 = max(iGridNpt(2)-1,1)
  ie3 = max(iGridNpt(3)-1,1)
  !write(6,*) 'vv', ie1, ie2, ie3
  iiCoord = 0
  iiiCoord = 0
  do i1=0,ie1
    do i2=0,ie2
      do i3=0,ie3
        pp(1) = GridOrigin(1)+GridAxis1(1)*i1/ie1+GridAxis2(1)*i2/ie2+GridAxis3(1)*i3/ie3
        pp(2) = GridOrigin(2)+GridAxis1(2)*i1/ie1+GridAxis2(2)*i2/ie2+GridAxis3(2)*i3/ie3
        pp(3) = GridOrigin(3)+GridAxis1(3)*i1/ie1+GridAxis2(3)*i2/ie2+GridAxis3(3)*i3/ie3
        !write(6,'(3f8.4)') pp
        ishow = 0
        do ii=1,nAtoms
          iaia = 0
          if (abs(Work(ipCoor+ii*3-3)-pp(1)) < CutOff) iaia = iaia+1
          if (abs(Work(ipCoor+ii*3-3+1)-pp(2)) < CutOff) iaia = iaia+1
          if (abs(Work(ipCoor+ii*3-3+2)-pp(3)) < CutOff) iaia = iaia+1
          if (iaia == 3) ishow = 1
        end do
        if (ishow == 1) then
          iiCoord = iiCoord+1
          iWork(ipCutOff+iiiCoord) = 1
        else
          iWork(ipCutOff+iiiCoord) = 0
        end if
        iiiCoord = iiiCoord+1
      end do
    end do
  end do
  write(6,*) nCoor-iiCoord,' points are eliminated'
  !write(6,*) 'old=',nCoor,' New=', iiCoord
  !nCoor = iiCoord
end if
if (isTheOne == 1) nCoor = int(OneCoor(7)+0.3)
if (isAtom == 1) nCoor = nAtoms
if (isUserGrid == 1) nCoor = nGridPoints
write(6,*) ' Number of grid points in file:  ',nCoor
!call iXML('nPoints',nCoor)

!***********************************************************************
! And now we had to choose orbitals to draw.
!
! Sometime we had to make an automatic guess....
!***********************************************************************

call PickOrb(ipNz,ipSort,ipGref,ipSort_ab,ipGref_ab,ipVol,ipE,ipOcc,ipE_ab,ipOcc_ab,nShowMOs,nShowMOs_ab,isener,nMOs,myTitle,ipType)

!---- Start run over sets of grid points

!if (isBinary == 1) then
nShowMOs = nShowMOs+isDensity+isSphere+isColor
!else
!  nShowMOs=1
!end if
if (isUHF == 1) nShowMOs_ab = nShowMOs_ab+isDensity+isSphere+isColor
nShowMOs2 = nShowMOs+nShowMOs_ab
write(6,*)
write(6,*) ' Total number of MOs               :',nMOs
!call iXML('nMOs',nMOs)
write(6,*) ' Number MOs for grid               :',nShowMOs2
write(6,*) ' Batches processed in increments of:',nInc
write(6,*)
iPrintCount = 0

call GetMem('MOValue','ALLO','REAL',ipMO,nInc*nMOs)
call GetMem('DOValue','ALLO','REAL',ipOut,nInc)

!if (imoPack .ne. 0) then
!  call GetMem('PackedBlock','ALLO','INTE',ipPBlock,nInc)
!else
call GetMem('PackedBlock','ALLO','INTE',ipPBlock,1)
iWork(ipPBlock) = 0
!end if

!... Allocate memory for the some grid points

call GetMem('Coor','ALLO','REAL',ipC,nInc*3)
iSphrDist = ip_Dummy
iSphrColor = ip_Dummy

if (isSphere == 1) then
  call GetMem('SpDi','ALLO','REAL',iSphrDist,nInc)
end if
if (isColor == 1) then
  call GetMem('SpCo','ALLO','REAL',iSphrColor,nInc)
end if

! check grids to calculate

do i=0,nMOs-1
  iWork(ipDoIt+i) = 0
  if (isUHF == 1) iWork(ipDoIt_ab+i) = 0
end do

if (NoOrb == 1) goto 550
do i=1,nShowMOs-isDensity-isSphere-isColor
  iWork(ipDoIt+iWork(ipGRef+i-1)-1) = 1
end do
if (isUHF == 1) then
  do i=1,nShowMOs_ab-isDensity-isSphere-isColor
    iWork(ipDoIt_ab+iWork(ipGRef_ab+i-1)-1) = 1
  end do
end if
ifpartial = 1-isTotal
if (isTotal == 1) then
  do i=0,nMOs-1
    if (abs(Work(ipOcc+i)) > Zero) then
      iWork(ipDoIt+i) = 1
    end if
    if (isUHF == 1) then
      if (abs(Work(ipOcc_ab+i)) > Zero) then
        iWork(ipDoIt_ab+i) = 1
      end if
    end if
  end do
else
  do i=0,nMOs-1
    if (Work(ipOcc+i) > Zero .and. iWork(ipDoIt+i) /= 1) then
      ifpartial = 1
    end if
    if (isUHF == 1) then
      if (Work(ipOcc_ab+i) > Zero .and. iWork(ipDoIt_ab+i) /= 1) then
        ifpartial = 1
      end if
    end if
  end do
end if
550 continue

! If plotting VB orbitals : count active orbitals and make sure
! all are included in DOIT :
!
! nothing to do with UHF

iCRSIZE = 1
NBYTES = 10
NINLINE = 10
!write(6,*) 'prepare  header '

call PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks)

!write(6,*) 'HERE header isdone'

LuVal_ = LuVal
if (isLuscus == 1) LuVal_ = LID
call PrintTitles(LuVal_,nShowMOs,isDensity,nMOs,iWork(ipGRef),isEner,Work(ipOcc),iWork(ipType),Crypt,iWork(ipNZ),Work(ipE),VBocc, &
                 ifpartial,isLine,isSphere,isColor,ISLUSCUS,ncoor,nBlocks,nInc)
if (isUHF == 1) then
  LuVal_ab_ = LuVal_ab
  if (isLuscus == 1) LuVal_ab_ = LID_ab
  call PrintTitles(LuVal_ab_,nShowMOs_ab,isDensity,nMOs,iWork(ipGRef_ab),isEner,Work(ipOcc_ab),iWork(ipType),Crypt,iWork(ipNZ), &
                   Work(ipE_ab),VBocc,ifpartial,isLine,isSphere,isColor,ISLUSCUS,ncoor,nBlocks,nInc)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Loop over grid points in batches of nInc

ive3 = max(iGridNpt(3)-1,1)
ive2 = max(iGridNpt(2)-1,1)
ive1 = max(iGridNpt(1)-1,1)
iv3 = 0
iv2 = 0
iv1 = 0
!if (isAtom== 1) goto 6000

iiiCoord = 0
!if (isCutOff == 1) nCoor = iiCoord
!ccccccccccccc  main loop starts here  ccccccccccccccccccccccccccccccccc
iShiftCut = 0
do iSec=1,nCoor,nInc
  mCoor = min(nInc,nCoor-iSec+1)
  !write(status,'(a,i8,a,i8)') ' batch ',iSec,' out of ',nCoor/nInc
  !call StatusLine('grid_it: ',status)
  ! Generate next portion of points..
  if (isTheOne == 0) then
    ! general case: we have a CUBIC box.
    ipPO = 0
    667 continue
    iiiCoord = iiiCoord+1
    gv3 = 1.0d+0*iv3/ive3
    gv2 = 1.0d+0*iv2/ive2
    gv1 = 1.0d+0*iv1/ive1

    if (isUserGrid == 0) then
      if (isCutOff == 0) then
        Work(ipC+ipPO*3) = GridOrigin(1)+GridAxis1(1)*gv1+GridAxis2(1)*gv2+GridAxis3(1)*gv3
        Work(ipC+ipPO*3+1) = GridOrigin(2)+GridAxis1(2)*gv1+GridAxis2(2)*gv2+GridAxis3(2)*gv3
        Work(ipC+ipPO*3+2) = GridOrigin(3)+GridAxis1(3)*gv1+GridAxis2(3)*gv2+GridAxis3(3)*gv3
      else
        ! using ipCutOff
        Work(ipC+ipPO*3) = 40
        Work(ipC+ipPO*3+1) = 40
        Work(ipC+ipPO*3+2) = 40
        if (iWork(ipCutOff+iiiCoord-1) == 1) then
          Work(ipC+ipPO*3) = GridOrigin(1)+GridAxis1(1)*gv1+GridAxis2(1)*gv2+GridAxis3(1)*gv3
          Work(ipC+ipPO*3+1) = GridOrigin(2)+GridAxis1(2)*gv1+GridAxis2(2)*gv2+GridAxis3(2)*gv3
          Work(ipC+ipPO*3+2) = GridOrigin(3)+GridAxis1(3)*gv1+GridAxis2(3)*gv2+GridAxis3(3)*gv3
        end if
      end if
    else
      Work(ipC+ipPO*3+1-1) = Work(ipGrid+(iSec+ipPO-1)*3+1-1)
      Work(ipC+ipPO*3+2-1) = Work(ipGrid+(iSec+ipPO-1)*3+2-1)
      Work(ipC+ipPO*3+3-1) = Work(ipGrid+(iSec+ipPO-1)*3+3-1)
    end if
    ! make a local copy of the weights of the corresponding grid points:
    iv3 = iv3+1
    if (iv3 > ive3) then
      iv3 = 0
      iv2 = iv2+1
    end if
    if (iv2 > ive2) then
      iv2 = 0
      iv1 = iv1+1
    end if
    ipPO = ipPO+1
    !write(6,*) 'ipo',ipP0
    if (ipPO <= mCoor-1) goto 667
    ! end of CUBIC box
  else
    ! coords are given in specific points
    ! coords for DEBUG mode
    if (isLine == 0) then
      do i=0,nCoor-1
        Work(ipC+i*3) = OneCoor(1)+OneCoor(4)*i
        Work(ipC+1+i*3) = OneCoor(2)+OneCoor(5)*i
        Work(ipC+2+i*3) = OneCoor(3)+OneCoor(6)*i
      end do
    else
      ! LINE keyword
      do i=0,nCoor-1
        Work(ipC+i*3) = OneCoor(1)+(OneCoor(4)-OneCoor(1))*i/(nCoor-1)
        Work(ipC+1+i*3) = OneCoor(2)+(OneCoor(5)-OneCoor(2))*i/(nCoor-1)
        Work(ipC+2+i*3) = OneCoor(3)+(OneCoor(6)-OneCoor(3))*i/(nCoor-1)
      end do
    end if
  end if
  ! end of coordinates.
  !VV: FIXME; separate color and sphere
  !if (isSphere == 1 .and. isColor == 1) then
  !  call Sphr_Grid(Work(ipCoor),mCoor,Work(ipC),Work(iSphrDist),Work(iSphrColor))
  !end if

  if (NoOrb == 1) then
    nDrv = 0
    call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipPab),nMOs,iWork(ipDoIt),nDrv,1)
  else
    nDrv = 0
    call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO),nCMO,iWork(ipDoIt),nDrv,1)
  end if

  !... Write out values

  nLine = min(nShowMOs,20)
  nSLine = nLine*mCoor
  if (isLine == 0) then
    nSLine = 1
    nLine = 1
  end if
  call GetMem('Line','ALLO','REAL',ipLine,nSLine)
  if (levelprint < 2) iPrintCount = 100
  !VV BUG Update ipcutOFF

  call DumpM2Msi(iRun,Luval_,LID,nShowMOs,isDensity,nMOs,iWork(ipGRef),Work(ipOcc),Work(ipMO),Work(ipOut),mCoor,iGauss,nInc, &
                 imoPack,iWork(ipPBlock),cMoBlock,nBytesPackedVal,dnorm,Crypt,VbOcc,isTheOne,isLine,isBinary,isEner,iWork(ipType), &
                 iWork(ipNZ),Work(ipE),Work(ipLine),nLine,Work(ipC),iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut), &
                 isSphere,Work(iSphrDist),isColor,Work(iSphrColor),isLuscus,NBYTES,NINLINE)

    !if (isXField == 1) call dcopy_(mCoor,Work(ipOut),1,Work(ipOutXF),1)
  if (isUHF == 1) then
    !VV:
    nDrv = 0
    call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO_ab),nCMO,iWork(ipDoIt_ab),nDrv,1)

    !... Write out values

    call DumpM2Msi(iRun,Luval_ab_,LID_ab,nShowMOs_ab,isDensity,nMOs,iWork(ipGRef_ab),Work(ipOcc_ab),Work(ipMO),Work(ipOut),mCoor, &
                   iGauss,nInc,imoPack,iWork(ipPBlock),cMoBlock,nBytesPackedVal,dnorm,Crypt,VbOcc,isTheOne,isLine,isBinary,isEner, &
                   iWork(ipType),iWork(ipNZ),Work(ipE_ab),Work(ipLine),nLine,Work(ipC),iPrintCount,isDebug,isCutOff, &
                   iWork(ipCutOff+iShiftCut),isSphere,Work(iSphrDist),isColor,Work(iSphrColor),isLuscus,NBYTES,NINLINE)

  end if ! end if (isUHF == 1)

  call GetMem('Line','FREE','REAL',ipLine,nSLine)
  iShiftCut = iShiftCut+mCoor
  if (isVirt == 1) then
    do iiMO=1,nMOs
      do jjMO=1,nMOs
        if (Work(ipOcc+iiMO-1) > 1.1 .and. Work(ipOcc+jjMO-1) < 0.9) then
          !  here if this is a pair Occ-Virt

          do i=1,nMOs
            Work(ipOoo+i-1) = 0
            if (i == iiMO) Work(ipOoo+i-1) = Two-Virt
            if (i == jjMO) Work(ipOoo+i-1) = Virt
          end do

          call outmo(0,2,Work(ipMO),Work(ipOoo),Work(ipOut),mCoor,nMOs)
          dd = 0
          do j=1,mCoor
            dd = dd+Work(ipOut+j-1)
          end do
          !ddNorm = ddNorm+dd
          call save_ddNorm(dd,iiMO,jjMO,Work(ipdd),nMOs)
        end if
      end do
    end do
  end if
end do !iSec
!ccccccccccccc  main loop ends here  ccccccccccccccccccccccccccccccccccc

write(6,*)
! Check norms

if (isAtom /= 1) then

  det3 = GridAxis1(1)*(GridAxis2(2)*GridAxis3(3)-GridAxis2(3)*GridAxis3(2))- &
         GridAxis1(2)*(GridAxis2(1)*GridAxis3(3)-GridAxis2(3)*GridAxis3(1))+ &
         GridAxis1(3)*(GridAxis2(1)*GridAxis3(2)-GridAxis2(2)*GridAxis3(1))
  det3 = det3/((iGridNpt(1)-1)*(iGridNpt(2)-1)*(iGridNpt(3)-1))
  dNorm = dNorm*det3

  !write(6,*) 'dNorm=',dNorm

  if (isVirt == 1) then
    call print_ddNorm(nMOs,Work(ipdd),det3)

  !write(6,*) 'ddNorm=',ddNorm*det3

  end if

  call Add_Info('GRIDIT_NORM',[dNorm],1,6)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!call DAClos(LuVal)
if (ISLUSCUS == 1) then
  LINE = ' '
  call PRINTLINE(LUVAL_,LINE,1,0)
  LINE = ' </DENSITY>'
  call PRINTLINE(LUVAL_,LINE,11,0)
  LINE = ' <BASIS>'
  call PRINTLINE(LUVAL_,LINE,8,0)
  LINE = ' </BASIS>'
  call PRINTLINE(LUVAL_,LINE,9,0)
  LINE = ' <INPORB>'
  call PRINTLINE(LUVAL_,LINE,9,0)
end if

if (isLine == 0) then
  write(str,'(9i8)') nIrrep,(nBas(i),i=0,nIrrep-1)
  call PrintLine(LuVal_,STR,72,0)
  !if (isBinary == 0) then
  !  write(LuVal_,'(a)') str
  !else
  !  write(LuVal_) str
  !end if
  ! Well, to avoid rewritting of Cerius2 we use old INPORB format temporary!
  LuOrb = isFreeUnit(46)
  call molcas_open_ext2(luorb,INPORB,'sequential','formatted',iostat,.false.,irecl,'old',is_error)
  !4001 read(LuOrb,'(a)',err=5000,end=5000) str
  !if (str /= '#ORB') goto 4001
  5001 read(LuOrb,'(a)',err=5000,end=5000) str
  !if (str(1:1) == '#') goto 5001
  iLen = len(str)
  do i=iLen,1,-1
    if (str(i:i) /= ' ') goto 5060
  end do
  5060 continue
  call PrintLine(LuVal_,STR,I,0)
  !5060 if (isBinary == 0) then
  !  write(LuVal_,'(a)') str(1:i)
  !else
  !  write(LuVal_) str(1:i)
  !end if
  goto 5001
  5000 close(unit=LuOrb)
end if ! isLine
close(unit=LuVal)
if (isUHF == 1) close(unit=LuVal_ab)
if (isTheOne == 1) then
  MM = mCoor-1
  if (MM > 10) MM = 10
  call Add_Info('GRIDIT_ONE',Work(ipOut),MM,6)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!... Epilogue, end

!6000  continue

if (isAtom == 1) then
  mCoor = nCoor
  nDrv = 0
  call MOEval(Work(ipMO),nMOs,mCoor,Work(ipCoor),Work(ipCMO),nCMO,iWork(ipDoIt),nDrv,1)
  call outmo(0,2,Work(ipMO),Work(ipOcc),Work(ipOut),nCoor,nMOs)
  write(6,'(60a1)') ('*',i=1,60)
  if (ifpartial == 0) then
    write(6,'(a5,3a10,a20)') 'Atom','x','y','z','Density'
  else
    write(6,'(a5,3a10,a20)') 'Atom','x','y','z','Density (partial)'
  end if
  do i=0,nAtoms-1
    write(6,'(a5,3f10.3,e20.10)') AtomLbl(i+1),Work(ipCoor+i*3),Work(ipCoor+i*3+1),Work(ipCoor+i*3+2),Work(ipOut+i)
  end do
  call Add_Info('GRIDIT_ATOM',Work(ipOut),nAtoms,6)

  call GetMem('Coor','FREE','REAL',ipCoor,3*nAtoms)

end if

if (ISLUSCUS == 1) then
  call PRTLUSENDGRID(LID)
  if (isUHF == 1) call PRTLUSENDGRID(LID_ab)
end if

if (isSphere == 1) then
  call GetMem('SpDi','FREE','REAL',iSphrDist,nInc)
end if
if (isColor == 1) then
  call GetMem('SpCo','FREE','REAL',iSphrColor,nInc)
end if
call GetMem('Coor','FREE','REAL',ipC,nInc*3)

call GetMem('DOValue','FREE','REAL',ipOut,nInc)
call GetMem('MOValue','FREE','REAL',ipMO,nInc*nMOs)

if (isUHF == 1) then
  call GetMem('CMO_ab','FREE','REAL',ipCMO_ab,nCMO)
  call GetMem('Ener_ab','FREE','REAL',ipE_ab,nMOs)
  call GetMem('Occu_ab','FREE','REAL',ipOcc_ab,nMOs)
  call GetMem('Sort_ab','FREE','INTE',ipSort_ab,nMOs)
  call GetMem('NRef_ab','FREE','INTE',ipGRef_ab,nMOs)
  call GetMem('DoIt_ab','FREE','INTE',ipDoIt_ab,nMOs)
end if
if (isUserGrid == 1) call GetMem('Grid','FREE','REAL',ipGrid,nGridPoints*3)
!if (imoPack /= 0) then
!  call GetMem('PackedBlock','FREE','INTE',ipPBlock,nInc)
!else
call GetMem('PackedBlock','FREE','INTE',ipPBlock,1)
!end if

call GetMem('Pab','FREE','REAL',ipPab,nMOs)
call GetMem('DoIt','FREE','INTE',ipDoIt,nMOs)

call GetMem('NRef','FREE','INTE',ipGRef,nMOs)
call GetMem('Nzer','FREE','INTE',ipNZ,nMOs*2)
call GetMem('Sort','FREE','INTE',ipSort,nMOs)
call GetMem('Vol','FREE','REAL',ipVol,nMOs)
call GetMem('iTyp','FREE','INTE',ipType,nMOs)
if (isVirt == 1) then
  call GetMem('ddNo','FREE','REAL',ipdd,nMOs*nMOs)
end if

call GetMem('Occ2','FREE','REAL',ipOoo,nMOs)
call GetMem('Occu','FREE','REAL',ipOcc,nMOs)
call GetMem('Ener','FREE','REAL',ipE,nMOs)

!if (isWDW == 1) call GetMem('WDW','FREE','REAL',ipWdW,nCenter)
call GetMem('CMO','FREE','REAL',ipCMO,nCMO)

if (isAtom == 0 .and. isXField == 0) call GetMem('Coor','FREE','REAL',ipCoor,3*nAtoms)

if (isCutOff == 1) call GetMem('CUTFL','FREE','INTE',ipCutOff,nCoor)

return

end subroutine DrvMO
