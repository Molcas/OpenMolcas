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
use grid_it_globals, only: AtomLbl, Coor, CutOff, Grid, GridAxis1, GridAxis2, GridAxis3, GridOrigin, iGauss, iGridNpt, isAtom, &
                           iBinary, isColor, isCurDens, isCutOff, isDebug, isDensity, iDerivative, isLine, isLuscus, isMOPack, &
                           isSphere, isTheOne, isTotal, isUHF, isUserGrid, isVirt, levelprint, LID, LID_ab, LuVal, LuVal_ab, &
                           nAtoms, nBytesPackedVal, nGridPoints, NoOrb, OneCoor, Virt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iRun
character(len=*), intent(in) :: INPORB
integer(kind=iwp) :: i, i1, i2, i3, iaia, iCRSIZE, idum(1), ie1, ie2, ie3, iErr, ii, iiCoord, iiiCoord, iiMO, iIrrep, iLen, &
                     istatus, ipPO, iPrintCount, irecl, iSec, iShiftCut, ishow, iv1, iv2, iv3, ive1, ive2, ive3, iWFtype, j, jj, &
                     jjMO, LuOrb, LuVal_, LuVal_ab_, mCoor, MM, nBlocks, NBYTES, nCMO, nCoor, nDrv, nInc, NINLINE, nLine, nMOs, &
                     nShowMOs, nShowMOs2, nShowMOs_ab, nTypes(7)
real(kind=wp) :: dd, det3, dNorm, dum(1), gv1, gv2, gv3, pp(3)
logical(kind=iwp) :: ifpartial, is_error, isEner
character(len=128) :: line, str
character(len=80) :: myTitle
integer(kind=iwp), allocatable :: DoIt(:), DoIt_ab(:), GRef(:), GRef_ab(:), iCutOff(:), iType(:), NZ(:), PBlock(:), Sort(:), &
                                  Sort_ab(:)
real(kind=wp), allocatable :: C(:,:), CMO(:), CMO_ab(:), ddNo(:,:), E(:), E_ab(:), MO(:), SLine(:,:), Occ(:), Occ_ab(:), Ooo(:), &
                              DOut(:), Pab(:), SphrColor(:), SphrDist(:)
character, allocatable :: cMoBlock(:)
character(len=*), parameter :: Crypt = 'fi123sd'
!---- Set size of batches
integer(kind=iwp), parameter :: nIncPack = 18*1024
integer(kind=iwp), external :: isFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
!... Prologue
nInc = nIncPack
isEner = .true.

dNorm = Zero
!ddNorm = Zero
if ((iRun == 1) .and. (levelprint >= 3)) then
  write(u6,*)
  write(u6,'(A,8I5)') 'Irreps  : ',(i,i=1,nIrrep)
  write(u6,'(A,8I5)') 'Basis   : ',(nBas(i),i=0,nIrrep-1)
  write(u6,'(A,3I5)') 'Grid Net: ',(iGridNpt(i),i=1,3)
  write(u6,*)
end if

!... Compute the size of the densities

nCMO = 0
nMOs = 0
do iIrrep=0,nIrrep-1
  nCMO = nCMO+nBas(iIrrep)**2
  nMOs = nMOs+nBas(iIrrep)
end do

if (isUHF .and. ((iDerivative /= 0) .or. isCurDens)) then
  write(u6,*) 'ERROR - Current density or derivatives not implemented for UHF!!!'
  call Abend()
end if
if (NoOrb .and. ((iDerivative /= 0) .or. isCurDens)) then
  write(u6,*) 'ERROR - Current density or derivatives not implemented with the NOORB keyword!!!'
  call Abend()
end if
if (isAtom .and. ((iDerivative /= 0) .or. isCurDens)) then
  write(u6,*) 'ERROR - Current density or derivatives not implemented with the ATOM keyword!!!'
  call Abend()
end if

call mma_allocate(CMO,nCMO,label='CMO')

call mma_allocate(E,nMOs,label='Ener')
call mma_allocate(Occ,nMOs,label='Occu')
call mma_allocate(Ooo,nMOs,label='Occ2')
if (isVirt) then
  call mma_allocate(ddNo,nMOs,nMOs,label='ddNo')
  ddNo(:,:) = Zero
end if
call mma_allocate(iType,nMOs,label='iTyp')
call mma_allocate(Sort,nMOs,label='Sort')
call mma_allocate(NZ,nMOs**2,label='Nzer')
call mma_allocate(GRef,nMOs,label='NRef')
call mma_allocate(DoIt,nMOs,label='DoIt')
call mma_allocate(Pab,nMOs,label='Pab')
GRef(:) = -1

! Read information from INPORB file

LuOrb = isFreeUnit(46)

if (isUHF) then

  ! allocate memory for extra arrays.
  call mma_allocate(CMO_ab,nCMO,label='CMO_ab')
  call mma_allocate(E_ab,nMOs,label='Ener_ab')
  call mma_allocate(Occ_ab,nMOs,label='Occ_ab')
  call mma_allocate(Sort_ab,nMOs,label='Sort_ab')
  call mma_allocate(GRef_ab,nMOs,label='NRef_ab')
  call mma_allocate(DoIt_ab,nMOs,label='DoIt_ab')

  call RdVec_(INPORB,LuOrb,'COE',1,nIrrep,NBAS,NBAS,CMO,CMO_ab,Occ,Occ_ab,E,E_ab,idum,myTitle,0,iErr,iWFtype)
  ! it can be only after SCF, so we do not need TypeIndex info
  iType(:) = 0

else !RHF case

  ! these arrays are only used for UHF, but need to be allocated anyway
  call mma_allocate(E_ab,0,label='Ener_ab')
  call mma_allocate(Occ_ab,0,label='Occ_ab')
  call mma_allocate(Sort_ab,0,label='Sort_ab')
  call mma_allocate(GRef_ab,0,label='NRef_ab')

  call RdVec(INPORB,LuOrb,'COE',nIrrep,NBAS,NBAS,CMO,Occ,E,idum,myTitle,0,iErr)

  ! construct Pab
  if (NoOrb) then
    !write(u6,*) 'nCMO,nMOs', nCMO,nMOs
    call makePab(CMO,Occ,Pab,nMOs,nMOs,nIrrep,nBas)
    !write(u6,*) 'Pab=',Pab(:)
  end if
  if (iErr == 1) then
    E(:) = Zero
  end if
  call RdVec(INPORB,LuOrb,'I',nIrrep,NBAS,NBAS,dum,dum,dum,iType,myTitle,0,iErr)
  if (iErr == 1) then
    iType(:) = 0
  end if

end if
do j=1,7
  nTypes(j) = 0
end do
do j=1,nMOs
  jj = iType(j)
  if (jj > 0) nTypes(jj) = nTypes(jj)+1
end do

close(LuOrb)

! Calculate net.

if (iRun == 1) then
  write(u6,'(A)') '   Input vectors read from INPORB'
  write(u6,'(A,A)') '   Orbital file label: ',trim(myTitle)
end if
do j=1,3
  iGridNpt(j) = iGridNpt(j)+1
end do

nCoor = iGridNpt(1)*iGridNpt(2)*iGridNpt(3)
iiCoord = nCoor

!---- if isCutOff is in used - recalculate nCoor

if (isCutOff) then
  if (isUserGrid) then
    write(u6,*) 'Not implemented'
    call Quit_OnUserError()
  end if

  call mma_allocate(iCutOff,nCoor,label='CUTFL')
  ie1 = max(iGridNpt(1)-1,1)
  ie2 = max(iGridNpt(2)-1,1)
  ie3 = max(iGridNpt(3)-1,1)
  !write(u6,*) 'vv',ie1,ie2,ie3
  iiCoord = 0
  iiiCoord = 0
  do i1=0,ie1
    do i2=0,ie2
      do i3=0,ie3
        pp(1) = GridOrigin(1)+GridAxis1(1)*i1/ie1+GridAxis2(1)*i2/ie2+GridAxis3(1)*i3/ie3
        pp(2) = GridOrigin(2)+GridAxis1(2)*i1/ie1+GridAxis2(2)*i2/ie2+GridAxis3(2)*i3/ie3
        pp(3) = GridOrigin(3)+GridAxis1(3)*i1/ie1+GridAxis2(3)*i2/ie2+GridAxis3(3)*i3/ie3
        !write(u6,'(3f8.4)') pp
        ishow = 0
        do ii=1,nAtoms
          iaia = 0
          if (abs(Coor(1,ii)-pp(1)) < CutOff) iaia = iaia+1
          if (abs(Coor(2,ii)-pp(2)) < CutOff) iaia = iaia+1
          if (abs(Coor(3,ii)-pp(3)) < CutOff) iaia = iaia+1
          if (iaia == 3) ishow = 1
        end do
        iiiCoord = iiiCoord+1
        if (ishow == 1) then
          iiCoord = iiCoord+1
          iCutOff(iiiCoord) = 1
        else
          iCutOff(iiiCoord) = 0
        end if
      end do
    end do
  end do
  write(u6,*) nCoor-iiCoord,' points are eliminated'
  !write(u6,*) 'old=',nCoor,' New=', iiCoord
  !nCoor = iiCoord
else
  call mma_allocate(iCutOff,1,label='CUTFL')
end if
if (isTheOne) nCoor = int(OneCoor(7)+0.3_wp)
if (isAtom) nCoor = nAtoms
if (isUserGrid) nCoor = nGridPoints
write(u6,*) ' Number of grid points in file:  ',nCoor
!call iXML('nPoints',nCoor)

!***********************************************************************
! And now we had to choose orbitals to draw.
!
! Sometime we had to make an automatic guess....
!***********************************************************************

call PickOrb(Nz,Sort,Gref,Sort_ab,GRef_ab,E,Occ,E_ab,Occ_ab,nShowMOs,nShowMOs_ab,isEner,nMOs,myTitle,iType)

!---- Start run over sets of grid points

!if (iBinary == 1) then
nShowMOs = nShowMOs+merge(1,0,isDensity)+merge(1,0,isSphere)+merge(1,0,isColor)
!else
!  nShowMOs=1
!end if
if (isUHF) nShowMOs_ab = nShowMOs_ab+merge(1,0,isDensity)+merge(1,0,isSphere)+merge(1,0,isColor)
nShowMOs2 = nShowMOs+nShowMOs_ab
write(u6,*)
write(u6,*) ' Total number of MOs               :',nMOs
!call iXML('nMOs',nMOs)
write(u6,*) ' Number MOs for grid               :',nShowMOs2
write(u6,*) ' Batches processed in increments of:',nInc
write(u6,*)
iPrintCount = 0

call mma_allocate(MO,nInc*nMOs,label='MOValue')
call mma_allocate(DOut,nInc,label='DOValue')

!if (isMOPack) then
!  call mma_allocate(PBlock,nInc,label='PackedBlock')
!else
call mma_allocate(PBlock,1,label='PackedBlock')
!end if
PBlock(:) = 0

!... Allocate memory for the some grid points

call mma_allocate(C,3,nInc,label='Coor')

if (isSphere) then
  call mma_allocate(SphrDist,nInc,label='SpDi')
else
  call mma_allocate(SphrDist,1,label='SpDi')
end if
if (isColor) then
  call mma_allocate(SphrColor,nInc,label='SpCo')
else
  call mma_allocate(SphrColor,1,label='SpCo')
end if

! check grids to calculate

DoIt(:) = 0
if (isUHF) DoIt_ab(:) = 0

if (.not. NoOrb) then
  do i=1,nShowMOs-merge(1,0,isDensity)-merge(1,0,isSphere)-merge(1,0,isColor)
    DoIt(GRef(i)) = 1
  end do
  if (isUHF) then
    do i=1,nShowMOs_ab-merge(1,0,isDensity)-merge(1,0,isSphere)-merge(1,0,isColor)
      DoIt_ab(GRef_ab(i)) = 1
    end do
  end if
  ifpartial = .not. isTotal
  if (isTotal) then
    do i=1,nMOs
      if (abs(Occ(i)) > Zero) then
        DoIt(i) = 1
      end if
      if (isUHF) then
        if (abs(Occ_ab(i)) > Zero) then
          DoIt_ab(i) = 1
        end if
      end if
    end do
  else
    do i=1,nMOs
      if ((Occ(i) > Zero) .and. (DoIt(i) /= 1)) then
        ifpartial = .true.
      end if
      if (isUHF) then
        if ((Occ_ab(i) > Zero) .and. (DoIt_ab(i) /= 1)) then
          ifpartial = .true.
        end if
      end if
    end do
  end if
end if

! If plotting VB orbitals : count active orbitals and make sure
! all are included in DOIT :
!
! nothing to do with UHF

iCRSIZE = 1
NBYTES = 10
NINLINE = 10
!write(u6,*) 'prepare  header '

call PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks)

!write(u6,*) 'HERE header isdone'

LuVal_ = LuVal
if (isLuscus) LuVal_ = LID
call PrintTitles(LuVal_,nShowMOs,isDensity,nMOs,GRef,isEner,Occ,iType,Crypt,NZ,E,ifpartial,isLine,isSphere,isColor,isLuscus,ncoor, &
                 nBlocks,nInc)
if (isUHF) then
  LuVal_ab_ = LuVal_ab
  if (isLuscus) LuVal_ab_ = LID_ab
  call PrintTitles(LuVal_ab_,nShowMOs_ab,isDensity,nMOs,GRef_ab,isEner,Occ_ab,iType,Crypt,NZ,E_ab,ifpartial,isLine,isSphere, &
                   isColor,isLuscus,ncoor,nBlocks,nInc)
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
!if (isAtom) then

call mma_allocate(cMoBlock,2*nIncPack,label='cMoBlock')

iiiCoord = 0
!if (isCutOff) nCoor = iiCoord
!ccccccccccccc  main loop starts here  ccccccccccccccccccccccccccccccccc
iShiftCut = 1
do iSec=1,nCoor,nInc
  mCoor = min(nInc,nCoor-iSec+1)
  !write(status,'(a,i8,a,i8)') ' batch ',iSec,' out of ',ceiling(real(nCoor,kind=wp)/real(nInc,kind=wp))
  !call StatusLine('grid_it: ',status)
  ! Generate next portion of points..
  if (isTheOne) then
    ! coords are given in specific points
    ! coords for DEBUG mode
    if (isLine) then
      ! LINE keyword
      do i=1,nCoor
        C(:,i) = OneCoor(1:3)+(OneCoor(4:6)-OneCoor(1:3))*(i-1)/(nCoor-1)
      end do
    else
      do i=1,nCoor
        C(:,i) = OneCoor(1:3)+OneCoor(4:6)*(i-1)
      end do
    end if
  else
    ! general case: we have a CUBIC box.
    do ipPO=1,mCoor
      iiiCoord = iiiCoord+1
      gv3 = real(iv3,kind=wp)/real(ive3,kind=wp)
      gv2 = real(iv2,kind=wp)/real(ive2,kind=wp)
      gv1 = real(iv1,kind=wp)/real(ive1,kind=wp)

      if (isUserGrid) then
        C(:,ipPO) = Grid(:,iSec+ipPO-1)
      else if (isCutOff) then
        ! using iCutOff
        C(:,ipPO) = 40
        if (iCutOff(iiiCoord) == 1) then
          C(:,ipPO) = GridOrigin(:)+GridAxis1(:)*gv1+GridAxis2(:)*gv2+GridAxis3(:)*gv3
        end if
      else
        C(:,ipPO) = GridOrigin(:)+GridAxis1(:)*gv1+GridAxis2(:)*gv2+GridAxis3(:)*gv3
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
      !write(u6,*) 'ipo',ipP0
    end do
    ! end of CUBIC box
  end if
  ! end of coordinates.
  !VV: FIXME; separate color and sphere
  !if (isSphere .and. isColor) then
  !  call Sphr_Grid(Coor,mCoor,C,SphrDist,SphrColor)
  !end if

  if (NoOrb) then
    nDrv = 0
    call MOEval(MO,nMOs,mCoor,C,Pab,nMOs,DoIt,nDrv,1)
  else
    nDrv = 0
    call MOEval(MO,nMOs,mCoor,C,CMO,nCMO,DoIt,nDrv,1)
  end if

  !... Write out values

  nLine = min(nShowMOs,20)
  if (.not. isLine) nLine = 0
  call mma_allocate(SLine,nLine,mCoor,label='Line')
  if (levelprint < 2) iPrintCount = 100
  !VV BUG Update iCutOff

  call DumpM2Msi(iRun,Luval_,LID,nShowMOs,isDensity,nMOs,GRef,Occ,MO,DOut,mCoor,iGauss,nInc,isMOPack,PBlock,cMoBlock, &
                 nBytesPackedVal,dnorm,Crypt,isTheOne,isLine,iBinary,isEner,iType,NZ,E,SLine,nLine,C,iPrintCount,isDebug,isCutOff, &
                 iCutOff(iShiftCut),isSphere,SphrDist,isColor,SphrColor,isLuscus,NBYTES,NINLINE)
  !if (isXField == 1) call dcopy_(mCoor,DOut,1,DOutXF,1)
  if (isUHF) then
    !VV:
    nDrv = 0
    call MOEval(MO,nMOs,mCoor,C,CMO_ab,nCMO,DoIt_ab,nDrv,1)

    !... Write out values

    call DumpM2Msi(iRun,Luval_ab_,LID_ab,nShowMOs_ab,isDensity,nMOs,GRef_ab,Occ_ab,MO,DOut,mCoor,iGauss,nInc,isMOPack, &
                   PBlock,cMoBlock,nBytesPackedVal,dnorm,Crypt,isTheOne,isLine,iBinary,isEner,iType,NZ,E_ab,SLine,nLine,C, &
                   iPrintCount,isDebug,isCutOff,iCutOff(iShiftCut),isSphere,SphrDist,isColor,SphrColor,isLuscus,NBYTES,NINLINE)

  end if

  call mma_deallocate(SLine)
  if (isCutOff) iShiftCut = iShiftCut+mCoor
  if (isVirt) then
    do iiMO=1,nMOs
      do jjMO=1,nMOs
        if ((Occ(iiMO) > 1.1_wp) .and. (Occ(jjMO) < 0.9_wp)) then
          !  here if this is a pair Occ-Virt

          do i=1,nMOs
            Ooo(i) = Zero
            if (i == iiMO) Ooo(i) = Two-Virt
            if (i == jjMO) Ooo(i) = Virt
          end do

          call outmo(0,2,MO,Ooo,DOut,mCoor,nMOs)
          dd = Zero
          do j=1,mCoor
            dd = dd+DOut(j)
          end do
          !ddNorm = ddNorm+dd
          call save_ddNorm(dd,iiMO,jjMO,ddNo,nMOs)
        end if
      end do
    end do
  end if
end do !iSec
!ccccccccccccc  main loop ends here  ccccccccccccccccccccccccccccccccccc

call mma_deallocate(cMoBlock)

write(u6,*)
! Check norms

if (.not. isAtom) then

  det3 = GridAxis1(1)*(GridAxis2(2)*GridAxis3(3)-GridAxis2(3)*GridAxis3(2))- &
         GridAxis1(2)*(GridAxis2(1)*GridAxis3(3)-GridAxis2(3)*GridAxis3(1))+ &
         GridAxis1(3)*(GridAxis2(1)*GridAxis3(2)-GridAxis2(2)*GridAxis3(1))
  det3 = det3/((iGridNpt(1)-1)*(iGridNpt(2)-1)*(iGridNpt(3)-1))
  dNorm = dNorm*det3

  !write(u6,*) 'dNorm=',dNorm

  if (isVirt) then
    call print_ddNorm(nMOs,ddNo,det3)

    !write(u6,*) 'ddNorm=',ddNorm*det3

  end if

  call Add_Info('GRIDIT_NORM',[dNorm],1,6)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!call DAClos(LuVal)
if (isLuscus) then
  LINE = ' '
  call PRINTLINE(LUVAL_,LINE,1,.false.)
  LINE = ' </DENSITY>'
  call PRINTLINE(LUVAL_,LINE,11,.false.)
  LINE = ' <BASIS>'
  call PRINTLINE(LUVAL_,LINE,8,.false.)
  LINE = ' </BASIS>'
  call PRINTLINE(LUVAL_,LINE,9,.false.)
  LINE = ' <INPORB>'
  call PRINTLINE(LUVAL_,LINE,9,.false.)
end if

if (.not. isLine) then
  write(str,'(9i8)') nIrrep,(nBas(i),i=0,nIrrep-1)
  call PrintLine(LuVal_,STR,72,.false.)
  !if (iBinary == 0) then
  !  write(LuVal_,'(a)') str
  !else
  !  write(LuVal_) str
  !end if
  ! Well, to avoid rewritting of Cerius2 we use old INPORB format temporary!
  LuOrb = isFreeUnit(46)
  call molcas_open_ext2(luorb,INPORB,'sequential','formatted',istatus,.false.,irecl,'old',is_error)
  !do
  !  read(LuOrb,'(a)',iostat=istatus) str
  !  if (istatus /= 0) exit
  !  if (str == '#ORB') exit
  !end do
  do
    read(LuOrb,'(a)',iostat=istatus) str
    if (istatus /= 0) exit
    !if (str(1:1) == '#') cycle
    iLen = len(str)
    do i=iLen,1,-1
      if (str(i:i) /= ' ') exit
    end do
    call PrintLine(LuVal_,STR,I,.false.)
    !if (iBinary == 0) then
    !  write(LuVal_,'(a)') str(1:i)
    !else
    !  write(LuVal_) str(1:i)
    !end if
  end do
  close(LuOrb)
end if ! isLine
close(LuVal)
if (isUHF) close(LuVal_ab)
if (isTheOne) then
  MM = mCoor-1
  if (MM > 10) MM = 10
  call Add_Info('GRIDIT_ONE',DOut,MM,6)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!... Epilogue, end

!end if ! isAtom

if (isAtom) then
  mCoor = nCoor
  nDrv = 0
  call MOEval(MO,nMOs,mCoor,Coor,CMO,nCMO,DoIt,nDrv,1)
  call outmo(0,2,MO,Occ,DOut,nCoor,nMOs)
  write(u6,'(a)') repeat('*',60)
  if (ifpartial) then
    write(u6,'(a5,3a10,a20)') 'Atom','x','y','z','Density (partial)'
  else
    write(u6,'(a5,3a10,a20)') 'Atom','x','y','z','Density'
  end if
  do i=1,nAtoms
    write(u6,'(a5,3f10.3,es20.10)') AtomLbl(i),Coor(:,i),DOut(i)
  end do
  call Add_Info('GRIDIT_ATOM',DOut(1:nAtoms),nAtoms,6)

end if

if (isLuscus) then
  call PRTLUSENDGRID(LID)
  if (isUHF) call PRTLUSENDGRID(LID_ab)
end if

call mma_deallocate(C)
call mma_deallocate(SphrDist)
call mma_deallocate(SphrColor)

call mma_deallocate(MO)
call mma_deallocate(DOut)

if (isUHF) then
  call mma_deallocate(CMO_ab)
  call mma_deallocate(DoIt_ab)
end if
call mma_deallocate(E_ab)
call mma_deallocate(Occ_ab)
call mma_deallocate(Sort_ab)
call mma_deallocate(GRef_ab)
if (isUserGrid) call mma_deallocate(Grid)
call mma_deallocate(PBlock)

call mma_deallocate(Pab)
call mma_deallocate(DoIt)

if (isVirt) call mma_deallocate(ddNo)
call mma_deallocate(iType)
call mma_deallocate(Sort)
call mma_deallocate(NZ)
call mma_deallocate(GRef)

call mma_deallocate(E)
call mma_deallocate(Occ)
call mma_deallocate(Ooo)

call mma_deallocate(CMO)

call mma_deallocate(AtomLbl)
call mma_deallocate(Coor)

call mma_deallocate(iCutOff)

return

end subroutine DrvMO
