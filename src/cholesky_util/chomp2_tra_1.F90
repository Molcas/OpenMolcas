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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Tra_1(COcc,CVir,Diag,DoDiag,Wrk,lWrk,iSym)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: transform Cholesky vectors to (ai) MO basis for symmetry
!          block iSym. Files are assumed open.
!          If requested (DoDiag=.true.), compute (ai|ai) integral
!          diagonal.

use Cholesky, only: InfVec, nnBstR, NumCho
use ChoMP2, only: lUnit_F, nT1am, nT1AOT
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: COcc(*), CVir(*)
real(kind=wp), intent(inout) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
integer(kind=iwp), intent(in) :: lWrk, iSym
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: iAdr, iBat, iLoc, iOpt, irc, iRed, iRedC, iVec, iVec1, iVec2, jNum, jVec, jVec1, kChoAO, kChoMO, kEnd0, &
                     kHlfTr, kOff, kOffMO, lChoAO, lChoMO, lHlfTr, lRead, lWrk0, lWrk1, mUsed, nMOVec, NumBat, NumV
character(len=*), parameter :: SecNam = 'ChoMP2_Tra_1'
integer(kind=iwp), external :: Cho_lRead

! Check if anything to do.
! ------------------------

if ((NumCho(iSym) < 1) .or. (nT1am(iSym) < 1)) return

! Initialize Diag (if needed).
! ----------------------------

if (DoDiag) Diag(1:nT1am(iSym)) = Zero

! Allocate memory for half-transformed vector.
! --------------------------------------------

lHlfTr = nT1AOT(iSym)

kHlfTr = 1
kEnd0 = kHlfTr+lHlfTr
lWrk0 = lWrk-kEnd0+1
if (lWrk0 < (nT1am(iSym)+nnBstR(iSym,1))) call SysAbendMsg(SecNam,'insufficient memory','[0]')

! Reserve memory for reading AO vectors.
! --------------------------------------

lRead = Cho_lRead(iSym,lWrk0)
if (lRead < 1) then
  write(u6,*) SecNam,': memory error: lRead = ',lRead
  call SysAbendMsg(SecNam,'memory error',' ')
  lWrk1 = 0 ! to avoid compiler warnings...
else
  lWrk1 = lWrk0-lRead
  if (lWrk1 < nT1am(iSym)) then
    lWrk1 = nT1am(iSym)
    lRead = lWrk-nT1am(iSym)
  end if
end if

! Set up batch.
! -------------

nMOVec = min(lWrk1/nT1am(iSym),NumCho(iSym))
if (nMOVec < 1) call SysAbendMsg(SecNam,'insufficient memory','[1]')
NumBat = (NumCho(iSym)-1)/nMOVec+1

! Set reduced set handles.
! ------------------------

iRedC = -1
iLoc = 3

! Transform each batch of vectors and compute diagonal contributions
! (if requested).
! ------------------------------------------------------------------

do iBat=1,NumBat

  if (iBat == NumBat) then
    NumV = NumCho(iSym)-nMOVec*(NumBat-1)
  else
    NumV = nMOVec
  end if
  iVec1 = nMOVec*(iBat-1)+1
  iVec2 = iVec1+NumV-1

  lChoMO = nT1am(iSym)*NumV

  kChoMO = kEnd0
  kChoAO = kChoMO+lChoMO
  lChoAO = lWrk0-kChoAO+1

  kOffMO = kChoMO
  jVec1 = iVec1
  do while (jVec1 <= iVec2)

    jNum = 0
    call Cho_VecRd(Wrk(kChoAO),lChoAO,jVec1,iVec2,iSym,jNum,iRedC,mUsed)
    if (jNum < 1) call SysAbendMsg(SecNam,'insufficient memory','[2]')

    kOff = kChoAO
    do jVec=1,jNum
      iVec = jVec1+jVec-1
      iRed = InfVec(iVec,2,iSym)
      if (iRedC /= iRed) then
        irc = 0
        call Cho_X_SetRed(irc,iLoc,iRed)
        if (irc /= 0) call SysAbendMsg(SecNam,'error in Cho_X_SetRed',' ')
        iRedC = iRed
      end if
      call ChoMP2_TraVec(Wrk(kOff),Wrk(kOffMO),COcc,CVir,Wrk(kHlfTr),lHlfTr,iSym,1,1,iLoc)
      kOff = kOff+nnBstR(iSym,iLoc)
      kOffMO = kOffMO+nT1am(iSym)
    end do

    jVec1 = jVec1+jNum

  end do

  iOpt = 1
  iAdr = nT1am(iSym)*(iVec1-1)+1
  call ddaFile(lUnit_F(iSym,1),iOpt,Wrk(kChoMO),lChoMO,iAdr)

  if (DoDiag) then
    do iVec=1,NumV
      kOff = kChoMO+nT1am(iSym)*(iVec-1)-1
      Diag(1:nT1am(iSym)) = Diag(1:nT1am(iSym))+Wrk(kOff+1:kOff+nT1am(iSym))**2
    end do
  end if

end do

end subroutine ChoMP2_Tra_1
